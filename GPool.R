###########
# GPool() #
###########
# Last edited Oct 26, 2012 MJ
### Purpose ###
# Produce a global species pool, drawing values for demographic and
# transmission parameters for some number of species.

### Input ###
# Nglobal: number of species in the global pool

# globmeth: the way that species traits are assigned
#  "random" - traits are drawn from a uniform distribution
#  "allom" - traits are produced using allometry (Dobson 2004 & ref. therein)
#  "allomr" - first produce a range of values based on allometry, then 
#             assign each trait randomly among species
# Bmeth: how to assign intraspecific transmission terms (Bii)
#  "allomr" - random assignment of traits among species, 
#             using allometry to generate values
#  "beta" - randomly draw from a beta distribution
#  "random" - randomly draw from uniform distribution
#  "scalewd" - scale transmission with death rates
#  "fixed" - user specifies a value to assign to all species
# Bii: if Bmeth = "fixed", the value of Bii given to all species
# a: allometric scaling parameter, influences the range of weights produced
# m: pathogen virulence. More virulent pathogens have larger values of m. 
#   Infected individuals suffer an m-fold reduction of life expectancy.
# R0ii: Intraspecific R0, used in calculating Bii values with allometry
# min.weight: weight of the smallest species in the global pool

### Output ###
# A list of control parameters (Nglobal, Bmeth, Bii, globmeth, a, m, R0ii)
# and species traits (global.pool)

GPool <- function(Nglobal=20, Bmeth="allomr", Bii=0.3, globmeth="allom",
                  a=1.15, m=1.5, R0ii=1.5, min.weight=1){
  species <- mat.or.vec(Nglobal, 10) # initialize matrix
  species[,1] <- c(1:Nglobal) # assign species ID's
  colnames(species) <- c("ID","b","d","DD","K","v","g","Bii","weight","R0ii")
  
  if (globmeth == "random"){
    for(i in 1:Nglobal){
      species[i,2] <- runif(1,min=0,max=1) #Randomly assign birth rate
      species[i,6] <- runif(1,min=0,max=1) #Randomly assign death rate due to pathogen
      species[i,7] <- runif(1,min=0,max=1) # Randomly assign recovery rate
    }
    
    for (i in 1:Nglobal){ # assign death rates based on birth rates
      if (species[i,2] <= .5){
        species[i,3] <- runif(1, min=0, max=species[i,2])
      } else {
        species[i,3] <- runif(1, min=0.5, max=species[i,2])
      }
      species[i,4] <- 0.01 # assign strength of density dependence (DD = 0.1)
      # finally assign carrying capacities according to K=(b-d)/DD
      species[i,5] <- ceiling((species[i,2] - species[i,3])/species[i,4]) 
    }
    
    if(Bmeth == "random"){
      for (i in 1:Nglobal){
        species[i,8] <- runif(1, min=0, max=1)
      }
    }
    if (Bmeth == "beta"){
      for (i in 1:Nglobal){
        species[i,8] <- rbeta(1, 1, 5)
      }
    } #Can alter shape & scale parameters if necessary
    if (Bmeth == "scalewd"){
      for (i in 1:Nglobal){
        species[i,8] <- species[i, 3]
      }
    } # Intraspecific transmission is equal to death rate - long-lived species transmit less
    if (Bmeth == "fixed"){
      species[,8] <- Bii
    }
    
    # calculate intraspecific Ro as (K*Beta)/(d+recovery+virulence) (De Leo and Dobson 1996)
    species[,10] <- (species[,5] * species[,8]) / (species[,6] + species[,3] + species[,7]) 
  }
  
  if (globmeth == "allom") { # follows De Leo & Dobson 1996
    for (i in 1:Nglobal) {
      # Weight (kg) = ws*a^i; i is species ID
      species[i, 9] <- min.weight * a ^ (i - 1) 
      # Birth rate = r + mu = 0.6w^-0.27 + 0.4w^-0.26 
      species[i, 2] <- 0.6 * species[i, 9] ^ -0.27 + 0.4 * species[i, 9] ^ -0.26
      # Death rate = 0.4w^-0.26
      species[i, 3] <- 0.4 * species[i, 9] ^ -0.26 
      species[i, 4] <- 0.01 # density dependence assumed equal across species
      species[i, 5] <- ((species[i, 2] - species[i, 3]) / species[i, 4])  # K=(b-d)/DD
      species[i, 6] <- (m - 1) * species[i, 3] # death rate from infection: alpha=(m-1)*mu 
      species[i, 7] <- species[i, 3] # recovery rate assumed to scale with body size as mu
      # Intraspecific transmission rate Bii=0.0247*Ro(m+(recovery rate/mu))*w^0.44
      species[i, 8] <- R0ii * 0.0247 * (m + species[i, 7] / species[i, 3]) * species[i, 9]^0.44 
      species[i, 10] <- R0ii # store R0ii for each species
    }
  }
  
  if (globmeth == "allomr"){
    weights <- min.weight * a ^ (c(1:Nglobal) - 1)
    birth.rates <- 0.6 * weights ^ -0.27 + 0.4 * weights ^ -0.26
    death.rates <- 0.4 * weights ^ -0.26
    recov.rates <- death.rates
    infect.death <- (m - 1) * death.rates
    intra.transmission <- R0ii * 0.0247 * (m + recov.rates / death.rates) * weights ^ 0.44
    dens.dep <- 0.01
    bd.pair <- cbind(birth.rates, death.rates) # pair birth and death rates
    species[, 2:3] <- bd.pair[sample(1:nrow(bd.pair), 
                                     size = nrow(bd.pair), replace = F), ] # randomize order of birth and death rate pairs
    species[,6] <- sample(infect.death, replace=F, size = Nglobal) # assign infect. induced mortality
    species[,7] <- sample(recov.rates, replace = F, size = Nglobal) # assign recovery rates
    species[,4] <- dens.dep # assume  strength of density dependence is equal across species
    species[,5] <- ceiling((species[, 2] - species[, 3]) / species[, 4]) # K=(b-d)/DD
    species[,9] <- sample(weights, replace = F, size = Nglobal)
    
    # assign Bii for each species
    if (Bmeth == "random") {
      species[,8]<-runif(Nglobal,min=0,max=1)
    }
    if(Bmeth=="allomr") {
      species[,8]<-sample(intra.transmission,replace=F,size=Nglobal)
    }
    if(Bmeth=="beta") {
      species[,8]<-rbeta(Nglobal,1,5) #can alter shape & scale parameters
    }
    if(Bmeth=="scalewd") {
      species[,8]<-species[,3] # Bii = death rate, long-lived species transmit less
    }
    if(Bmeth=="fixed") {
      species[,8]<-Bii
    }
    
    # calculate intraspecific Ro as (K*Beta)/(d+recovery+virulence) (De Leo and Dobson 1996)
    species[,10] <- (species[,5] * species[,8]) / (species[,6] + species[,3] + species[,7])
  }
  
  stopifnot (species[, 8] >= 0 & species[, 8] <= 1) # Bii must be between 0 and 1
  stopifnot (species[, 2] > 0) # birth rate has to be >0
  stopifnot (species[, 3] > 0) # death rate has to be >0
  stopifnot (species[, 6] >= 0) # pathogen-induced mortality can't be negative (pathogens have to reduce survival)
  stopifnot (species[, 7] >= 0) # recovery rate can't be negative
  
  list(Nglobal = Nglobal,
       Bmeth = Bmeth,
       Bii = Bii,
       globmeth = globmeth,
       a = a,
       m = m,
       R0ii = R0ii,
       global.pool = data.frame(species))
}
