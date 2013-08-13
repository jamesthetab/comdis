###########
# GPool() #
###########
# Last edited Nov 6, 2012 MJ deleted globmeth="random", Bmeth="scalewd",
# included scaling parameter for recovery rate (kRecovery), and fixed 
# calculation of Bii's so that each species gets the same value of R0ii, 
# eliminated the for-loop in globmeth="allom" with vectorization, and 
# included option to draw R0ii values from truncated gamma distribution
# as in HybridGPool()

### Purpose ###
# Produce a global species pool, drawing values for demographic and
# transmission parameters for some number of species.

### Input ###
# Nglobal: number of species in the global pool

# globmeth: the way that species traits are assigned
#  "allom" - traits are produced using allometry (Dobson 2004 & ref. therein)
#  "allomr" - first produce a range of values based on allometry, then 
#             assign each trait randomly among species
# Bmeth: how to assign intraspecific transmission terms (Bii)
#  "allomr" - random assignment of traits among species, 
#             using allometry to generate values
#  "beta" - randomly draw from a beta distribution
#  "random" - randomly draw from uniform distribution
#  "fixed" - user specifies a value to assign to all species
# Bii: if Bmeth = "fixed", the value of Bii given to all species
# a: allometric scaling parameter, influences the range of weights produced
# m: pathogen virulence. More virulent pathogens have larger values of m. 
#   Infected individuals suffer an m-fold reduction of life expectancy.
# R0ii: Intraspecific R0, used in calculating Bii values with allometry
#   can be set to one number, applied across all species, or set to "TG", where
#   R0ii values are drawn from a truncated gamma distribution. By default, 
#   small species receive highest values of R0ii if globmeth="allom"
# min.weight: weight of the smallest species in the global pool
# kRecovery: determines the scaling of infectious period relative to average
#   lifespan. If the avg infectious period is 1/10th the average lifespan,
#   kRecovery=10 (it's the denominator)

### Output ###
# A list of control parameters (Nglobal, Bmeth, Bii, globmeth, a, m, R0ii)
# and species traits (global.pool)

GPool <- function(Nglobal=20, Bmeth="allomr", Bii=0.3, globmeth="allom",
                  a=1.15, m=1.5, R0ii="TG", min.weight=1, kRecovery=10, seed=runif(1, 1, 10000),
                  k = 0.5, w = 5) {
	stopifnot(R0ii == "TG" | (class(R0ii) == "numeric" & length(R0ii) == 1))
  species <- mat.or.vec(Nglobal, 10) # initialize matrix
  species[,1] <- c(1:Nglobal) # assign species ID's
  colnames(species) <- c("ID","b","d","DD","K","v","g","Bii","weight","R0ii")
  
  if (R0ii == "TG"){
  	library(distr)
    #Shape and scale parameters of the truncated gamma. Can alter if necessary.
  	TG <- Truncate(Gammad(scale=k, shape=w), lower=0, upper=10)
    set.seed(seed)
  	R0ii <- r(TG)(Nglobal)
  	R0ii <- R0ii[rev(order(R0ii))]
  }
  
  if (globmeth == "allom") { # follows De Leo & Dobson 1996
    species[, 9] <- min.weight * a ^ (species[, 1] - 1)
    # Birth rate = r + mu = 0.6w^-0.27 + 0.4w^-0.26 
    species[, 2] <- 0.6 * species[, 9] ^ -0.27 + 0.4 * species[, 9] ^ -0.26
    # Death rate = 0.4w^-0.26
    species[, 3] <- 0.4 * species[, 9] ^ -0.26 
    species[, 4] <- 0.01 # density dependence assumed equal across species
    species[, 5] <- ((species[, 2] - species[, 3]) / species[, 4])  # K=(b-d)/DD
    species[, 6] <- (m - 1) * species[, 3] # death rate from infection: alpha=(m-1)*mu 
    species[, 7] <- kRecovery*species[, 3] # recovery rate = k*d 
    # Intraspecific transmission rate Bii=Ro*(d+v+g)/K
    species[, 8] <-  (R0ii * (species[, 3] + species[, 6] + species[, 7])) / species[, 5] # direct calculation from R0
    species[, 10] <- (species[, 5] * species[, 8]) / (species[, 6] + species[, 3] + species[, 7]) # R0ii
  }
  
  if (globmeth == "allomr"){
    weights <- min.weight * a ^ (c(1:Nglobal) - 1)
    birth.rates <- 0.6 * weights ^ -0.27 + 0.4 * weights ^ -0.26
    death.rates <- 0.4 * weights ^ -0.26
    recov.rates <- kRecovery*death.rates
    infect.death <- (m - 1) * death.rates
    dens.dep <- 0.01
    bd.pair <- cbind(birth.rates, death.rates) # pair birth and death rates
    species[, 2:3] <- bd.pair[sample(1:nrow(bd.pair), 
                                     size = nrow(bd.pair), replace = F), ] # randomize order of birth and death rate pairs
    species[, 6] <- sample(infect.death, replace=F, size = Nglobal) # assign infect. induced mortality
    species[, 7] <- sample(recov.rates, replace = F, size = Nglobal) # assign recovery rates
    species[, 4] <- dens.dep # assume  strength of density dependence is equal across species
    species[, 5] <- ceiling((species[, 2] - species[, 3]) / species[, 4]) # K=(b-d)/DD
    species[, 9] <- sample(weights, replace = F, size = Nglobal)
    intra.transmission <- R0ii * (death.rates + infect.death + recov.rates) / species [, 5]
    
    # assign Bii for each species
    if (Bmeth == "random") {
      species[, 8]<-runif(Nglobal,min=0,max=1)
    }
    if(Bmeth=="allomr") {
      species[, 8]<-sample(intra.transmission,replace=F,size=Nglobal)
    }
    if(Bmeth=="beta") {
      species[, 8]<-rbeta(Nglobal,1,5) #can alter shape & scale parameters
    }
    if(Bmeth=="fixed") {
      species[, 8]<-Bii
    }
    
    # calculate intraspecific Ro as (K*Beta)/(d+recovery+virulence) (De Leo and Dobson 1996)
    species[, 10] <- (species[, 5] * species[, 8]) / (species[, 6] + species[, 3] + species[, 7])
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

# Examples
pool.tg <- GPool(globmeth="allom", R0ii="TG", seed=2, k=.1, w=10)
pool.fixed <- GPool(globmeth="allom", R0ii=1)

pairs(pool.tg$global.pool)

ggplot(pool.tg$global.pool) + geom_line(aes(x=weight, y=K)) + 
  theme_bw() + xlab("Weight") + ylab("Carrying capacity")

ggplot(pool.tg$global.pool, aes(x=weight, y=R0ii)) + geom_point() + stat_smooth(color="darkgreen") +
  theme_bw() + xlab("Weight") + ylab(expression(paste("Host competence (", R[0],")", sep="")))