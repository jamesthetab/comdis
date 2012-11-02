# HybridGPool() Function
# Last updated: 29 Oct. 2012 by JRM

#Description:
#Generating a global species pool using a hybrid approach of Roche et al. 2012 and 
#Dobson 2004. 


HybridGPool <- function(Y=10, modalO=3, z=0.1, a=2, b=1, m=1.5){
  
  #Start by generating a community based on Preston's law
  from <- 1
  to <- 15
  P <- c(from:to) #Sequence of octaves
  
  S <- Y*exp(-(z*(P-modalO)^2)) #Modify the modal octave with "mode"
  S <- round(S) #Round "S" to get whole number species
  Plog2 <- 2^(P) #Octaves in terms of individuals
  Preston <- data.frame(S, P, Plog2)
  
  Abund <- NULL
  Rank <- NULL
  for(i in 1:length(S)){
    Abund <- c(Abund, rep(Plog2[i], S[i]))
    Rank <- c(Rank, rep(P[i], S[i]))
  }
  #Prepare masses, based on ranks (octaves)
  logM <- a - b*log(Rank)
  M <- 10^logM
  
  #Assemble Species Traits:
  species <- mat.or.vec(length(Abund), 10) # initialize matrix
  species[,1] <- c(1:length(Abund)) # assign species ID's
  colnames(species) <- c("ID","b","d","DD","K","v","g","Bii","weight","R0ii")
  
  #Create a truncated gamma and sample R0ii. Ordered so smaller species have larger R0ii
  library(distr)
  k <- 0.3; w <- 3 #Shape and scale parameters of the truncated gamma. Can alter if necessary.
  TG <- Truncate(Gammad(scale=k, shape=w), lower=0, upper=2)
  R0ii <- r(TG)(length(Abund))
  R0ii <- R0ii[order(R0ii)]
  #Follow De Leo & Dobson 1996 & some Roche et al. 2012
    for (i in 1:length(Abund)) {
      #Carrying Capacity:
      species[i, 5] <- Abund[i]
      species[i, 4] <- 0.01 # density dependence assumed equal across species
      # Weights from Preston's law (octaves):
      species[i, 9] <- M[i]
      # Birth rate = r + mu = 0.6w^-0.27 + 0.4w^-0.26 
      species[i, 2] <- 0.6 * species[i, 9] ^ -0.27 + 0.4 * species[i, 9] ^ -0.26
      # Death rate = 0.4w^-0.26
      species[i, 3] <- 0.4 * species[i, 9] ^ -0.26  
      species[i, 6] <- (m - 1) * species[i, 3] # virulence: alpha=(m-1)*mu 
      species[i, 7] <- species[i, 3] # recovery rate assumed to scale with body size as mu
      # Intraspecific transmission rate Bii=0.0247*Ro(m+(recovery rate/mu))*w^0.44
      species[i, 8] <- R0ii[i] * 0.0247 * (m + species[i, 7] / species[i, 3]) * species[i, 9]^0.44 
      species[i, 10] <- R0ii[i] # store R0ii for each species
    }
  
  
  stopifnot (species[, 8] >= 0 & species[, 8] <= 1) # Bii must be between 0 and 1
  stopifnot (species[, 2] > 0) # birth rate has to be >0
  stopifnot (species[, 3] > 0) # death rate has to be >0
  stopifnot (species[, 6] >= 0) # pathogen-induced mortality can't be negative (pathogens have to reduce survival)
  stopifnot (species[, 7] >= 0) # recovery rate can't be negative
  
  list(Nglobal = length(Abund),
       Y = Y,
       z = z,
       modalO = modalO,
       a = a,
       b = b,
       m = m,
       global.pool = data.frame(species))
}

##### EXAMPLE #####
poolA<-HybridGPool()
run<-ComSim(poolA, iter=1000, comsizes = c(2:poolA$Nglobal), mode="freq", Kmeth="fixed",
            SampleMass=T)
plot(run, method="raw")
plot(run, method="adj")
plot(run, "even")
plot(run, "rich")