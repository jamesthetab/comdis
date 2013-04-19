# Recreating the individual based model of Allouche et al. 2012
# with discrete time steps

## BACKGROUND
# This model is spatially implicit, with A "sites" that are all equally
# connected. Each site falls on an environmental condition axis, receiving
# some value "E" that characterizes local conditions. The amount of 
# habitat heterogeneity among sites is represented by Emin-Emax, with
# E ~ Unif(Emin, Emax). The local range of environmental conditions is a
# subset of the global range (ERmin, ERmax). There are N species in the
# regional pool that can colonize the habitat patches. Each species has 
# some environmental optimum mu.i, and some niche width sigma.i, which together
# define a Gaussian function that defines the probability of establishment 
# given some environmental condition (see plot below or Allouche supplement).
# Other than different environmental optima that relate to the probability of 
# establishment, all species are considered equivalent, having the same 
# Pr(death), Pr(reproduction), Pr(attempt to colonize from regional pool).

# The individual based model uses discrete timesteps. Each habitat patch/site
# can have only one individual at a time. At each timestep, individuals that 
# occupy a habitat patch die with some probability, reproduce with some probability
# (number of offspring is currently Poisson distributed), or do nothing. Offspring
# colonize empty sites. If there are fewer empty sites than there are offspring, 
# the successful colonist is randomly selected from the pool of attempting colonists. 
# Individuals from the regional pool attempt to colonize each empty habitat patch at
# every timestep, with some low probability of attempting to colonize. If an attempt
# is made, then the individual successfully colonizes with Pr(establishment), a function
# of local environmental conditions. 

## PARAMETERS
# A:        Number of sites
# ERmin:    Global environmental minimum
# ERmax:    Global environmental maximum
# Emin:     Local environmental minimum
# Emax:     Local environmental maximum
# N:        Number of species
# sigma.i:  Niche width of each species
# pM:       Pr(mortality) at each timestep
# pR:       Pr(reproducing) at each timestep
# R:        Mean # offspring produced w/ reproduction (# babies|repro ~ pois(R))
# I:        Per-patch Pr(attempting to colonize from regional pool)

#-------------#
# Example run #
#-------------#
source("hostIBM.R")

out <- hostIBM(pM=.1, Emin=-20, Emax=20, sig=1, R=3, N=100, A=300, I=.2)
plot(out$p.occ, type="l", xlab="Timestep", ylab="Proportion of patches occupied")
plot(out$richness, type="l", xlab="Timestep", ylab="Species richness")

# overlay more trajectories
for (i in 1:3){
  out <- hostIBM()
  lines(x=1:length(out$richness), y=out$richness, col="blue")
}

require(ggplot2)
ggplot(out$niche.d, aes(x=E, y=Pr.estab, group=species)) + geom_line() +
  geom_point(aes(x=E, y=-0.005), alpha=.2)

# Plot the abundance of each species over time
timesteps <- length(out$richness)
plot(x=NULL, y=NULL, xlim=c(0, timesteps), ylim=c(0, 40), type="n",
     xlab="Time", ylab="Abundance")
for (i in 1:N){
  # sum abundance across all sites
  abundance <- apply(out$state[,,i], MARGIN=1, sum)
  lines(x=1:timesteps, y=abundance)
  text(x=timesteps*1.02, y=abundance[timesteps], labels=as.character(i))
}

# disperse born individuals to all patches randomly, not just unoccupied patches
# allow for displacement of parasites by other parasites
# include commensal dynamics & virulence

#----------------------------------------------------------------------#
# Investigate the effect of resource heterogeneity on species richness #
#----------------------------------------------------------------------#
# 1. generate random intervals that are subsets of the global interval
global.median <- median(c(ERmin, ERmax))
n.intervals <- 10
lower.limits <- seq(ERmin, global.median-.5, length.out=n.intervals)
upper.limits <- seq(ERmax, global.median, length.out=n.intervals)
cbind(lower.limits, upper.limits)
hab.het <- upper.limits - lower.limits # habitat heterogeneity width
hist(hab.het)
n.iter <- 1 # number of iterations per interval
# We want a dataframe with one row per interval and columns:
# mean(richness|stability), interval width, variance(richness|stability)
mean.rich <- rep(NA, n.intervals)
var.rich <- rep(NA, n.intervals)
end.rich <- array(NA, dim=c(n.intervals, n.iter))
for (i in 1:n.intervals){
  for (iter in 1:n.iter){
    Emin <- lower.limits[i]
    Emax <- upper.limits[i]
    out <- hostIBM(pM=.8, Emin=Emin, Emax=Emax, sig=3, pR=1, R=1, N=50, A=50, I=.2, timesteps=200)
    timesteps <- length(out$richness)
    end.rich[i, iter] <- mean(out$richness[timesteps-100:timesteps])
  }
  print(i)
}
# put data into long form
end.richness <- c(end.rich)
interval <- rep(hab.het, each=n.iter) # width of interval
plot(interval, end.rich)
