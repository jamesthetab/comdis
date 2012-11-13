################
### RegSim() ###
################ 

### Purpose ###
# Simulate communities, varying both regional and local species richness.
# This differs from ComSim in that species lost from the regional pool are
# done so in succession, either randomly or based on carrying capacity.

### Input ###
# regmax: maximum regional species richness
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
# risk: "K" means species with lowest K is lost from pool, "random" is 
#   random sequential loss of species from pool
# m: pathogen virulence. More virulent pathogens have larger values of m. 
#   Infected individuals suffer an m-fold reduction of life expectancy.
# R0ii: Intraspecific R0, used in calculating Bii values with allometry
# min.weight: weight of the smallest species in the global pool
# iter: number of simulated local communities for each regional pool
# comsizes: range of species richnesses for local communities
# Kmeth: Either "free" (total community carrying capacity is free to vary
#  as species are added), or "fixed" (community K is fixed, and abundances are
#  adjusted to sum to community K, while relative abundances do not change).

### Output ###
# A matrix with information from each simulated community.

RegSim <- function(mode, regmax = 10, globmeth = "allomr", Bmeth = "fixed",
                   Bii = 0.3, a = 1.15, risk = "random", m = 1.5, R0ii = 1.5,
                   min.weight = 1, iter = 10, comsizes = NULL, Kmeth = "free") {
  regdat <- NULL
  pool <- NULL
  if (risk == "random") {
    for (i in 1:(regmax - 1)) {
      if (i == 1) {
        pool <- GPool(Nglobal = regmax, globmeth = globmeth, 
                      Bmeth = Bmeth, Bii = Bii, a = a, m = m, 
                      R0ii = R0ii, min.weight = min.weight)
      } else {
        del <- sample(pool$global.pool[, 1], size = 1) # select species to delete
        pool$global.pool <- pool$global.pool[-del, ]
      }
      run <- ComSim(pool, mode = mode, comsizes = c(2:nrow(pool$global.pool)), 
                    iter = iter, Kmeth = Kmeth)
      run$data <- cbind(run$data, regional.richness = rep(regmax - (i - 1), iter))
      regdat <- rbind(regdat, run$data)
    }
  }
  
  if (risk == "K") {
    for (i in 1:(regmax - 1)) {
      if (i == 1) {
        pool <- GPool(Nglobal = regmax, globmeth = globmeth, Bmeth = Bmeth, 
                      Bii = Bii, a = a, m = m, R0ii = R0ii, min.weight = min.weight)
      } else {
        del <- which.min(pool$global.pool[, 5]) # delete species with lowest K
        pool$global.pool <- pool$global.pool[-del, ]
      }
      run <- ComSim(pool, mode = mode, comsizes = c(2:nrow(pool$global.pool)), 
                    iter = iter, Kmeth = Kmeth)
      run$data <- cbind(run$data, regional.richness = rep(regmax - (i - 1), iter))
      regdat <- rbind(regdat, run$data)
    }
  }  
  class(regdat) <- "RegSim"
  regdat
}

#------------------------------------------------------------------------------
# Plotting function
detach("package:gam")
plot.RegSim <- function(run, method = "2d") { # Plots "RegSim" objects 
  if (method == "2d") {
    require(ggplot2)
    class(run) <- "data.frame"
    iter <- nrow(run)
    chart_title <- substitute(paste("Simulation results: ", iter, " iterations", sep=""),
                              list(iter=iter))
    return(ggplot(run, aes(regional.richness, richness, color=Ro)) +
      geom_jitter(alpha = ifelse(iter <= 1000, 1, exp(-0.00015 * iter) + .05)) +
      labs(x = "Regional richness", y = "Local richness") +
      lab(title = chart_title) +
      theme(panel.background = element_rect(colour = "black", fill = "white"),
            panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
              stat_abline(intercept = 1, slope = 0, lty = 2, col = "black"))
  }
  if (method == "3dsmooth") {  # Note that this requires more data points
    require(RcmdrPlugin.HH)
    return(scatter3d(run$regional.richness, run$Ro, run$richness, residuals = T,
                     model.summary = T, point.col = "black", fit = "smooth",
                     xlab = "Regional richness", ylab = "Ro", 
                     zlab = "Local richness", bg.col = "white"))
  }
}

#------------------------------------------------------------------------------
# Examples

# Allometric assignment of species traits
regA1 <- RegSim(regmax=10, globmeth="allom", Kmeth="fixed", iter=100,
              risk="random", mode="dens")
plot(regA1, "3dsmooth") # Local dilution, no regional effect

regA2 <- RegSim(regmax=10, globmeth="allom", Kmeth="fixed", iter=100,
                risk="K", mode="dens")
plot(regA2, "3dsmooth") 
# Local dilution, regional amplification, order of local extinctions important
# I.e. highly saturated communities may experience regional amplification effects

# Allometric random (a better null model than plain old random)
regAR1 <- RegSim(regmax=10, globmeth="allomr", Kmeth="fixed", iter=100,
                 risk="random", Bmeth="allomr", mode="dens")
plot(regAR1, "3dsmooth") # Local dilution, no regional effect

regAR2 <- RegSim(regmax=20, globmeth="allomr", Kmeth="fixed", iter=100,
                 risk="K", Bmeth="allomr", mode="dens")
plot(regAR2, "3dsmooth") # Local dilution, no regional effect