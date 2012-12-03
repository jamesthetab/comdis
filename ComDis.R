################
### ComDis() ###
################
# Last edited Nov 1 to allow use of kWeightPenalty
### Purpose ###
# Iteratively disassemble communities until they reach 2 species.

### Input ###
# globalpool: output from the function GPool()
# mode: density dependent "dens" or frequency dependent "freq" transmission
# iter: how many iterated communities to produce
# Kmeth: Either "free" (total community carrying capacity is free to vary
#  as species are added), or "fixed" (community K is fixed, and abundances are
#  adjusted to sum to community K, while relative abundances do not change).
# start: maximum species richness
# kWeightPenalty: A constant that determines the strength of weight-dependence 
#  in the sampling of communities. If kWeightPenalty > 0, small bodied species
#  are more likely to be sampled. If kWeightPenalty = 0, then there is no bias.

### Output ###
# A list containing the Kmethod used, relevant output for each iteration, 
# and community composition for each iteration.

ComDis <- function(globalpool, mode,
                   iter = 100, Kmeth = "free", 
									 kWeightPenalty = 0) {
  richness <- NULL
  Ro <- NULL
  delta.Ro <- NULL
  shannondiv <- NULL
  density <- NULL
  outer.iteration <- NULL
  inner.iteration <- NULL
  J <- NULL
  globalpool <- globalpool$global.pool
  start <- nrow(globalpool)
  all.data <- NULL
  all.composition <- NULL
  
  for (j in 1:iter) {
    com <- sample(1:nrow(globalpool), size = start, replace = F)
    out <- matrix(NA, ncol = length(com), nrow = length(com) - 1)
    # Inner for-loop disassembles the community until there are 2 species
    for (i in 1:(length(com) - 1)) {
      if (i == 1) { # for the first iteration, we start with out starting community
        out[i, 1:(length(com))] <- com
      } else { # otherwise, delete one species 
        out[i, 1:(length(com) - (i - 1))] <- sample(out[(i - 1), which(out[(i - 1), ] != "NA")], 
                                                    size = (length(com) - (i - 1)), replace = F,
        																						prob = traits$weight ^ - kWeightPenalty)
      }
      ids <- out[i, which(out[i, ] != "NA")] # Store species ID's
      traits <- globalpool[ids, ] # extract species traits
      if (Kmeth == "fixed"){
      	Kcom =2 * max(globalpool[, 5])
      	Kscale <- Kcom / sum(traits[, 5]) # caclulate K scaling parameter
      	traits[, 5] <- traits[, 5] * Kscale # scale abundances
      }
      pshan <- traits[, 5] / sum(traits[, 5]) # relative abundance
      shannondiv[i] <- -(sum((pshan) * log(pshan))) 
      richness[i] <- nrow(traits)
      density[i] <- sum(traits[, 5])
      inner.iteration[i] = i  
      outer.iteration = rep(j, length(inner.iteration))
      Ro[i] <- CommunityR0(traits, mode)
      delta.Ro[i] <- ifelse(i == 1, NA, Ro[i-1] - Ro[i])
      df <- cbind(outer.iteration, inner.iteration, richness, density, shannondiv, Ro, delta.Ro)
    }
    all.data <- rbind(all.data, df)
    all.composition <- rbind(all.composition, out)   
  }
  
  J <- all.data[, 5] / log(all.data[, 3]) # evenness
  all.data <- cbind(all.data, J)
  d <- list(Kmeth = Kmeth, pool=globalpool, all.data = all.data, 
            all.composition = all.composition)
  class(d) <- "ComDis"
  d
}


#------------------------------------------------------------------------------
# Plotting function

plot.ComDis <- function(run) { # Plots disassembly trajectories for "ComDis" objects
  require(ggplot2)
  iter <- max(run$all.data[, 1])
  chart_title <- substitute(paste("Simulation results: ", iter, " iterations", sep=""),
                            list(iter = iter))
  ggplot(data.frame(run$all.data)) + 
    geom_line(aes(richness, Ro, group = outer.iteration),
    					alpha = ifelse(iter <= 50, 1, exp(-0.005 * iter) + .05)) +
    labs(x = "Species richness", y = "Ro") +
    labs(title = chart_title) +
    coord_trans(y = "log10") +
    theme(panel.background = element_rect(colour = "black", fill = "white"),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
            stat_abline(intercept = 1, slope = 0, lty = 2, col = "black")
}

#------------------------------------------------------------------------------
# Example
poolA <- GPool(Nglobal=6, globmeth="allom")
dis1<-ComDis(poolA,iter=3,Kmeth="free", mode="freq")
plot(dis1)