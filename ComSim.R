################
### ComSim() ###
################

### Purpose ###
# Randomly draw species from a global pool to produce simulated communities and
# calculate community R0 for each simulated community. 

### Input ###
# globalpool: output from the function GPool()
# mode: density dependent "dens" or frequency dependent "freq" transmission
# iter: how many iterated communities to produce
# comsizes: All possible values of species richness for simulated communities.
# Kmeth: Either "free" (total community carrying capacity is free to vary
#  as species are added), or "fixed" (community K is fixed, and abundances are
#  adjusted to sum to community K, while relative abundances do not change).

### Output ###
# A list containing relevant info on inputs including which global pool used,
# number of iterations, range of community richness values, community K method,
# and a data frame with info on each iteration. 

ComSim <- function(globalpool, mode, iter=1000, comsizes=c(2:20), Kmeth="free"){
  richness <- NULL
  Ro <- NULL
  shannondiv <- NULL
  density <- NULL
  iteration <- NULL
  J <- NULL
  composition <- matrix(NA, nrow = iter, ncol = max(comsizes))
  globalpool <- globalpool$global.pool # Extract species data from globalpool list
  
  if (Kmeth == "free") {
    for (j in 1:iter) {
      #Assign community members randomly
      comsize <- ifelse(length(comsizes) != 1, sample(comsizes, size=1), comsizes) 
      com <- c(sample(c(1:nrow(globalpool)), size = comsize, replace = F)) # ID's
      comtraits <- globalpool[com, ] # extract species traits
      pshan <- comtraits[, 5] / sum(comtraits[, 5]) # relative abundance
      shannondiv[j] <- -(sum((pshan) * log(pshan))) 
      richness[j] <- nrow(comtraits)
      density[j] <- sum(comtraits[, 5])
      iteration[j] = j  
      Ro[j] <- CommunityR0(comtraits, mode)
      composition[j, 1:length(com)] <- com
    }
  }
  
  if (Kmeth == "fixed") {
    for (j in 1:iter) {
      Kcom = 2 * max(globalpool[, 5])  # Community K = 2*max(K) for global pool
      #Assign community members randomly
      comsize <- ifelse(length(comsizes) != 1, sample(comsizes, size=1), comsizes)
      com <- c(sample(c(1:nrow(globalpool)), size = comsize, replace = F)) # ID's
      comtraits <- globalpool[com, ] # Get species traits for each community member
      # Adjust K's to make the sum of species K = community K, with abundances 
      # proportional to their relative abundances
      Kscale <- Kcom / sum(comtraits[, 5]) # Caclulate K scaling parameter
      comtraits[,5 ] <- comtraits[, 5] * Kscale # adjust densities
      pshan <- comtraits[, 5] / sum(comtraits[, 5]) # relative abundances
      shannondiv[j] <- -(sum((pshan) * log(pshan))) # still returning NA's?
      richness[j] <- nrow(comtraits)
      density[j] <- sum(comtraits[, 5])
      iteration[j] = j
      Ro[j] <- CommunityR0(comtraits, mode)
      composition[j, 1:length(com)] <- com
    }
  }
  
  # Start with all species, substract randomly until community K reached or passed
  if (Kmeth == "thresh") { 
    Kcom = 2 * max(globalpool[, 5])  # community K = 2*max(K) for global pool
    for (j in 1:iter) {
      comsize <- max(comsizes)
      com <- c(sample(c(1:nrow(globalpool)), size = comsize, replace = F))
      comtraits <- globalpool[com, ] # starting community
      while (sum(comtraits[, 5]) >= Kcom) { 
        del <- sample(nrow(comtraits), size = 1, replace = F) # choose a row
        comtraits <- comtraits[-del, ] # delete the row
      }
      pshan <- comtraits[, 5] / sum(comtraits[, 5]) # relative abundance
      shannondiv[j] <- -(sum((pshan) * log(pshan))) # still returning NA's?
      richness[j] <- nrow(comtraits)
      density[j] <- sum(comtraits[, 5])
      iteration[j] <- j
      Ro[j] <- CommunityR0(comtraits, mode)
      composition[j,1:length(com)] <- com
    }
  }
  
  
  Hmax <- log(richness)  
  J <- shannondiv / Hmax # again, returns NA's where shannondiv=NA. Need to fix
  dens.adj.rich <- resid(lm(richness ~ density)) # store density adjusted richness
  L <- list(pool = globalpool, 
            iter = iter,
            comsizes = comsizes, 
            Kmeth = Kmeth,
            data = as.data.frame(cbind(richness, Ro, shannondiv, J, density, iteration, 
                                       dens.adj.rich)), 
            composition = as.data.frame(composition))
  class(L) <- "ComSim"
  L
}
#--------------------------------------------------------------------------------------




# Plotting function
plot.ComSim <- function(run, method = "rich") { # plots "ComSim" objects
  # condenses the previous four functions: plot.raw(), plot.adj(), plot.even(), and plot.rich()
  require(ggplot2)
  iter <- max(run$data[, 6])
  chart_title <- substitute(paste("Simulation results: ", iter, " iterations", sep=""), 
                            list(iter=iter))
  
  if (method == "raw") {
    return(ggplot(run$data, aes(richness, density, color = Ro)) + 
      geom_jitter(alpha = ifelse(iter <= 1000, 1, exp(-0.00015 * iter) + .05)) +
      theme(panel.background = element_rect(colour = "black", fill = "white"), 
            panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
              labs(x = "Species richness", y = "Total density at carrying capacity") + 
              labs(title = chart_title) + 
              scale_color_gradient(low = "blue", high = "red"))
  }
  
  if (method == "adj") {
    return(ggplot(run$data, aes(dens.adj.rich, Ro, color = J)) +
      geom_jitter(alpha = ifelse(iter <= 1000, 1, exp(-0.00015 * iter) + .05))
           + geom_smooth() +
             theme(panel.background = element_rect(colour = "black", fill = "white"),
                   panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
                     scale_color_continuous(name = "Evenness") +
                     labs(x = "Density-adjusted species richness", y = "Ro") +
                     labs(title = chart_title))
  }
  
  if (method == "even") {
    return(ggplot(run$data, aes(J, Ro, color = dens.adj.rich)) +
      geom_jitter(alpha = ifelse(iter <= 1000, 1, exp(-0.00015 * iter) + .05)) +
      geom_smooth()+
      theme(panel.background = element_rect(colour = "black", fill = "white"), 
            panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
              scale_color_continuous(name = "Density-adjusted richness") + 
              labs(x = "Pielou's index of evenness", y = "Ro")+
              labs(title = chart_title))
  }
  
  if (method == "rich") {
    return(ggplot(run$data, aes(richness, Ro, color = J)) +
      geom_jitter(alpha = ifelse(iter <= 1000, 1, exp(-0.00015 * iter) + .05)) +
      labs(x = "Species richness", y = "Ro") +
      labs(title = chart_title) +
      scale_color_gradient(low = "black", high = "red") +
      coord_trans(y = "log10") +
      theme(panel.background = element_rect(colour = "black", fill = "white"),
            panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
              stat_abline(intercept = 1, slope = 0, lty = 2, col = "black"))
  }
}


#------------------------------------------------------------------------------
# Examples

poolA<-GPool(globmeth="allom")
run<-ComSim(poolA, iter=1000, mode="freq", Kmeth="fixed")
plot(run, method="raw")
plot(run, method="adj")
plot(run, "even")
plot(run, "rich")
