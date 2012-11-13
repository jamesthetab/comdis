################
### ComSim() ###
################
# Edited Nov. 11 MJ: Added calculation of maximum prevalence when birth rates =
# death rates (as with HybridGPool), progress bar, cij as function argument.

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
# SampleMass: Can be TRUE or FALSE, determines whether sampling of species occurs
#   with probabilities proportional to abundance, with more abundant species
#   more likely to be sampled. 
# cij: a scaling parameter that determines strength of interspecific transmission
#  where Bij = cij * (Bii+Bjj)/2 following Dobson 2004

### Output ###
# A list containing relevant info on inputs including which global pool used,
# number of iterations, range of community richness values, community K method,
# whether SampleMass was T or F, and a data frame with info on each iteration. 

ComSim <- function(globalpool, mode, iter=1000, comsizes=c(2:globalpool$Nglobal), 
									 Kmeth="free", SampleMass=F, cij=0.05){
  richness <- NULL
  Ro <- NULL
  shannondiv <- NULL
  density <- NULL
  iteration <- NULL
  J <- NULL
  maxprev <- NULL
  composition <- matrix(NA, nrow = iter, ncol = max(comsizes))
  globalpool <- globalpool$global.pool # Extract species data from globalpool list
  bar <- txtProgressBar (min = 0, max = iter, style = 3)
  
  if (Kmeth == "free") {
    for (j in 1:iter) {
      #Assign community members randomly
      comsize <- ifelse(length(comsizes) != 1, sample(comsizes, size=1), comsizes) 
      if (SampleMass == TRUE){ #sample with probability proportional to abundance
        com <- sample(c(1:nrow(globalpool)), size = comsize, replace = F, 
        							prob = (globalpool[,5]/100))
      } else {
        com <- sample(c(1:nrow(globalpool)), size = comsize, replace = F) # ID's
      }
      comtraits <- globalpool[com, ] # extract species traits
      pshan <- comtraits[, 5] / sum(comtraits[, 5]) # relative abundance
      shannondiv[j] <- -(sum((pshan) * log(pshan))) 
      richness[j] <- nrow(comtraits)
      density[j] <- sum(comtraits[, 5])
      iteration[j] = j  
      Ro[j] <- CommunityR0(comtraits, mode, cij=cij)
      maxprev[j] <- ifelse(all(comtraits[, 2] != comtraits[, 3]), # only valid when b=d
      										 NA, MatSIR(comtraits, mode, cij=cij))
      composition[j, 1:length(com)] <- com
      setTxtProgressBar (bar, j)
    }
  }
  
  if (Kmeth == "fixed") {
    for (j in 1:iter) {
      Kcom = 3 * max(globalpool[, 5])  # Community K = 2*max(K) for global pool
      #Assign community members randomly
      comsize <- ifelse(length(comsizes) != 1, sample(comsizes, size=1), comsizes)
      if (SampleMass == TRUE){ # sample with probability proportional to abundance
        com <- sample(c(1:nrow(globalpool)), size = comsize, replace = F, 
        							prob = (globalpool[,5]/100))
      } else {
        com <- sample(c(1:nrow(globalpool)), size = comsize, replace = F) # ID's
      } 
      comtraits <- globalpool[com, ] # Get species traits for each community member
      # Adjust K's to make the sum of species K = community K, with abundances 
      # proportional to their relative abundances
      Kscale <- Kcom / sum(comtraits[, 5]) # Caclulate K scaling parameter
      comtraits[, 5] <- comtraits[, 5] * Kscale # adjust densities
      pshan <- comtraits[, 5] / sum(comtraits[, 5]) # relative abundances
      shannondiv[j] <- -(sum((pshan) * log(pshan))) # still returning NA's?
      richness[j] <- nrow(comtraits)
      density[j] <- sum(comtraits[, 5])
      iteration[j] = j
      Ro[j] <- CommunityR0(comtraits, mode, cij=cij)
      maxprev[j] <- ifelse(all(comtraits[, 2] != comtraits[, 3]), # only valid when b=d
      										 NA, MatSIR(comtraits, mode, cij=cij))
      composition[j, 1:length(com)] <- com
      setTxtProgressBar (bar, j)
    }
  }
  
  if (Kmeth == "saturate.1" | Kmeth == "saturate.2" | Kmeth == "saturate.3") { # curvilinear saturating relationship between
  	for (j in 1:iter) {      # species richness & community density
  		#Assign community members randomly
  		comsize <- ifelse(length(comsizes) != 1, sample(comsizes, size=1), comsizes)
  		if (SampleMass==TRUE){ #sample with probability proportional to abundance
  			com <- sample(c(1:nrow(globalpool)), size = comsize, replace = F, 
  										prob = (globalpool[,5]/100))
  		} else {
  			com <- sample(c(1:nrow(globalpool)), size = comsize, replace = F) # ID's
  		} 
  		comtraits <- globalpool[com, ] # Get species traits for each community member
  		# Adjust K's to make the sum of species K = community K, with abundances 
  		# proportional to their relative abundances
      if (Kmeth == "saturate.1"){
        Kcom <- 500 - 3100 / (comsize + 5)
      }
      if (Kmeth == "saturate.2"){
        Kcom = 500 / (1 + 50*exp(-0.15*(comsize+10)))
      }
      if (Kmeth == "saturate.3"){
        Kcom = 500 / (1 + 200*exp(-0.15*(comsize+10)))
      }
        
      if ( sum(comtraits[, 5]) < Kcom ){
         comtraits[, 5] <- comtraits[, 5]
      } else {
        Kscale <- Kcom / sum(comtraits[, 5]) # Caclulate K scaling parameter
        comtraits[,5 ] <- comtraits[, 5] * Kscale # adjust densities
      }
      
  		pshan <- comtraits[, 5] / sum(comtraits[, 5]) # relative abundances
  		shannondiv[j] <- -(sum((pshan) * log(pshan))) # still returning NA's?
  		richness[j] <- nrow(comtraits)
  		density[j] <- sum(comtraits[, 5])
  		iteration[j] = j
  		Ro[j] <- CommunityR0(comtraits, mode, cij=cij)
  		maxprev[j] <- ifelse(all(comtraits[, 2] != comtraits[, 3]), # only valid when b=d
  												 NA, MatSIR(comtraits, mode, cij=cij))
  		composition[j, 1:length(com)] <- com
  		setTxtProgressBar (bar, j)
  	}
  }
    
  Hmax <- log(richness)  
  J <- shannondiv / Hmax # again, returns NA's where shannondiv=NA. Need to fix
  if(Kmeth=="saturate"){
  	dens.adj.rich <- resid(lm(density ~ richness + I(richness^2)))
  } else {
  	dens.adj.rich <- resid(lm(density ~ richness))
  } # store density adjusted richness
  
  L <- list(pool = globalpool, 
            iter = iter,
            comsizes = comsizes, 
            Kmeth = Kmeth,
            SampleMass = SampleMass,
            data = as.data.frame(cbind(richness, Ro, maxprev, shannondiv, J, density, iteration, 
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

poolA <- HybridGPool()
run1 <- ComSim(poolA, iter=100, mode="dens", Kmeth="free", SampleMass=F, cij=0.05)

# pre-programmed plotting functionality
plot(run1, method="raw")
plot(run1, method="adj")
plot(run1, "even")
plot(run1, "rich")

# comparing Ro and maximum prevalence
ggplot(run1$data, aes(x=Ro, y=maxprev)) + 
	geom_point() + geom_smooth() + theme_bw()

# plotting richness against maximum prevalence
ggplot(run1$data, aes(x=richness, y=maxprev)) + 
	geom_point() + geom_smooth() + theme_bw()


run2 <- ComSim(poolA, iter=100, mode="freq", Kmeth="free", SampleMass=F)
plot(run2)

ggplot(run2$data, aes(x=Ro, y=maxprev)) + 
	geom_point() + geom_smooth() + theme_bw()

ggplot(run2$data, aes(x=richness, y=maxprev)) + 
	geom_point() + geom_smooth() + theme_bw()
