#############
# ComPerm() #
#############
# Nov 1 MJ updated to allow larger global pools. Note that this tends to 
# make permutations that have intermediate numbers of inversions, so the 
# transitions from non-randomness to randomness are more abrupt. Perhaps
# it's best to stick with small global pool sizes (up to 6 or 7 species).
# Also added a plotting function with 3d or facetted option.

### Purpose ###
# ComPerm repeatedly permutes a global species pool and simulates communities
# from each permutation. Intraspecific transmission (Bii) values are permuted, 
# and shuffled among species to generate a range of possibilities. When there 
# are 0 inversions, then Bii values match allometry as in Dobson (2004), when
# there are a maximal number of inversions, then Bii values are reversed, so 
# that small species (low weight) have the highest Bii values, rather than the
# lowest as in Dobson (2004).

### Inputs ###
# pool: An object that contains the output of GPool()
# kWeightPenalty: A constant that determines the strength of weight-dependence 
#  in the sampling of communities. If kWeightPenalty > 0, small bodied species
#  are more likely to be sampled. If kWeightPenalty = 0, then there is no bias.
# mode: Either frequency ("freq") or density dependent ("dens") transmission.
# comsizes: All possible values of species richness for simulated communities.
# iter: # of iterated communities to produce for each global pool permutation
# Kmeth: Either "free" (total community carrying capacity is free to vary
#  as species are added), or "fixed" (community K is fixed, and abundances are
#  adjusted to sum to community K, while relative abundances do not change).

### Outputs ###
# An object of class "ComPerm" that is basically a list containing
#  pool: all of the global pool permutations
#  iter: the number of iterated communities per global pool permutation
#  comsizes: the range of species richnesses used
#  K meth: The method used to assign abundances
#  data: calculated results for each iteration
#  all.comp: the composition of species for each iteration (species ID's)

ComPerm <- function(pool, kWeightPenalty, mode,
										comsizes = 2:pool$Nglobal, iter = 1000, Kmeth = "free",
										Nperms = 1000) {   
	stopifnot(Kmeth == "fixed" | Kmeth == "free")
	# Generate all possible Bii permutations from a global pool
	require(gtools)
	if (pool$Nglobal < 7){
		perms <- permutations(pool$Nglobal, pool$Nglobal) # Bii permutation indices
	} else {
		perms <- matrix(NA, nrow=Nperms, ncol=pool$Nglobal) # can control # perms
		perms[1, ] <- 1:pool$Nglobal # first permutation has 0 inversions
		perms[Nperms, ] <- pool$Nglobal:1 # last has max(# inversions)
		perms[2:(Nperms - 1),] <- permute(1:pool$Nglobal) # rest are random perms
	}
	Bii.vec <- pool$global.pool[, 8] #allometric scaling Bii vector
	Bii.perms <- matrix(Bii.vec[perms[1:length(perms[, 1]), ]], 
											ncol = pool$Nglobal)# permutations of the Bii vector
	
	inversions <- function(init) { # counts inversions in a permutation
		inv_count <- 0
		for(i in 1:(length(init) - 1)) {
			for (j in (i + 1):(length(init - 1))) {
				if (init[i] > init[j]) {
					inv_count <- inv_count + 1
				}
			}
		}
		inv_count
	}
	
	Bii.inversions <- apply(Bii.perms, 1, inversions) # calculate number of inversions
	
	# For each permutation except the first, create new global pool
	global.perms <- pool$global.pool
	for (i in 2:length(perms[, 1])) { 
		perm.i <- pool$global.pool
		perm.i[, 8] <- Bii.perms[i, ]
		global.perms <- rbind(global.perms, perm.i) # create all global permutations
	}
	
	global.perms$permutation <- rep(1:length(perms[, 1]), # store permutation ID
																	each = pool$Nglobal, 
																	length.out = length(global.perms[, 1]))
	global.perms$inversions <- rep(Bii.inversions, each = pool$Nglobal, 
																 length.out = length(global.perms[, 1])) 
	
	# Simulate communities and Ro for each global pool permutation
	number.perms <- max(global.perms$permutation) # total number of permutations
	all.richness <- NULL
	all.Ro <- NULL
	all.shannondiv <- NULL
	all.density <- NULL
	all.iteration <- NULL
	all.J <- NULL
	all.dens.adj.rich <- NULL
	all.comp <- NULL
	perm <- NULL
	invers<-NULL
	bar <- txtProgressBar (min = 0, max = number.perms, style = 3)
	
	for (p in 1:number.perms) { # for each permutation
		richness <- NULL
		Ro <- NULL
		shannondiv <- NULL
		density <- NULL
		iteration <- NULL
		J <- NULL
		Hmax <- NULL
		composition <- matrix(NA, nrow = iter, ncol = max(comsizes))
		globalpool <- subset(global.perms, permutation == p) 
		
		for (j in 1:iter) {
			comsize <- ifelse(length(comsizes) != 1, sample(comsizes, size = 1), comsizes) 
			# Choose which species will be in the community based on weight 
			# (if kWeightPenalty != 0, then sampling is based on weight)
			com.weights <- sample(globalpool$weight, size = comsize, replace = F, 
														prob = globalpool$weight ^ - kWeightPenalty) 
			comtraits <- subset(globalpool, globalpool$weight %in% com.weights)
			if (Kmeth == "fixed") {
				# Adjust abundances to make the sum of species K = community K, 
				# without changing relative abundances
				Kcom = 2 * max(globalpool[, 5])  # community K = 2*max(K)
				Kscale <- Kcom / sum(comtraits[, 5]) # K scaling parameter
				comtraits[, 5] <- comtraits[, 5] * Kscale # adjust abundances
			}
			pshan <- comtraits[, 5] / sum(comtraits[, 5]) # relative abundance
			shannondiv[j] <- -(sum((pshan) * log(pshan))) 
			richness[j] <- nrow(comtraits)
			density[j] <- sum(comtraits[, 5])
			iteration[j] <- j  
			Ro[j] <- CommunityR0(comtraits, mode)
			composition[j, 1:length(comtraits[, 1])] <- comtraits[, 1]
		}
		
		# Gather output
		Hmax <- c(Hmax, log(richness))  
		J <- c(J, shannondiv / Hmax) # evenness
		dens.adj.rich <- resid(lm(richness ~ density)) # density adjusted richness
		all.richness <- c(all.richness, richness)
		all.Ro <- c(all.Ro, Ro)
		all.shannondiv <- c(all.shannondiv, shannondiv)
		all.J <- c(all.J, J)
		all.density <- c(all.density, density)
		all.iteration <- c(all.iteration, iteration)
		all.dens.adj.rich <- c(all.dens.adj.rich, dens.adj.rich)
		all.comp <- rbind(all.comp, composition)
		perm <- c(perm, rep(p, iter)) # permutation number
		invers <- c(invers, rep(Bii.inversions[p], iter)) # number of inversions
		setTxtProgressBar (bar, p)
	}
	
	data = data.frame(all.richness, all.Ro, all.shannondiv, all.J, all.density, 
										all.iteration, all.dens.adj.rich, perm, invers)
	
	L <- list(pool = global.perms,
						iter = iter,
						comsizes = comsizes, 
						Kmeth = Kmeth, 
						data = data, 
						all.comp = as.data.frame(all.comp))
	
	class(L) <- "ComPerm"
	L
}

#------------------------------------------------------------------------------
# Plotting function
plot.ComPerm <- function(run, method="facet"){
	if (method == "3d"){
		require(RcmdrPlugin.HH); detach("package:gam")
		p <- scatter3d(x=run$data$invers, y=run$data$all.Ro, z=run$data$all.richness, 
									 fit="smooth", xlab="Relationship between extinction risk and Bii",
									 ylab="Community R0",
									 zlab="Species Richness")	
	}
	
	if (method == "facet"){
		require(ggplot2)
		iter <- run$iter
		p <- ggplot(run$data) + theme_bw() +
			geom_point(aes(x=all.richness, y=all.Ro, color=as.factor(perm)), 
								 alpha = ifelse(iter <= 22, 1, 
								 							 (exp(-0.01 * iter) + .2))) +
			geom_smooth(aes(x=all.richness, y=all.Ro, color=as.factor(perm)), alpha=0.1) +
			facet_wrap(~invers) + theme(legend.position="none") 
	}
	return(suppressWarnings(suppressMessages(print(p))))
}

#------------------------------------------------------------------------------
# Examples

poolA<-GPool(Nglobal=6, globmeth="allom", a=2)
test<-ComPerm(poolA, kWeightPenalty=10, iter=20, Kmeth="free", mode="freq", Nperms=20)
plot(test)