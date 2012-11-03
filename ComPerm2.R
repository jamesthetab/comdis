##############
# ComPerm2() #
##############

### Purpose ###
# ComPerm2 repeatedly permutes a global species pool and simulates community
# disassembly trajectories from each permutation. Intraspecific transmission 
# (Bii) values are permuted, and shuffled among species to generate a range 
# of possibilities. When there are 0 inversions, then Bii values match 
# allometry as in Dobson (2004), when there are a maximal number of inversions, 
# then Bii values are reversed, so that small species (low weight) have the 
# highest Bii values. This differs from ComPerm() in that it produces 
# disassembly trajectories consisting of sequential extinctions, rather than 
# Simply producing simulated communities.

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
# Nperms: The number of global pool permutations produced if Nglobal<6. 

### Outputs ###
# An object of class "ComPerm2" that is basically a list containing
#  pool: all of the global pool permutations
#  iter: the number of iterated communities per global pool permutation
#  comsizes: the range of species richnesses used
#  K meth: The method used to assign abundances
#  data: calculated results for each iteration
#  all.comp: the composition of species for each iteration (species ID's)

ComPerm2 <- function(pool, kWeightPenalty, mode,
										comsizes = 2:pool$Nglobal, iter = 10, Kmeth = "free",
										Nperms = 1000) {   
	stopifnot(Kmeth == "fixed" | Kmeth == "free")
	write("generating permutations...", "")
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
	data <- NULL
	all.comp <- NULL
	write("running simulations...", "")
	bar <- txtProgressBar (min = 0, max = number.perms, style = 3)
	
	for (p in 1:number.perms){
		globalpool <- subset(global.perms, permutation == p)
		dummylist <- list(global.pool = globalpool) # because ComDis extracts the pool from a list
		runp <- ComDis(dummylist, mode, iter, Kmeth, kWeightPenalty)
		permutation <- rep(p, nrow(runp$all.data))
		inversions <- rep(Bii.inversions[p], nrow(runp$all.data))
		rund <- cbind(runp$all.data, permutation, inversions)
		data <- rbind(data, rund)
		all.comp <- rbind(all.comp, runp$all.composition)
		setTxtProgressBar (bar, p)
	}
	
	L <- list(pool = global.perms,
						iter = iter,
						comsizes = comsizes, 
						Kmeth = Kmeth, 
						data = data.frame(data), 
						all.comp = as.data.frame(all.comp))
	
	class(L) <- "ComPerm2"
	L
}

#------------------------------------------------------------------------------
# Plotting function
plot.ComPerm2 <- function(run){
	require(ggplot2)
	iter <- max(run$data[, 1])
	chart_title <- substitute(paste("Disassembly trajectories: ", iter, " per permutation", sep=""),
														list(iter = iter))
	require(ggplot2)
	ggplot(run$data) + theme_bw() +
		geom_line(aes(richness, Ro, group = interaction(outer.iteration, permutation),
									color=as.factor(permutation)),
							alpha=ifelse(max(run$data$outer.iteration) <= 50, 1, (exp(-0.01 * max(run$data$outer.iteration)) + .05))) +
								facet_wrap(~inversions) +
								theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
								theme(legend.position="none") +
								ggtitle(chart_title)
}

#------------------------------------------------------------------------------
# Example
pool1 <- GPool(Nglobal=6, globmeth="allom", a=2)
test <- ComPerm2(pool1, mode="freq", kWeightPenalty=2, iter=10, Kmeth="free")
plot(test)