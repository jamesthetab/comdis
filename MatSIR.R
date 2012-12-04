##############
## MatSIR() ##
##############
# Last edited Nov 11, 2012 MJ: Included virulence term, turned the whole thing 
# into a function with comtraits as an input and maximal 
# disease prevalence as an output.

### Purpose ### 
# Run a multi-host species SIR model out through time.

### Input ### 
# comtraits: A comtraits matrix
# mode: density dependent "dens" or frequency dependent "freq" transmission
# cij: a scaling parameter that determines strength of interspecific transmission
#  where Bij = cij * (Bii+Bjj)/2 following Dobson 2004
# index.case: which species has the 1 infected individual initially, defaults to 
#   random species (with equal probabilities)
# times: specification of timesteps for the ode solver
 
### Output ### 
# Equilibrium and maximal disease prevalance as indicated by
# numerical solution to the system of differential equations that 
#   comprise the SIR model, represented in matrix form as:
# dX <- (B %*% X) + X[vec]*(A %*% X)

# e.g. in a 2-species case
# X is a vector of classes for each species (S1, I1, R1, S2, I2, R2)
# X[vec] is a vector of susceptible classes (S1, S1, S1, S2, S2, S2)

#     [0 -B11  0  0  -B12  0]
#     [0  B11  0  0   B12  0]
# A = [0  0    0  0   0    0] (transmission)
#     [0 -B21  0  0  -B22  0]
#     [0  B21  0  0   B22  0]

#     [b1-d1      b1       b1      0         0     0] (b=births, d=deaths)
#     [0     -(d1+g1+v1)    0      0         0     0] (g=recovery)
#     [0          g1      -d1      0         0     0] (v=virulence)
# B = [0           0        0  b2-d2        b2    b2]
#     [0           0     0      0   -(d2+g2+v2)    0]
#     [0           0        0      0        g2   -d2]

### Notes ###
# The current formulation is only valid when birth rates equal death rates

MatSIR <- function(comtraits, mode,
									 cij=0.05,
									 index.case=sample(1:nrow(comtraits), size=1),
									 times=seq(0, 500, length=100)){
	require(deSolve)
	if(mode == "freq"){
		# First correct Bii values to account for independence from density
		comtraits[, 8] <- comtraits[, 8] * comtraits[, 5] # cancel out density: Bii = (K*(R0*(d + g +V)))/K
	}
	if (all(comtraits[, 2] != comtraits[, 3])) {
		comtraits[, 3] <- comtraits[, 2]
		print("Death rates set to equal birth rates")
	} # b must = d in the current formulation

	SIRmatrix <- function(t, X, parms){
  	with(parms, {
    	dX <- (B %*% X) + X[vec]*((A*p) %*% X)
    	return(list(c(dX)))
  	})
	}

	B.1 <- matrix(rep(comtraits[, 8], nrow(comtraits)), 
  	            nrow=nrow(comtraits), ncol=nrow(comtraits))
	B.2 <- matrix(rep(comtraits[, 8], each=nrow(comtraits)), 
    	          nrow=nrow(comtraits), ncol=nrow(comtraits))
	B.m <- cij * (B.1 + B.2) / 2
	diag(B.m) <- comtraits[, 8] # specify Bii terms along the diagonal
	stopifnot(B.m >= 0) # Beta can't be negative

	# Insert Beta values into matrix A
	A <- matrix(0, nrow = nrow(comtraits) * 3, 
  	         ncol = nrow(comtraits) * 3, byrow = T)
	A.cs <- seq(2, nrow(A), by=3) # indices for rows to be filled
	A.rs <- seq(1, ncol(A)) # indices for columns to be filled
	A.rs <- A.rs[- which(A.rs %% 3 == 0)] # Remove the rows for the recovered

	A[A.rs, A.cs]<-rep(B.m, each=2) # Fill the A matrix with Beta values
	sign <- rep(c(-1, 1, 0),length.out=ncol(A)) # -1 for S, 1 for I, 0 for R
	A <- sign * A 

	B <- matrix(0, nrow=nrow(comtraits)*3, ncol=nrow(comtraits)*3, byrow=T)
	B.init <- NULL
	B.vec <- NULL
 
	for(i in 1:nrow(comtraits)){
  	B.vec <- c(B.vec, 
  						 comtraits[i, 2] - comtraits[i, 3], 
      	       comtraits[i, 2], 
  						 comtraits[i, 2],
          	   0, 
  						 -(comtraits[i, 3] + comtraits[i, 7] + comtraits[i, 6]), 
  						 0,
           	   0, 
  						 comtraits[i, 7], 
  						 - comtraits[i, 3])
	}

	B.init <- matrix(B.vec, nrow=3 * nrow(comtraits), ncol=3, byrow=T)
	strt <- seq(1, 3 * nrow(comtraits), by=3) # index demographic parameters
	stp<-seq(3, 3 * nrow(comtraits), by=3) 

	for(i in 1:nrow(comtraits)){
  	B[strt[i]:stp[i], strt[i]:stp[i]] <- B.init[strt[i]:stp[i], 1:3]
	}

	parms <- list(A, B)

	## Specify initial conditions
	X <- rep(NA, length.out=nrow(comtraits) * 3) # vector of abundances at time=0

	S <- comtraits[, 5] # carrying capacities for each species (So)
	I <- rep(0, length.out=nrow(comtraits)) 
	I[index.case] <- 1 
	S[index.case] <- S[index.case] - 1 # adjust S so that K is unchanged
	R <- rep(0, length.out=nrow(comtraits)) # (naiive populations)

	X <- c(rbind(S, I, R)) # S1, I1, R1, S2, I2, R2,..., Sn, In, Rn

	# initial condition vector names
	spec.num <- 1:nrow(comtraits)
	class.nam <- c("S", "I", "R") 
	names(X) <- apply(expand.grid(class.nam, spec.num), 1, paste, sep="", collapse="") 

	# indices of X to multiply by A %*% X to get the "S" in the "SI"
	vec <- rep(seq(from=1, to=length(X), by=3),
           each=3, length=length(X)) # indices of the "S" classes 
	# specify a p term to allow for density or frequency dep transmission
	p <- ifelse(mode=="freq", 1/sum(X), 1)

	out <- ode(X, times, SIRmatrix, parms) # run model
	
	# calculate maximal prevalance
	all.cols <- 1:ncol(out)
	I.cols <- all.cols[-which(all.cols %% 3 != 0)] 	#col indices for I classes
	N.infected <- apply(out[, I.cols], 1, sum)
	N.total <- apply(out[, -1], 1, sum) # all columns except 1st (time)
	prev <- N.infected / N.total # prevalence at each time
	end.prev <- prev[length(prev)]
	max.prev <- max(prev) # maximum prevalence across all species
	return(c(max.prev, end.prev))
}