#############################################
## Multihost SIR model, matrix formulation ##
#############################################
# Last edited Oct 26, 2012 MJ

# Purpose: Run a multi-host species SIR model out through time.
# Input: A GPool object. We will probably want to decide what type of input 
#   is ideal for our purposes (perhaps we want to run the model on a 
#   simulated community from comsim())
# Output: Numerical solution to the system of differential equations that 
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

#     [b1-d1      b1    b1      0         0     0] (b=births, d=deaths)
#     [0     -(d1+g1)    0      0         0     0] (g=recovery)
#     [0          g1   -d1      0         0     0]
# B = [0           0     0  b2-d2        b2    b2]
#     [0           0     0      0   -(d2+g2)    0]
#     [0           0     0      0        g2   -d2]

# Notes: 1) This script requires GPool
#        2) How to incorportate frequency dependent transmission? 
#             Divide transmission term by density of all species?

library(deSolve)

# Produce species parameters with GPool() 
pool <- GPool(Nglobal=2, globmeth="allomr")
pool$global.pool[,3]=pool$global.pool[,2] # d=b according to Roche et al. 2012

# At some point, we may decide to assign contact rates and susceptibilities 
#   according to Roche et al. 2012 (see Joe's script)
# Current version uses Bii terms from the gpool() function

SIRmatrix <- function(t, X, parms){
  with(parms, {
    dX <- (B %*% X) + X[vec]*(A %*% X)
    return(list(c(dX)))
  })
}

comtraits <- pool$global.pool
cij <- 0.00005 # Assign strength of interspecific transmission
B.1 <- matrix(rep(comtraits[, 8], nrow(comtraits)), 
              nrow=nrow(comtraits), ncol=nrow(comtraits))
B.2 <- matrix(rep(comtraits[, 8], each=nrow(comtraits)), 
              nrow=nrow(comtraits), ncol=nrow(comtraits))
B.m=cij * (B.1 + B.2) / 2
diag(B.m) <- comtraits[, 8] # specify Bii terms along the diagonal
stopifnot(B.m <= 1 & B.m >= 0) # Beta has to be between 0 and 1

# Insert Beta values into matrix A
A <- matrix(0, nrow = pool$Nglobal * 3, 
           ncol=pool$Nglobal*3, byrow=T)
A.cs <- seq(2, nrow(A), by=3) # indices for rows to be filled
A.rs <- seq(1, ncol(A)) # indices for columns to be filled
A.rs <- A.rs[- which(A.rs %% 3 == 0)] # Remove the rows for the recovered

A[A.rs, A.cs]<-rep(B.m, each=2) # Fill the A matrix with Beta values
sign <- rep(c(-1, 1, 0),length.out=ncol(A)) # -1 for S, 1 for I, 0 for R
A <- sign * A 

B <- matrix(0, nrow=pool$Nglobal*3, ncol=pool$Nglobal*3, byrow=T)
B.init <- NULL
B.vec <- NULL
 
for(i in 1:pool$Nglobal){
  B.vec <- c(B.vec, pool$global.pool[i, 2] - pool$global.pool[i, 3], 
             pool$global.pool[i, 2], pool$global.pool[i, 2],
             0, -(pool$global.pool[i, 3] + pool$global.pool[i, 7]), 0,
             0, pool$global.pool[i, 7], - pool$global.pool[i, 3])
}

B.init <- matrix(B.vec, nrow=3 * pool$Nglobal, ncol=3, byrow=T)
strt <- seq(1, 3 * pool$Nglobal, by=3) # index demographic parameters
stp<-seq(3, 3 * pool$Nglobal, by=3) 

for(i in 1:pool$Nglobal){
  B[strt[i]:stp[i], strt[i]:stp[i]] <- B.init[strt[i]:stp[i], 1:3]
}

parms <- list(A, B)

## Specify initial conditions
index.case <- 1 # which species has the one infected individual? (species ID)
X <- rep(NA, length.out=pool$Nglobal * 3) # vector of abundances at time=0

S <- pool$global.pool[, 5] # carrying capacities for each species (So)
I <- rep(0, length.out=pool$Nglobal) 
I[index.case] <- 1 
S[index.case] <- S[index.case] - 1 # adjust S so that K is unchanged
R <- rep(0, length.out=pool$Nglobal) # (naiive populations)

X <- c(rbind(S, I, R)) # S1, I1, R1, S2, I2, R2,..., Sn, In, Rn

# initial condition vector names
spec.num <- pool$global.pool[, 1]
class.nam <- c("S", "I", "R") 
names(X) <- apply(expand.grid(class.nam, spec.num), 1, paste, collapse = "", sep = "") 

# indices of X to multiply by A %*% X to get the "S" in the "SI"
vec <- rep(seq(from=1, to=length(X), by=3),
           each=3, length=length(X)) # indices of the "S" classes 

## Define vector of timesteps 
times  <- seq(0, 10, length=100)


######################################
# Run the model and plot the results #
######################################
out <- ode(X, times, SIRmatrix, parms)
head(out)
plot(out)

