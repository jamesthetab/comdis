##################
# CommunityR0()  #
##################
### Purpose ###
# Calculate community-level reproductive ratio (R0) of a multihost pathogen

### Input ###
# comtraits: matrix of community traits
# mode: density dependent "dens" or frequency dependent "freq" transmission
# cij: parameter that determines the strength of interspecific transmission
#   relative to intraspecific transmission

### Output ###
# Dominant eigenvalue of the next generation matrix G (following Dobson 2004)

CommunityR0 <- function(comtraits, mode, cij = 0.05){
  if (mode == "freq"){
  	# First correct Bii values to account for independence from density
  	# comtraits[, 8] <- comtraits[, 8] * comtraits[, 5] # cancel out density: Bii = (K*(R0*(d + g +V)))/K
    # above line commented out because adjustment is now made when first assigning traits
  	B.1 <- matrix(rep(comtraits[, 8], nrow(comtraits)), nrow = nrow(comtraits), ncol = nrow(comtraits))
  	B.2 <- matrix(rep(comtraits[, 8], each = nrow(comtraits)), nrow = nrow(comtraits), ncol=nrow(comtraits))
  	B.m = cij * (B.1 + B.2) / 2
  	diag(B.m) <- comtraits[, 8] # Specifies intraspecific transmission terms along the diagonal
  	stopifnot (B.m >= 0) # Beta has to be > or == 0
    Gsetup <- matrix(rep(1 / (comtraits[, 6] + comtraits[, 7] + comtraits[, 3]), nrow(comtraits)),
                     nrow=nrow(comtraits), ncol=nrow(comtraits))
    Pdens <- matrix(rep(comtraits[, 5], nrow(comtraits)), 
                    nrow = nrow(comtraits), ncol = nrow(comtraits))
    Pcont <- matrix(rep(comtraits[, 5] / sum(comtraits[, 5]), each = nrow(comtraits)), 
                    nrow = nrow(comtraits), ncol = nrow(comtraits))
    Pnew <- matrix(rep(comtraits[, 5] / sum(comtraits[, 5]), nrow(comtraits)), 
                   nrow = nrow(comtraits), ncol = nrow(comtraits))
    P <- Pdens * Pcont # Combine density and proportion contact matrices to get the P matrix
    G <- B.m * Gsetup * Pnew # (New formulation? Matches Allen et al. 2011 better...)
  }
  
  if (mode == "dens"){
  	B.1 <- matrix(rep(comtraits[, 8], nrow(comtraits)), nrow = nrow(comtraits), ncol = nrow(comtraits))
  	B.2 <- matrix(rep(comtraits[, 8], each = nrow(comtraits)), nrow = nrow(comtraits), ncol=nrow(comtraits))
  	B.m = cij * (B.1 + B.2) / 2
  	diag(B.m) <- comtraits[, 8] # Specifies intraspecific transmission terms along the diagonal
  	stopifnot (B.m >= 0) # Beta has to be > or == 0
    Gsetup <- matrix(rep(1 / (comtraits[, 6] + comtraits[, 7] + comtraits[, 3]), nrow(comtraits)), 
                     nrow = nrow(comtraits), ncol = nrow(comtraits), byrow=T)
    Pdens <- matrix(rep(comtraits[, 5], nrow(comtraits)), nrow = nrow(comtraits), ncol = nrow(comtraits))
    G <- B.m * Gsetup * Pdens
  }
  eigen(G)$values[1] # Consistent with package "popbio"'s calculation of dominant eigenvalues
}