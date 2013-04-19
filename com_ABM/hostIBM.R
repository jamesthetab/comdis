hostIBM <- function(A=100, N=100, ERmin=-30, ERmax=30, Emin=-30, Emax=30, 
                    sig=5, timesteps=1000, pM=.1, pR=1, R=2, I=0.1){
  E <- runif(A, Emin, Emax) # habitat patch environment type
  mu.i <- runif(N, ERmin, ERmax) # optimum environment
  sigma <- rep(sig, N) # runif(N, 1, 10) # # niche width
  Z <- array(dim=c(N)) # normalization constant
  
  # now calculating all elements of Z coefficient 
  # (ensures all species have equal Pr(establishment) in regional pool)
  for (i in 1:N){
    integrand <- function(E) {
      exp(-((E - mu.i[i]) ^ 2) / (2 * sigma[i] ^ 2))
    }
    res <- integrate(integrand, lower=ERmin, upper=ERmax)
    Z[i] <- 1 / res$value
  }
  
  # probability of establishment
  Pcol <- array(dim=c(A, N))
  for (i in 1:A){
    for (j in 1:N){
      Pcol[i, j] <- Z[j] * exp(-((E[i] - mu.i[j]) ^ 2) / (2*sigma[j] ^ 2))
    }
  }

  # store niche data
  species <- rep(1:N, each=A)
  E <- rep(E, N)
  Pr.estab <- c(Pcol)
  niche.d <- data.frame(species, E, Pr.estab)
  niche.d <- niche.d[with(niche.d, order(species, E)),]
  
  #for (i in 1:N){
  #  tempd <- data.frame(Pcol = Pcol[,i], E=E)
  #  tempd <- tempd[order(E),]
  #  Pr.estab <- c(Pr.estab, tempd$Pcol)
  #}
  #niche.d <- data.frame(species, E, Pr.estab)
  
  # initialize output objects
  state <- array(0, dim=c(timesteps, A, N))
  richness <- rep(NA, timesteps)
  richness[1] <- 0
  p.occ <- rep(NA, timesteps)
  p.occ[1] <- 0
  
  for (t in 2:timesteps){
    state[t,,] <- state[t-1,,]
    
    ## DEATHS ##
    deaths <- array(rbinom(A*N, 1, c(state[t,,])*pM), dim=c(A, N))
    state[t,,] <- state[t,,] - deaths
    
    ## BIRTHS ##
    pot.fecundity <- array(rpois(A * N, lambda = c(state[t,,] * R)), dim=c(A, N)) # potential number of offspring
    repro <- array(rbinom(A*N, 1, pR), dim=c(A, N)) # whether reproduction actually occurs
    fecundity <- repro * pot.fecundity
    sum.fec <- apply(fecundity, 2, sum) # number of offspring per species
    
    ## OFFSPRING COLONIZE EMPTY SITES ##
    occupancy <- apply(state[t,,], 1, max)
    if (sum(occupancy) < A & sum(sum.fec) > 0){ # if empty sites & new offspring
      empty.sites <- which(occupancy == 0)
      # assign one empty site randomly to each individual offspring
      if (length(empty.sites) == 1){ # if there's only one empty site
        col.sites <- rep(empty.sites, sum(sum.fec)) # all offspring attempt to colonize
      }else{
        col.sites <- sample(empty.sites, sum(sum.fec), replace=T) # else, randomly choose
      }
      col.spec <- rep(1:N, times=sum.fec) # how many of each species colonizing across sites
      colonizing.offspring <- array(0, dim=c(A, N)) # how many of each species colonizing each site
      for(i in 1:length(col.sites)){
        colonizing.offspring[col.sites[i], col.spec[i]] <- colonizing.offspring[col.sites[i], col.spec[i]] + 1 
      }
      # which colonizing offspring can actually establish?
      binom.mat <- ifelse(colonizing.offspring > 0, 1, 0)
      colonists <- array(rbinom(n = A * N, 
                                size = c(colonizing.offspring),
                                prob = c(binom.mat * Pcol)),
                         dim=c(A, N))
      
      # are there colonization conflicts (> 1 individual trying to colonize each site?)
      attempting <- apply(colonists, 1, sum) # number indiv. attempting to colonize each site
      if (any(attempting > 1)){
        # resolve colonization conflicts
        conflicts <- which(attempting > 1) # which sites have conflicts
        for (k in conflicts){ # for each conflict
          # how many of individuals of each species are attempting to simultaneously colonize?
          n.attempting <- rep(1:N, times = colonists[k,])
          # randomly select one successful from those attempting
          successful <- sample(n.attempting, size=1)
          new.row <- rep(0, length.out=N)
          new.row[successful] <- 1
          colonists[k,] <- new.row # individual becomes the only colonist
        }
      }
      # add successful colonists
      state[t,,] <- state[t,,] + colonists
    } 
    # end of offspring colonization
    
    ## IMMIGRANTS COLONIZE EMPTY SITES ##
    occupancy <- apply(state[t,,], 1, max)
    if(sum(occupancy) < A){
      empty.sites <- which(occupancy == 0)
      # which species immigrate to each site?
      immigration <- array(rbinom(length(empty.sites)*N,
                                  1, I), dim=c(length(empty.sites), N))
      # which immigrants establish?
      Pest <- immigration * Pcol[empty.sites, ]
      establishment <- array(rbinom(length(Pest), 1, c(Pest)),
                             dim=c(length(empty.sites), N))
      
      # resolve conflicts arising from simultaneous colonization
      col.attempts <- apply(establishment, 1, sum)
      if (any(col.attempts > 1)){ # if individuals are trying to simultaneously colonize
        conflicts <- which(col.attempts > 1) # which empty sites have conflicts
        for (k in conflicts){ # for each conflict
          # how many of individuals of each species are attempting to simultaneously colonize?
          attempting <- rep(1:N, times = establishment[k,])
          # successful individual randomly selected from those attempting
          successful <- sample(attempting, size=1)
          new.row <- rep(0, length.out=N)
          new.row[successful] <- 1
          establishment[k,] <- new.row
        }
      }
      # add establishing immigrants
      state[t, empty.sites, ] <- state[t, empty.sites,] + establishment
    }
    richness[t] <- sum(apply(state[t,,], 2, max))
    p.occ[t] <- length(which(state[t,,] == 1)) / A
  }
  # check to confirm there are not multiple indivduals in any sites
  stopifnot(all(apply(state[t,,], 1, sum) < 2))
  return(list(richness=richness, p.occ = p.occ, state=state, niche.d=niche.d))
}