# Deterministic & stochastic extinctions plotted together
source("GPool.R")
source("CommunityR0.R")
source("ComDis.R")
source("ComPerm2.R")

niter = 1000
k <- .5
w <- 1.5

### Generate data with deterministic extinctions ###
determ.d <- array(dim=c(1, 4)) 
colnames(determ.d) <- c("permutation", "inversions", "delta.Ro", "percent.dR0")

for (i in 1:niter){
  pool1 <- GPool(globmeth="allom", R0ii="TG", Nglobal=5, a=3, k=k, w=w)
  #### freq/free
  test1 <- ComPerm2(pool1, mode="freq", exmeth="deterministic", kWeightPenalty=1, 
                    iter=1, Kmeth="free", cij=.1)
  sub1 <- subset(test1$data, !is.na(delta.Ro))
  sub1$percent.dR0 <- with(sub1, Ro / (Ro - delta.Ro))
  determ.d <- rbind(determ.d, sub1[, c(9, 10, 7, 11)])
  #### dens/free
  test2 <- ComPerm2(pool1, mode="dens", exmeth="deterministic", kWeightPenalty=1, 
                    iter=1, Kmeth="free", cij=.1)
  sub2 <- subset(test2$data, !is.na(delta.Ro))
  sub2$percent.dR0 <- with(sub2, Ro / (Ro - delta.Ro))
  determ.d <- rbind(determ.d, sub2[, c(9, 10, 7, 11)])
  #### freq/fixed
  test3 <- ComPerm2(pool1, mode="freq", exmeth="deterministic", kWeightPenalty=1, 
                    iter=1, Kmeth="fixed", cij=.1)
  sub3 <- subset(test3$data, !is.na(delta.Ro))
  sub3$percent.dR0 <- with(sub3, Ro / (Ro - delta.Ro))
  determ.d <- rbind(determ.d, sub3[, c(9, 10, 7, 11)])
  #### dens/fixed
  test4 <- ComPerm2(pool1, mode="dens", exmeth="deterministic", kWeightPenalty=1, 
                    iter=1, Kmeth="fixed", cij=.1)
  sub4 <- subset(test4$data, !is.na(delta.Ro))
  sub4$percent.dR0 <- with(sub4, Ro / (Ro - delta.Ro))
  determ.d <- rbind(determ.d, sub4[, c(9, 10, 7, 11)])
  print(i)
}

determ.d <- data.frame(determ.d)
determ.d <- determ.d[-1, ]
scen <- rep(c("freq/free", "dens/free", "freq/fixed", "dens/fixed"), each=nrow(sub1))
scen <- rep(scen, niter)
determ.d$scen <- scen

### generate data with stochastic extinctions ###
stoch.d <- array(dim=c(1, 4)) # 4 scenarios X 100 communities X number of rows for each subdataframe
colnames(stoch.d) <- c("permutation", "inversions", "delta.Ro", "percent.dR0")
kWP <- -1 # weight penalty
P.ex <- kWP*pool1$global.pool$weight #1/pool1$global.pool$weight^kWP
par(mai=c(rep(1.3, 4)))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
P.ex <- range01(P.ex)
par(mfrow=c(1,1))
plot(pool1$global.pool$weight, 1-P.ex, type="l", ylab="Pr(extirpation)")

for (i in 1:niter){
  pool1 <- GPool(globmeth="allom", R0ii="TG", Nglobal=5, a=3, k=k, w=w)
  #### freq/free
  test1 <- ComPerm2(pool1, mode="freq", exmeth="stochastic", kWeightPenalty=1, 
                    iter=1, Kmeth="free", cij=.1)
  sub1 <- subset(test1$data, !is.na(delta.Ro))
  sub1$percent.dR0 <- with(sub1, Ro / (Ro - delta.Ro))
  stoch.d <- rbind(stoch.d, sub1[, c(9, 10, 7, 11)])
  #### dens/free
  test2 <- ComPerm2(pool1, mode="dens", exmeth="stochastic", kWeightPenalty=1, 
                    iter=1, Kmeth="free", cij=.1)
  sub2 <- subset(test2$data, !is.na(delta.Ro))
  sub2$percent.dR0 <- with(sub2, Ro / (Ro - delta.Ro))
  stoch.d <- rbind(stoch.d, sub2[, c(9, 10, 7, 11)])
  #### freq/fixed
  test3 <- ComPerm2(pool1, mode="freq", exmeth="stochastic", kWeightPenalty=1, 
                    iter=1, Kmeth="fixed", cij=.1)
  sub3 <- subset(test3$data, !is.na(delta.Ro))
  sub3$percent.dR0 <- with(sub3, Ro / (Ro - delta.Ro))
  stoch.d <- rbind(stoch.d, sub3[, c(9, 10, 7, 11)])
  #### dens/fixed
  test4 <- ComPerm2(pool1, mode="dens", exmeth="stochastic", kWeightPenalty=1, 
                    iter=1, Kmeth="fixed", cij=.1)
  sub4 <- subset(test4$data, !is.na(delta.Ro))
  sub4$percent.dR0 <- with(sub4, Ro / (Ro - delta.Ro))
  stoch.d <- rbind(stoch.d, sub4[, c(9, 10, 7, 11)])
  print(i)
}


### generate data with random extinctions ###
rand.d <- array(dim=c(1, 4)) # 4 scenarios X 100 communities X number of rows for each subdataframe
colnames(rand.d) <- c("permutation", "inversions", "delta.Ro", "percent.dR0")
kWP <- 0 # weight penalty
P.ex <- kWP*pool1$global.pool$weight #1/pool1$global.pool$weight^kWP
par(mai=c(rep(1.3, 4)))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
P.ex <- range01(P.ex)
par(mfrow=c(1,1))
plot(pool1$global.pool$weight, 1-P.ex, type="l", ylab="Pr(extirpation)")

for (i in 1:niter){
  pool1 <- GPool(globmeth="allom", R0ii="TG", Nglobal=5, a=3, k=k, w=w)
  #### freq/free
  test1 <- ComPerm2(pool1, mode="freq", exmeth="stochastic", kWeightPenalty=0, 
                    iter=1, Kmeth="free", cij=.1)
  sub1 <- subset(test1$data, !is.na(delta.Ro))
  sub1$percent.dR0 <- with(sub1, Ro / (Ro - delta.Ro))
  rand.d <- rbind(rand.d, sub1[, c(9, 10, 7, 11)])
  #### dens/free
  test2 <- ComPerm2(pool1, mode="dens", exmeth="stochastic", kWeightPenalty=0, 
                    iter=1, Kmeth="free", cij=.1)
  sub2 <- subset(test2$data, !is.na(delta.Ro))
  sub2$percent.dR0 <- with(sub2, Ro / (Ro - delta.Ro))
  rand.d <- rbind(rand.d, sub2[, c(9, 10, 7, 11)])
  #### freq/fixed
  test3 <- ComPerm2(pool1, mode="freq", exmeth="stochastic", kWeightPenalty=0, 
                    iter=1, Kmeth="fixed", cij=.1)
  sub3 <- subset(test3$data, !is.na(delta.Ro))
  sub3$percent.dR0 <- with(sub3, Ro / (Ro - delta.Ro))
  rand.d <- rbind(rand.d, sub3[, c(9, 10, 7, 11)])
  #### dens/fixed
  test4 <- ComPerm2(pool1, mode="dens", exmeth="stochastic", kWeightPenalty=0, 
                    iter=1, Kmeth="fixed", cij=.1)
  sub4 <- subset(test4$data, !is.na(delta.Ro))
  sub4$percent.dR0 <- with(sub4, Ro / (Ro - delta.Ro))
  rand.d <- rbind(rand.d, sub4[, c(9, 10, 7, 11)])
  print(i)
}

rand.d <- data.frame(rand.d)
rand.d <- rand.d[-1, ]
rand.d$scen <- scen
rand.d$percent.dR0 <- rand.d$percent.dR0*100
rand.d$percent.dR0 <- rand.d$percent.dR0 - 100



stoch.d <- data.frame(stoch.d)
stoch.d <- stoch.d[-1, ]
stoch.d$scen <- scen
stoch.d$percent.dR0 <- stoch.d$percent.dR0*100
determ.d$percent.dR0 <- determ.d$percent.dR0*100
stoch.d$percent.dR0 <- stoch.d$percent.dR0 - 100
determ.d$percent.dR0 <- determ.d$percent.dR0 - 100
stoch.d$lpR0 <- log((stoch.d$percent.dR0 + 100))
determ.d$lpR0 <- log((determ.d$percent.dR0 + 100))

hist(stoch.d$percent.dR0 + 100)
### Plot results side by side for easy comparison ###
#ylims <- c(min(stoch.d$delta.Ro), max(stoch.d$delta.Ro))
#plot.v <- function(scenario){
#  scenario <<- scenario
#  require(vioplot)
#  red <- "firebrick1"
#  blue <-"lightskyblue"
#  botleft <- c(-10, -4)
#  topright <- c(100, 4)
#  xmax <- length(unique(determ.d$inversions)) * 2
#  xlimits <- seq(.5, xmax)
#  plot(x=xlimits, y=rep(ylims, length.out=length(xlimits)), 
#       type="n", ann=FALSE, axes=F)
#  lim <- par("usr")
#  rect(botleft[1], botleft[2], topright[1], topright[2], border=blue, col=blue)
#  rect(-10, 0, 100, -4, border=red, col=red)
#  determ.x <- which(1:(length(unique(determ.d$inversions))*2) %% 2 == 1) # every 2nd value
#  determ.x <<- determ.x
#  stoch.x <- determ.x + 1
#  for (i in unique(determ.d$inversions)){
#    vioplot(subset(determ.d, scen == scenario & inversions == i)$delta.Ro,
#            col="white", colMed="lightgrey", ylim=ylims,
#            at=determ.x[which(unique(determ.d$inversions) == i)], add=T)
#    vioplot(subset(stoch.d, scen == scenario & inversions == i)$delta.Ro,
#            col="yellow", colMed="lightgrey", ylim=ylims,
#            at=stoch.x[which(unique(stoch.d$inversions) == i)], add=T)
#  }
#}


# paneled violin plots
#par(mfrow=c(2,2), mai=rep(c(.5, .5), 2))
#plot.v("freq/free")
#title(scenario)
#axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))
#plot.v("dens/free")
#title(scenario)
#axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))
#plot.v("freq/fixed")
#title(scenario)
#axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))
#plot.v("dens/fixed")
#title(scenario)
#axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))


### With percentages
ylarge <- max(c(stoch.d$percent.dR0, determ.d$percent.dR0))
ysmall <- min(c(stoch.d$percent.dR0, determ.d$percent.dR0))
ylims <- c(ysmall, ylarge)

plot.v <- function(scenario, ylab=NULL){
  scenario <<- scenario
  require(vioplot)
  red <- "firebrick1"
  blue <-"lightskyblue"
  botleft <- c(-10, -4)
  topright <- c(100, 40000)
  xmax <- length(unique(determ.d$inversions)) * 2 + .5
  xlimits <- seq(.5, xmax)
  plot(x=xlimits, y=rep(ylims, length.out=length(xlimits)), 
       type="n", ann=FALSE, axes=F, ylab=ylab)
  lim <- par("usr")
  rect(botleft[1], botleft[2], topright[1], topright[2], border=blue, col=blue)
  rect(-10, 0, 100, -4000000, border=red, col=red)
  determ.x <- which(1:(length(unique(determ.d$inversions))*2) %% 2 == 1) # every 2nd value
  determ.x <<- determ.x
  stoch.x <- determ.x + 1
  for (i in unique(determ.d$inversions)){
    vioplot(subset(determ.d, scen == scenario & inversions == i)$percent.dR0,
            col="white", colMed="lightgrey", ylim=ylims,
            at=determ.x[which(unique(determ.d$inversions) == i)], add=T)
    vioplot(subset(stoch.d, scen == scenario & inversions == i)$percent.dR0,
            col="yellow", colMed="lightgrey", ylim=ylims,
            at=stoch.x[which(unique(stoch.d$inversions) == i)], add=T)
  }
}

# paneled violin plots
par(mfrow=c(2,2), mai=rep(c(.5, .5), 2))
plot.v("freq/free")
title(scenario)
axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))
plot.v("dens/free")
title(scenario)
axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))
plot.v("freq/fixed")
title(scenario)
axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))
plot.v("dens/fixed")
title(scenario)
axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))


# Adjusting for spacing
in.marg <- .05
nf <- layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=T))
layout.show(nf)
marg.size <- 6
par(oma=c(marg.size-2, marg.size, 2, 1))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v("freq/free")
axis(side=2)
mtext(expression(paste("Extirpation-induced %", Delta, R[0])), side=2, line=1.5, adj=0.0, cex=.9, at=c(.6), outer=TRUE)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v("dens/free")
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v("freq/fixed")
axis(side=2)
mtext(expression(paste("Extirpation-induced %", Delta, R[0])), side=2, line=1.5, adj=0.0, cex=.9, at=c(.1), outer=TRUE)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v("dens/fixed")

# Adding miscellaneous labels
# Left side
mtext("Variable density", side=2, line=4, cex=1.1, at=.75, outer=T)
mtext("Fixed density", side=2, line=4, cex=1.1, at=.25, outer=T)
mtext("Frequency-dependent transmission", side=3, line=0, outer=T, at=.25, cex=1.1)
mtext("Density-dependent transmission", side=3, line=0, outer=T, at=.75, cex=1.1)
mtext("Relationship between host competence and body size",
      side=1, line=2.5, outer=T, at=.5, cex=1.4)
mtext("Negative", side=1, line=0, outer=T, at=c(.05, .55))
mtext("Positive", side=1, line=0, outer=T, at=c(.45, .95))
mtext("None", side=1, line=0, outer=T, at=c(.25, .75))

### Plotting deterministic case over the stochastic case
plot.v2 <- function(scenario, ylab=NULL){
  scenario <<- scenario
  require(vioplot)
  red <- "firebrick1"
  blue <-"lightskyblue"
  botleft <- c(-10, -4)
  topright <- c(100, 40000)
  xmax <- length(unique(determ.d$inversions)) + .5
  xlimits <- seq(.5, xmax)
  plot(x=xlimits, y=rep(ylims, length.out=length(xlimits)), 
       type="n", ann=FALSE, axes=F, ylab=ylab)
  lim <- par("usr")
  rect(botleft[1], botleft[2], topright[1], topright[2], border=blue, col=blue)
  rect(-10, 0, 100, -4000000, border=red, col=red)
  determ.x <<- 1:length(unique(determ.d$inversions))
  for (i in unique(determ.d$inversions)){
    vioplot(subset(stoch.d, scen == scenario & inversions == i)$percent.dR0,
            col="yellow", colMed="lightgrey", ylim=ylims,
            at=determ.x[which(unique(determ.d$inversions) == i)], add=T)
    vioplot(subset(determ.d, scen == scenario & inversions == i)$percent.dR0,
            col="white", colMed="lightgrey", ylim=ylims,
            at=determ.x[which(unique(determ.d$inversions) == i)], add=T)
  }
}


# paneled violin plots
par(mfrow=c(2,2), mai=rep(c(.5, .5), 2))
plot.v2("freq/free")
title(scenario)
axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))
plot.v2("dens/free")
title(scenario)
axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))
plot.v2("freq/fixed")
title(scenario)
axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))
plot.v2("dens/fixed")
title(scenario)
axis(side=1, at=determ.x, labels=rev(unique(determ.d$inversions)))


# Adjusting for spacing
in.marg <- .05
nf <- layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=T))
layout.show(nf)
marg.size <- 6
par(oma=c(marg.size-2, marg.size, 2, 1))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v2("freq/free")
axis(side=2)
mtext(expression(paste("Extirpation-induced %", Delta, R[0])), side=2, line=1.5, adj=0.0, cex=.9, at=c(.6), outer=TRUE)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v2("dens/free")
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v2("freq/fixed")
axis(side=2)
mtext(expression(paste("Extirpation-induced %", Delta, R[0])), side=2, line=1.5, adj=0.0, cex=.9, at=c(.1), outer=TRUE)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v2("dens/fixed")

# Adding miscellaneous labels
# Left side
mtext("Variable density", side=2, line=4, cex=1.1, at=.75, outer=T)
mtext("Fixed density", side=2, line=4, cex=1.1, at=.25, outer=T)
mtext("Frequency-dependent transmission", side=3, line=0, outer=T, at=.25, cex=1.1)
mtext("Density-dependent transmission", side=3, line=0, outer=T, at=.75, cex=1.1)
mtext("Relationship between host competence and body size",
      side=1, line=2.5, outer=T, at=.5, cex=1.4)
mtext("Negative", side=1, line=0, outer=T, at=c(.05, .55))
mtext("Positive", side=1, line=0, outer=T, at=c(.45, .95))
mtext("None", side=1, line=0, outer=T, at=c(.25, .75))


### Plotting the distributions side by side
source("~/vioplot2.R")
### Plotting deterministic case over the stochastic case
plot.v3 <- function(scenario, ylab=NULL){
  scenario <<- scenario
  require(vioplot)
  red <- "firebrick1"
  blue <-"lightskyblue"
  botleft <- c(-10, -4)
  topright <- c(100, 40000)
  xmax <- length(unique(determ.d$inversions)) + .5
  xlimits <- seq(.5, xmax)
  plot(x=xlimits, y=rep(ylims, length.out=length(xlimits)), 
       type="n", ann=FALSE, axes=F, ylab=ylab)
  lim <- par("usr")
  rect(botleft[1], botleft[2], topright[1], topright[2], border=blue, col=blue)
  rect(-10, 0, 100, -4000000, border=red, col=red)
  determ.x <<- 1:length(unique(determ.d$inversions))
  for (i in unique(determ.d$inversions)){
    vioplot2(subset(stoch.d, scen == scenario & inversions == i)$percent.dR0,
            col="yellow", colMed="lightgrey", ylim=ylims,
            at=determ.x[which(unique(determ.d$inversions) == i)], 
            side="right", add=T)
    vioplot2(subset(determ.d, scen == scenario & inversions == i)$percent.dR0,
            col="white", colMed="lightgrey",
            at=determ.x[which(unique(determ.d$inversions) == i)], 
            side="left", add=T)
  }
}


# Adjusting for spacing
in.marg <- .05
nf <- layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=T))
layout.show(nf)
marg.size <- 6
par(oma=c(marg.size-2, marg.size, 2, 1))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("freq/free")
axis(side=2)
#mtext(expression(paste("Extirpation-induced %", Delta, R[0])), side=2, line=1.7, adj=0.0, cex=.9, at=c(.6), outer=TRUE)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("dens/free")
legend(x=4, y=300, legend=c("Deterministic extirpations", "Stochastic extirpations"),
       fill=c("white", "yellow"))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("freq/fixed")
axis(side=2)
#mtext(expression(paste("Extirpation-induced %", Delta, R[0])), side=2, line=1.7, adj=0.0, cex=.9, at=c(.1), outer=TRUE)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("dens/fixed")

# Adding miscellaneous labels
# Left side
mtext("Variable density", side=2, line=4, cex=1.1, at=.75, outer=T)
mtext("Fixed density", side=2, line=4, cex=1.1, at=.25, outer=T)
mtext("Frequency-dependent transmission", side=3, line=0, outer=T, at=.25, cex=1.1)
mtext("Density-dependent transmission", side=3, line=0, outer=T, at=.75, cex=1.1)
mtext("Relationship between host competence and body size",
      side=1, line=2.5, outer=T, at=.5, cex=1.4)
mtext("Negative", side=1, line=0, outer=T, at=c(.05, .55))
mtext("Positive", side=1, line=0, outer=T, at=c(.45, .95))
mtext("None", side=1, line=0, outer=T, at=c(.25, .75))
mtext(expression(paste("Extirpation-induced %", Delta, R[0])), side=2, line=1.7, adj=0.0, cex=.9, at=c(.35), outer=TRUE)
