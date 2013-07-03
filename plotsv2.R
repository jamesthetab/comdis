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
start <- seq(1, nrow(determ.d), 3)
mean.dR0 <- rep(NA, nrow(determ.d)/3)
sd.dR0 <- rep(NA, nrow(determ.d)/3)
perm.dR0 <- rep(NA, nrow(determ.d)/3)
inv.dR0 <- rep(NA, nrow(determ.d)/3)
scen.dR0 <- rep(NA, nrow(determ.d)/3)
for (i in 1:length(mean.dR0)){
  indx <- start[i]:(start[i]+2)
  mean.dR0[i] <- mean(determ.d$delta.Ro[indx])
  sd.dR0[i] <- sd(determ.d$delta.Ro[indx])
  perm.dR0[i] <- determ.d$permutation[start[i]]
  inv.dR0[i] <- determ.d$inversions[start[i]]
  scen.dR0[i] <- determ.d$scen[start[i]]
}
boxplot(mean.dR0 ~ inv.dR0)
plot(density(mean.dR0[which(inv.dR0 == 4) %in% which(scen.dR0 == "freq/free")], bw=.05))
lines(density(mean.dR0[which(inv.dR0 == 5) %in% which(scen.dR0 == "freq/free")], bw=.05), col="blue")
determ.d2 <- data.frame(mean.dR0=mean.dR0,
                        sd.dR0=sd.dR0,
                        perm.dR0=perm.dR0,
                        inv.dR0=inv.dR0,
                        scen.dR0=scen.dR0)

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

mean.dR0 <- rep(NA, nrow(determ.d)/3)
sd.dR0 <- rep(NA, nrow(determ.d)/3)
perm.dR0 <- rep(NA, nrow(determ.d)/3)
inv.dR0 <- rep(NA, nrow(determ.d)/3)
scen.dR0 <- rep(NA, nrow(determ.d)/3)
for (i in 1:length(mean.dR0)){
  indx <- start[i]:(start[i]+2)
  mean.dR0[i] <- mean(stoch.d$delta.Ro[indx])
  sd.dR0[i] <- sd(stoch.d$delta.Ro[indx])
  perm.dR0[i] <- stoch.d$permutation[start[i]]
  inv.dR0[i] <- stoch.d$inversions[start[i]]
  scen.dR0[i] <- stoch.d$scen[start[i]]
}

stoch.d2 <- data.frame(mean.dR0=mean.dR0,
                        sd.dR0=sd.dR0,
                        perm.dR0=perm.dR0,
                        inv.dR0=inv.dR0,
                        scen.dR0=scen.dR0)

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

mean.dR0 <- rep(NA, nrow(determ.d)/3)
sd.dR0 <- rep(NA, nrow(determ.d)/3)
perm.dR0 <- rep(NA, nrow(determ.d)/3)
inv.dR0 <- rep(NA, nrow(determ.d)/3)
scen.dR0 <- rep(NA, nrow(determ.d)/3)
for (i in 1:length(mean.dR0)){
  indx <- start[i]:(start[i]+2)
  mean.dR0[i] <- mean(rand.d$delta.Ro[indx])
  sd.dR0[i] <- sd(rand.d$delta.Ro[indx])
  perm.dR0[i] <- rand.d$permutation[start[i]]
  inv.dR0[i] <- rand.d$inversions[start[i]]
  scen.dR0[i] <- rand.d$scen[start[i]]
}

rand.d2 <- data.frame(mean.dR0=mean.dR0,
                       sd.dR0=sd.dR0,
                       perm.dR0=perm.dR0,
                       inv.dR0=inv.dR0,
                       scen.dR0=scen.dR0)



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



#
## Visualization
#####################################################
### Plot results side by side for easy comparison ###
#####################################################
# Deterministic vs. stochastic
source("~/vioplot2.R")
ylims <- c(min(stoch.d$delta.Ro), max(stoch.d$delta.Ro))
plot.v3 <- function(scenario, ylab=NULL){
  scenario <<- scenario
  require(vioplot)
  red <- "firebrick1"
  blue <-"lightskyblue"
  botleft <- c(-1, -4)
  topright <- c(12, 0)
  xmax <- length(unique(determ.d$inversions)) + .5
  xlimits <- seq(.5, xmax)
  plot(x=xlimits, y=rep(ylims, length.out=length(xlimits)), 
       type="n", ann=FALSE, axes=F, ylab=ylab)
  lim <- par("usr")
  rect(botleft[1], botleft[2], topright[1], topright[2], border=red, col=red)
  rect(-1, 0, 12, 4, border=blue, col=blue)
  determ.x <<- 1:length(unique(determ.d$inversions))
  for (i in unique(determ.d$inversions)){
    vioplot2(subset(stoch.d, scen == scenario & inversions == i)$delta.Ro,
             col="yellow", colMed="lightgrey", ylim=ylims,
             at=determ.x[which(unique(determ.d$inversions) == i)], 
             side="right", add=T)
    vioplot2(subset(determ.d, scen == scenario & inversions == i)$delta.Ro,
             col="white", colMed="lightgrey",
             at=determ.x[which(unique(determ.d$inversions) == i)], 
             side="left", add=T)
  }
}

cexadj <- .8
in.marg <- .05
nf <- layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=T))
layout.show(nf)
marg.size <- 6
par(oma=c(marg.size-2, marg.size, 4, 1))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("freq/free")
axis(side=2, cex.axis=cexadj)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("dens/free")
legend(x=6, y=2.5, legend=c("Deterministic extirpations", "Stochastic extirpations"),
       fill=c("white", "yellow"), cex=cexadj)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("freq/fixed")
axis(side=2)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("dens/fixed")
mtext("Variable density", side=2, line=4.5, cex=cexadj*1.1, at=.75, outer=T)
mtext("Fixed density", side=2, line=4.5, cex=cexadj*1.1, at=.25, outer=T)
mtext("Frequency-dependent transmission", side=3, line=0, outer=T, at=.25, cex=cexadj*1.1)
mtext("Density-dependent transmission", side=3, line=0, outer=T, at=.75, cex=cexadj*1.1)
mtext("Relationship between host competence and body size",
      side=1, line=2.5, outer=T, at=.5, cex=cexadj*1.4)
mtext("Negative", side=1, line=0, outer=T, at=c(.05, .55), cex=cexadj*.9)
mtext("Positive", side=1, line=0, outer=T, at=c(.45, .95), cex=cexadj*.9)
mtext("None", side=1, line=0, outer=T, at=c(.25, .75), cex=cexadj*.9)
mtext(expression(paste("Extirpation-induced ", Delta, R[0])), side=2, line=2, adj=0.0, cex=cexadj*.9, at=c(.35), outer=TRUE)
mtext("Figure 2: Effects of extirpation on disease risk",
      side=3, line=2.5, outer=T, at=.5, cex=cexadj*1.4)
# Save as 10 X 6 pdf


### With percentages ###
ylarge <- max(c(stoch.d$percent.dR0, determ.d$percent.dR0))
ysmall <- min(c(stoch.d$percent.dR0, determ.d$percent.dR0))
ylims <- c(ysmall, ylarge)

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
  rect(botleft[1], botleft[2], topright[1], topright[2], border=red, col=red)
  rect(-10, 0, 100, -4000000, border=blue, col=blue)
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

in.marg <- .05
nf <- layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=T))
layout.show(nf)
marg.size <- 6
par(oma=c(marg.size-2, marg.size, 2, 1))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("freq/free")
axis(side=2)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("dens/free")
legend(x=5, y=300, legend=c("Deterministic extirpations", "Stochastic extirpations"),
       fill=c("white", "yellow"))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("freq/fixed")
axis(side=2)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("dens/fixed")
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






#--------------------------#
# Deterministic vs. random #
#--------------------------#
ylims <- c(min(rand.d$delta.Ro), max(rand.d$delta.Ro))
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
  rect(botleft[1], botleft[2], topright[1], topright[2], border=red, col=red)
  rect(-10, 0, 100, -4000000, border=blue, col=blue)
  determ.x <<- 1:length(unique(determ.d$inversions))
  for (i in unique(determ.d$inversions)){
    vioplot2(subset(rand.d, scen == scenario & inversions == i)$delta.Ro,
             col="yellow", colMed="lightgrey", ylim=ylims,
             at=determ.x[which(unique(determ.d$inversions) == i)], 
             side="right", add=T)
    vioplot2(subset(determ.d, scen == scenario & inversions == i)$delta.Ro,
             col="white", colMed="lightgrey",
             at=determ.x[which(unique(determ.d$inversions) == i)], 
             side="left", add=T)
  }
}

in.marg <- .05
nf <- layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=T))
layout.show(nf)
marg.size <- 6
par(oma=c(marg.size-2, marg.size, 2, 1))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("freq/free")
axis(side=2)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("dens/free")
legend(x=5, y=5, legend=c("Deterministic extirpations", "Random extirpations"),
       fill=c("white", "yellow"))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("freq/fixed")
axis(side=2)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("dens/fixed")
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




# Plots based on the distribution of expected values (should be more normal)
# deterministic vs. stochastic
ylims <- c(min(stoch.d2$mean.dR0), max(stoch.d2$mean.dR0))
plot.v4 <- function(scenario, ylab=NULL){
  scenario <<- scenario
  require(beanplot)
  red <- "firebrick1"
  blue <-"lightskyblue"
  botleft <- c(-10, -4)
  topright <- c(100, 40000)
  xmax <- length(unique(determ.d$inversions)) + .5
  xlimits <- seq(.5, xmax)
  plot(x=xlimits, y=rep(ylims, length.out=length(xlimits)), 
       type="n", ann=FALSE, axes=F, ylab=ylab)
  lim <- par("usr")
  rect(botleft[1], botleft[2], topright[1], topright[2], border=red, col=red)
  rect(-10, 0, 100, -4000000, border=blue, col=blue)
  determ.x <<- 1:length(unique(determ.d$inversions))
  for (i in unique(determ.d$inversions)){
    vioplot2(subset(rand.d2, scen.dR0 == scenario & inv.dR0 == i)$mean.dR0,
             col="yellow", colMed="lightgrey", ylim=ylims,
             at=determ.x[which(unique(determ.d$inversions) == i)], 
             side="right", add=T)
    vioplot2(subset(determ.d2, scen.dR0 == scenario & inv.dR0 == i)$mean.dR0,
             col="white", colMed="lightgrey", ylim=ylims,
             at=determ.x[which(unique(determ.d$inversions) == i)], 
             side="left", add=T)
  }
}

in.marg <- .05
nf <- layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=T))
layout.show(nf)
marg.size <- 6
par(oma=c(marg.size-2, marg.size, 2, 1))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v4("freq/free")
axis(side=2)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v4("dens/free")
legend(x=5, y=5, legend=c("Deterministic extirpations", "Random extirpations"),
       fill=c("white", "yellow"))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v4("freq/fixed")
axis(side=2)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v4("dens/fixed")
mtext("Variable density", side=2, line=4, cex=1.1, at=.75, outer=T)
mtext("Fixed density", side=2, line=4, cex=1.1, at=.25, outer=T)
mtext("Frequency-dependent transmission", side=3, line=0, outer=T, at=.25, cex=1.1)
mtext("Density-dependent transmission", side=3, line=0, outer=T, at=.75, cex=1.1)
mtext("Relationship between host competence and body size",
      side=1, line=2.5, outer=T, at=.5, cex=1.4)
mtext("Negative", side=1, line=0, outer=T, at=c(.05, .55))
mtext("Positive", side=1, line=0, outer=T, at=c(.45, .95))
mtext("None", side=1, line=0, outer=T, at=c(.25, .75))
mtext(expression(paste("Expected ", Delta, R[0])), side=2, line=1.7, adj=0.0, cex=.9, at=c(.35), outer=TRUE)

#######################################################################
## Plot showing the basis for competence-extinction risk relationship #
#######################################################################
Nspec <- 12
sd1 <- 1
sd2 <- 2
lhpace <- rnorm(Nspec, 0, 3)
inv.def <- -lhpace + rnorm(Nspec, 0, sd1)
inv.def2 <- -lhpace + rnorm(Nspec, 0, sd2)
compet <- -inv.def + rnorm(Nspec, 0, sd1)
compet2 <- -inv.def2 + rnorm(Nspec, 0, sd2)
ext.risk <- -lhpace + rnorm(Nspec, 0, sd1)
ext.risk2 <- -lhpace + rnorm(Nspec, 0, sd2)
plot(lhpace, inv.def)
plot(inv.def, compet)
plot(lhpace, ext.risk)
plot(ext.risk, compet)
require(ggplot2)
theme(axis.ticks = element_blank(), axis.text.x = element_blank())
errcol <- "red"
sh <- 1
f1d <- data.frame(lhpace, inv.def, compet, ext.risk,
                  inv.def2, compet2, ext.risk2) 
panel1 <- ggplot(f1d, aes(x=lhpace, y=inv.def)) + theme_classic() +
  stat_smooth(aes(x=lhpace, y=inv.def2), method=lm, fill="lightgrey", color=errcol, linetype="dashed") +
  stat_smooth(method=lm, fill="darkgray", color="black") +
  geom_point(aes(x=lhpace, y=inv.def2), shape=sh, size=3, color=errcol) + 
  geom_point(shape=sh, size=3) +
  xlab("Life history pace") + ylab("Investment in defense")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
panel2 <- ggplot(f1d, aes(x=lhpace, y=ext.risk)) + theme_classic() + 
  stat_smooth(method=lm, fill="darkgrey", color="black")+
  stat_smooth(aes(x=lhpace, y=ext.risk2), method=lm, fill="lightgrey", color=errcol, linetype="dashed") +
  geom_point(aes(x=lhpace, y=ext.risk2), shape=sh, size=3, color=errcol) + 
  geom_point(shape=sh, size=3) +
  xlab("Life history pace") + ylab("Extirpation risk")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
panel3 <- ggplot(f1d, aes(x=inv.def, y=compet)) + theme_classic() +
  stat_smooth(method=lm, fill="darkgrey", color="black")+
  stat_smooth(aes(x=inv.def2, y=compet2), method=lm, fill="lightgrey", color=errcol, linetype="dashed") + 
  geom_point(aes(x=inv.def2, y=compet2), shape=sh, size=3, color=errcol) + 
  geom_point(shape=sh, size=3) +
  ylab("Host competence") + xlab("Investment in defense")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
panel4 <- ggplot(f1d, aes(x=ext.risk, y=compet)) + theme_classic() +
  stat_smooth(method=lm, fill="darkgrey", color="black")+
  stat_smooth(aes(x=ext.risk2, y=compet2), method=lm, fill="lightgrey", color=errcol, linetype="dashed") +
  geom_point(aes(x=ext.risk2, y=compet2), shape=sh, size=3, color=errcol) + 
  geom_point(shape=sh, size=3) + 
  xlab("Extirpation risk") + ylab("Host competence")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
require(gridExtra)
grid.arrange(panel1, panel2, panel3, panel4, ncol=2)