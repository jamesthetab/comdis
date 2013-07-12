# Revised Ecology Letters LH paper figures

################
### Figure 1 ###
################
slope <- - 1
Nspec <- 10
sd1 <- .5
sd2 <- 1
while(slope < 0){
  lhpace <- rnorm(Nspec, 0, 3)
  inv.def <- -lhpace + rnorm(Nspec, 0, sd1)
  inv.def2 <- -lhpace + rnorm(Nspec, 0, sd2)
  compet <- -inv.def + rnorm(Nspec, 0, sd1)
  compet2 <- -inv.def2 + rnorm(Nspec, 0, sd2)
  ext.risk <- -lhpace + rnorm(Nspec, 0, sd1)
  ext.risk2 <- -lhpace + rnorm(Nspec, 0, sd2)
  slope <- lm(compet2 ~ ext.risk2)$coefficients[2]
}
require(ggplot2)
errcol <- "red"
sh2 <- 1
sh1 <- 19
f1d <- data.frame(lhpace, inv.def, compet, ext.risk,
                  inv.def2, compet2, ext.risk2) 
panel1 <- ggplot(f1d, aes(x=lhpace, y=inv.def)) + theme_classic() +
  stat_smooth(aes(x=lhpace, y=inv.def2), method=lm, fill="lightgrey", color=errcol, linetype="dashed") +
  stat_smooth(method=lm, fill="darkgray", color="black") +
  geom_point(shape=sh1, size=3) +
  geom_point(aes(x=lhpace, y=inv.def2), shape=sh2, size=3, color=errcol) + 
  xlab("Life history pace") + ylab("Investment in defense")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  ggtitle('A') + theme(plot.title=element_text(hjust=0, vjust=2))
panel2 <- ggplot(f1d, aes(x=lhpace, y=ext.risk)) + theme_classic() + 
  stat_smooth(method=lm, fill="darkgrey", color="black")+
  stat_smooth(aes(x=lhpace, y=ext.risk2), method=lm, fill="lightgrey", color=errcol, linetype="dashed") +
  geom_point(shape=sh1, size=3) +
  geom_point(aes(x=lhpace, y=ext.risk2), shape=sh2, size=3, color=errcol) + 
  xlab("Life history pace") + ylab("Extirpation risk")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  ggtitle('B') + theme(plot.title=element_text(hjust=0, vjust=2))
panel3 <- ggplot(f1d, aes(x=inv.def, y=compet)) + theme_classic() +
  stat_smooth(method=lm, fill="darkgrey", color="black")+
  stat_smooth(aes(x=inv.def2, y=compet2), method=lm, fill="lightgrey", color=errcol, linetype="dashed") + 
  geom_point(shape=sh1, size=3) +
  geom_point(aes(x=inv.def2, y=compet2), shape=sh2, size=3, color=errcol) + 
  ylab("Host competence") + xlab("Investment in defense")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  ggtitle('C') + theme(plot.title=element_text(hjust=0, vjust=2))
panel4 <- ggplot(f1d, aes(x=ext.risk, y=compet)) + theme_classic() +
  stat_smooth(method=lm, fill="darkgrey", color="black")+
  stat_smooth(aes(x=ext.risk2, y=compet2), method=lm, fill="lightgrey", color=errcol, linetype="dashed") +
  geom_point(shape=sh1, size=3) +
  geom_point(aes(x=ext.risk2, y=compet2), shape=sh2, size=3, color=errcol) + 
  xlab("Extirpation risk") + ylab("Host competence")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  ggtitle('D') + theme(plot.title=element_text(hjust=0, vjust=2))
require(gridExtra)
grid.arrange(panel1, panel2, panel3, panel4, ncol=2)









################
### Figure 2 ###
################
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

# Deterministic vs. stochastic
source("~/vioplot2.R")
ylims <- c(min(stoch.d$delta.Ro), max(stoch.d$delta.Ro))
plot.v3 <- function(scenario, ylab=NULL){
  scenario <<- scenario
  require(vioplot)
  red <- "firebrick2"
  blue <- "steelblue"
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
             col="lightgrey", ylim=ylims,
             at=determ.x[which(unique(determ.d$inversions) == i)], 
             side="right", add=T)
    vioplot2(subset(determ.d, scen == scenario & inversions == i)$delta.Ro,
             col="white", 
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
       fill=c("white", "lightgrey"), cex=cexadj)
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
# Save as 11 X 6 pdf





################
### Figure 3 ###
################
pool1 <- GPool(globmeth="allom", R0ii="TG", Nglobal=5, a=3, k=k, w=w)
#27
pool1
#### dens/fixed
cij <- seq(0.01, .99, length.out=100)
y <- rep(cij, each=11) # values of cij, repeated 6 times, once for each inversion
x <- rep(0:10, length(cij)) # inversions
z <- matrix(NA, nrow=length(cij), ncol=11)

for (i in 1:length(cij)){
  d <- ComPerm2(pool1, mode="dens", iter=1, Kmeth="fixed", cij=cij[i], exmeth="deterministic")
  subd <- subset(d$data, !is.na(delta.Ro))
  for (inv in unique(subd$inversions)){
    rows <- which(subd$inversions == inv)
    mean.eff <- mean(subd$delta.Ro[rows], na.rm=T)
    z[i, inv+1] <- mean.eff
  }
}

z2 <- c(t(z))
df1 <- data.frame(x, y, z2)

#### freq/free
z <- matrix(NA, nrow=length(cij), ncol=11)

for (i in 1:length(cij)){
  d <- ComPerm2(pool1, mode="freq", iter=1, Kmeth="free", cij=cij[i], exmeth="deterministic")
  subd <- subset(d$data, !is.na(delta.Ro))
  for (inv in unique(subd$inversions)){
    rows <- which(subd$inversions == inv)
    mean.eff <- mean(subd$delta.Ro[rows], na.rm=T)
    z[i, inv+1] <- mean.eff
  }
}

z2 <- c(t(z))
df2 <- data.frame(x, y, z2)

#### dens/free
z <- matrix(NA, nrow=length(cij), ncol=11)
for (i in 1:length(cij)){
  d <- ComPerm2(pool1, mode="dens", iter=1, Kmeth="free", cij=cij[i], exmeth="deterministic")
  subd <- subset(d$data, !is.na(delta.Ro))
  for (inv in unique(subd$inversions)){
    rows <- which(subd$inversions == inv)
    mean.eff <- mean(subd$delta.Ro[rows], na.rm=T)
    z[i, inv+1] <- mean.eff
  }
}

z2 <- c(t(z))
df3 <- data.frame(x, y, z2)


#### freq/fixed
z <- matrix(NA, nrow=length(cij), ncol=11)

for (i in 1:length(cij)){
  d <- ComPerm2(pool1, mode="freq", iter=1, Kmeth="fixed", cij=cij[i], exmeth="deterministic")
  subd <- subset(d$data, !is.na(delta.Ro))
  for (inv in unique(subd$inversions)){
    rows <- which(subd$inversions == inv)
    mean.eff <- mean(subd$delta.Ro[rows], na.rm=T)
    z[i, inv+1] <- mean.eff
  }
}

z2 <- c(t(z))
df4 <- data.frame(x, y, z2)


## Combining all contour data and adding a factor identifying transmission & extinction dynamics
alld <- rbind(df1, df2, df3, df4)
alld$group <- rep(c("dens/fixed", "freq/free", "dens/free", "freq/fixed"), each=nrow(df1))

p5 <- ggplot(alld, aes(x=x, y=y, z=z2)) + 
  theme(panel.background = element_rect(color="grey", fill="white"),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  facet_wrap(~group) + scale_x_reverse(limits=c(10, 0)) +
  geom_tile(aes(x=x,y=y,fill=z2), width=2) +
  stat_contour(aes(colour = ..level..), color="black") +
  scale_colour_gradient(low="black", high="black") +
  scale_fill_gradient2(expression(paste(Delta, R[0])), 
                       low="firebrick1", high="blue") 

p5
require(directlabels)
direct.label(p5,list("top.points",cex=.7))

require(lattice)
filled.contour(x=df4$x, y=df4$y, z=df1$z2)
wireframe(z2 ~ x*y|group, data=alld)

ramp <- colorRamp(c("red", "white", "blue"))
contourplot(z2 ~ x*y|group, data=alld, region=T, 
            at=seq(-0.20, .2, .01),
            xlim = rev(extendrange(x)),
            col.regions=rgb(ramp(seq(0, 1, length = 1000)), max = 255))




#########################
### Figure S1 and S2 ####
#########################
# generate data
poold <- NULL
reps <- 1000
seeds <- runif(reps, 0, 100000)
for (i in 1:reps){
  pool <- GPool(globmeth="allom", R0ii="TG", Nglobal=5, a=3, k=k, w=w, seed=seeds[i])
  pool$global.pool$rep <- i
  poold <- rbind(poold, pool$global.pool)
}


### S1 ------------------------------------------------------------------------
test <- ComPerm2(pool, mode="dens", exmeth="deterministic",
                 iter=1, Kmeth="free", cij=.1)
test$pool$inversions <- test$pool$inversions*-1
test$pool$inversions <- test$pool$inversions + 10
test$pool$inversions <- paste(test$pool$inversions, " Inversions")
levs <- NULL
for (i in 1:11){
  levs[i] <- paste(i-1, " Inversions")
}
test$pool$inversions <- factor(test$pool$inversions,
                               levels = levs)
ggplot(test$pool, aes(x=weight, y=Bii, group=permutation)) + theme_classic() +  
  stat_smooth(method="lm", fill="lightgrey", col="black") + 
  geom_point(shape=1, size=2) + 
  xlab("Weight (kg)") + ylab(expression(beta[ii])) +
  theme(axis.title.y = element_text(angle = 0, hjust = 1)) +
  facet_wrap(~inversions, nrow=3) + 
  scale_x_log10() + 
  ggtitle(expression(paste("Figure S1: Relationship between host competence and body mass across ", R[0], " vector inversions"))) +
  theme(
    axis.title.y = element_text(hjust=-.5, size=15),
    axis.title.x = element_text(size=15),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  )


### S2 ------------------------------------------------------------------------
ggplot(poold, aes(x=weight, y=Bii, group=rep)) + theme_classic() +
  geom_line(alpha=0.05) + xlab("Weight (kg)") + ylab(expression(beta[ii])) +
  theme(axis.title.y = element_text(angle = 0, hjust = -.2, size=15),
        plot.margin = unit(c(2, 1, 1, 1), "lines"),
        plot.title = element_text(hjust=0.5, vjust=2)) + 
  scale_x_log10() + 
  labs(title = expression(paste("Figure S2: Simulated host competence, intraspecific ", R[0], ", and body mass", sep="")))
sub <- ggplot(poold, aes(x=R0ii, y=Bii, group=rep)) + theme_classic() +
  geom_line(alpha=0.1) + ylab(expression(beta[ii])) + xlab(expression(paste("Intraspecific ", R[0])))+
  theme(axis.title.y = element_text(angle = 0, hjust = -.2, size=12))
print(sub, vp=viewport(.7, .67, .5, .4))








#################
### Figure S3 ###
#################
#--------------------------#
# Deterministic vs. random #
#--------------------------#
ylims <- c(min(rand.d$delta.Ro), max(rand.d$delta.Ro))
plot.v3 <- function(scenario, ylab=NULL){
  scenario <<- scenario
  require(vioplot)
  red <- "firebrick2"
  blue <-"steelblue"
  botleft <- c(-10, -5)
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
             col="lightgrey", ylim=ylims,
             at=determ.x[which(unique(determ.d$inversions) == i)], 
             side="right", add=T)
    vioplot2(subset(determ.d, scen == scenario & inversions == i)$delta.Ro,
             col="white", 
             at=determ.x[which(unique(determ.d$inversions) == i)], 
             side="left", add=T)
  }
}

cexadj <- .8
in.marg <- .05
nf <- layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=T))
layout.show(nf)
marg.size <- 6
par(oma=c(marg.size-2, marg.size, 5, 1))
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("freq/free")
axis(side=2, cex.axis=cexadj)
par(mai=c(in.marg, in.marg, in.marg, in.marg))
plot.v3("dens/free")
legend(x=6, y=5.4, legend=c("Deterministic extirpations", "Random extirpations"),
       fill=c("white", "lightgrey"), cex=cexadj)
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
mtext("Figure S4: Effects of extirpation on disease risk with deterministic and random extirpations",
      side=3, line=2, outer=T, at=.5, cex=cexadj*1.2)
# Save as 11 X 6 pdf







#################
### Figure S4 ###
#################














#################
### Figure S5 ###
#################
library(distr)
k <- c(.1, .5, 1)
w <- c(1, 1.5, 2)

params <- data.frame(k, w)
params$mean <- rep(NA, nrow(params))
params$median <- rep(NA, nrow(params))
params$sd <- rep(NA, nrow(params))

# Run and store data for all 4 scenarios, recording delta.Ro and %delta.Ro

niter = 1000
d1 <- array(dim=c(1, 4)) # 4 scenarios X 100 communities X number of rows for each subdataframe
colnames(d1) <- c("permutation", "inversions", "delta.Ro", "percent.dR0")


for (i in 1:niter){
  pool1 <- GPool(globmeth="allom", R0ii="TG", Nglobal=5, a=3, k=k[1], w=w[1])
  #### freq/free
  test1 <- ComPerm2(pool1, mode="freq", exmeth="deterministic",
                    iter=1, Kmeth="free", cij=.1)
  sub1 <- subset(test1$data, !is.na(delta.Ro))
  # calculate percent change in R0
  sub1$percent.dR0 <- with(sub1, Ro / (Ro - delta.Ro))
  d1 <- rbind(d1, sub1[, c(9, 10, 7, 11)])
  #### dens/free
  test2 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="free", cij=.1)
  sub2 <- subset(test2$data, !is.na(delta.Ro))
  sub2$percent.dR0 <- with(sub2, Ro / (Ro - delta.Ro))
  d1 <- rbind(d1, sub2[, c(9, 10, 7, 11)])
  #### freq/fixed
  test3 <- ComPerm2(pool1, mode="freq", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub3 <- subset(test3$data, !is.na(delta.Ro))
  sub3$percent.dR0 <- with(sub3, Ro / (Ro - delta.Ro))
  d1 <- rbind(d1, sub3[, c(9, 10, 7, 11)])
  #### dens/fixed
  test4 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub4 <- subset(test4$data, !is.na(delta.Ro))
  sub4$percent.dR0 <- with(sub4, Ro / (Ro - delta.Ro))
  d1 <- rbind(d1, sub4[, c(9, 10, 7, 11)])
}

d1 <- data.frame(d1)
d1 <- d1[-1, ]
scen <- rep(c("freq/free", "dens/free", "freq/fixed", "dens/fixed"), each=nrow(sub1))
scen <- rep(scen, niter)
d1$scen <- scen
d1$case <- "k=0.1, w=1"



# second conditions
d2 <- array(dim=c(1, 4)) # 4 scenarios X 100 communities X number of rows for each subdataframe
colnames(d2) <- c("permutation", "inversions", "delta.Ro", "percent.dR0")

for (i in 1:niter){
  pool1 <- GPool(globmeth="allom", R0ii="TG", Nglobal=5, a=3, k=k[2], w=w[2])
  #### freq/free
  test1 <- ComPerm2(pool1, mode="freq", exmeth="deterministic",
                    iter=1, Kmeth="free", cij=.1)
  sub1 <- subset(test1$data, !is.na(delta.Ro))
  # calculate percent change in R0
  sub1$percent.dR0 <- with(sub1, Ro / (Ro - delta.Ro))
  d2 <- rbind(d2, sub1[, c(9, 10, 7, 11)])
  #### dens/free
  test2 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="free", cij=.1)
  sub2 <- subset(test2$data, !is.na(delta.Ro))
  sub2$percent.dR0 <- with(sub2, Ro / (Ro - delta.Ro))
  d2 <- rbind(d2, sub2[, c(9, 10, 7, 11)])
  #### freq/fixed
  test3 <- ComPerm2(pool1, mode="freq", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub3 <- subset(test3$data, !is.na(delta.Ro))
  sub3$percent.dR0 <- with(sub3, Ro / (Ro - delta.Ro))
  d2 <- rbind(d2, sub3[, c(9, 10, 7, 11)])
  #### dens/fixed
  test4 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub4 <- subset(test4$data, !is.na(delta.Ro))
  sub4$percent.dR0 <- with(sub4, Ro / (Ro - delta.Ro))
  d2 <- rbind(d2, sub4[, c(9, 10, 7, 11)])
}

d2 <- data.frame(d2)
d2 <- d2[-1, ]
scen <- rep(c("freq/free", "dens/free", "freq/fixed", "dens/fixed"), each=nrow(sub1))
scen <- rep(scen, niter)
d2$scen <- scen
d2$case <- "k=0.5, w=1.5"


# third conditions
d3 <- array(dim=c(1, 4)) # 4 scenarios X 100 communities X number of rows for each subdataframe
colnames(d3) <- c("permutation", "inversions", "delta.Ro", "percent.dR0")


for (i in 1:niter){
  pool1 <- GPool(globmeth="allom", R0ii="TG", Nglobal=5, a=3, k=k[3], w=w[3])
  #### freq/free
  test1 <- ComPerm2(pool1, mode="freq", exmeth="deterministic",
                    iter=1, Kmeth="free", cij=.1)
  sub1 <- subset(test1$data, !is.na(delta.Ro))
  # calculate percent change in R0
  sub1$percent.dR0 <- with(sub1, Ro / (Ro - delta.Ro))
  d3 <- rbind(d3, sub1[, c(9, 10, 7, 11)])
  #### dens/free
  test2 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="free", cij=.1)
  sub2 <- subset(test2$data, !is.na(delta.Ro))
  sub2$percent.dR0 <- with(sub2, Ro / (Ro - delta.Ro))
  d3 <- rbind(d3, sub2[, c(9, 10, 7, 11)])
  #### freq/fixed
  test3 <- ComPerm2(pool1, mode="freq", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub3 <- subset(test3$data, !is.na(delta.Ro))
  sub3$percent.dR0 <- with(sub3, Ro / (Ro - delta.Ro))
  d3 <- rbind(d3, sub3[, c(9, 10, 7, 11)])
  #### dens/fixed
  test4 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub4 <- subset(test4$data, !is.na(delta.Ro))
  sub4$percent.dR0 <- with(sub4, Ro / (Ro - delta.Ro))
  d3 <- rbind(d3, sub4[, c(9, 10, 7, 11)])
}

d3 <- data.frame(d3)
d3 <- d3[-1, ]
scen <- rep(c("freq/free", "dens/free", "freq/fixed", "dens/fixed"), each=nrow(sub1))
scen <- rep(scen, niter)
d3$scen <- scen
d3$case <- "k=1, w=2"

# plot showing the density distributions in each case
R0dist <- NULL
for (i in 1:nrow(params)){
  TG <- Truncate(Gammad(scale=k[i], shape=w[i]), lower=0, upper=10)
  Nglobal=1000000
  R0ii <- r(TG)(Nglobal)
  R0dist[[i]] <- R0ii
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

require(ggplot2)
dist1 <- ggplot(data.frame(R0i = R0dist[[1]]), aes(x=R0i)) + theme_classic() +
  geom_density(fill="white") + ylab("Density") + xlab(expression(R[0])) + 
  xlim(0, 10)+ annotate(geom = "text", 
                        label="Gamma distribution: shape=0.1, scale=1", x=6, y=13, size=4)
dist2 <- ggplot(data.frame(R0i=R0dist[[2]]), aes(x=R0i)) + theme_classic() +
  geom_density(fill="lightgray") + ylab("Density") + xlab(expression(R[0]))+ 
  xlim(0, 10) + annotate(geom = "text", 
                         label="Gamma distribution: shape=0.5, scale=1.5", x=6, y=.94, size=4)
dist3 <- ggplot(data.frame(R0i=R0dist[[3]]), aes(x=R0i)) + theme_classic() +
  geom_density(fill="#707070") + ylab("Density") + xlab(expression(R[0]))+ 
  xlim(0, 10) + annotate(geom = "text", 
                         label="Gamma distribution: shape=1, scale=2", x=6, y=.3575, size=4)
require(gridExtra)

ggplot(data.frame(R0i = R0dist[[1]])) + theme_classic() +
  geom_density( aes(x=R0i, y=..density..), fill="white") + ylab("Density") + xlab(expression(R[0])) + 
  geom_density(fill="grey", data=data.frame(R0i=R0dist[[2]]), aes(x=R0i, y=..density..)) +
  geom_density(fill="darkgrey", data=data.frame(R0i=R0dist[[3]]), aes(x=R0i, y=..density..)) +
  xlim(0, 10)+ annotate(geom = "text", 
                        label="Gamma distribution: shape=0.1, scale=1", x=6, y=13, size=4)

# Combine all data 
d4 <- rbind(d1, d2, d3)
d4$inv.fact <- NA
d4$inv.fact[which(d4$inversions==10)] <- "A"
d4$inv.fact[which(d4$inversions==9)] <- "B"
d4$inv.fact[which(d4$inversions==8)] <- "C"
d4$inv.fact[which(d4$inversions==7)] <- "D"
d4$inv.fact[which(d4$inversions==6)] <- "E"
d4$inv.fact[which(d4$inversions==5)] <- "F"
d4$inv.fact[which(d4$inversions==4)] <- "G"
d4$inv.fact[which(d4$inversions==3)] <- "H"
d4$inv.fact[which(d4$inversions==2)] <- "I"
d4$inv.fact[which(d4$inversions==1)] <- "J"
d4$inv.fact[which(d4$inversions==0)] <- "K"



d4$percent.dR0 <- d4$percent.dR0*100
d4$scenario <- NA
d4$scenario[which(d4$scen == "dens/fixed")] <- "DD transmission, compensatory extirpations"
d4$scenario[which(d4$scen == "dens/free")] <- "DD transmission, subtractive extirpations"
d4$scenario[which(d4$scen == "freq/fixed")] <- "FD transmission, compensatory extirpations"
d4$scenario[which(d4$scen == "freq/free")] <- "FD transmission, subtractive extirpations"
box1 <- ggplot(d4, aes(x=inv.fact, y=percent.dR0 - 100, fill=case)) + 
  theme_classic() + 
  geom_boxplot(outlier.size=0) + 
  facet_grid(~scenario) + 
  theme(legend.position="none") + 
  scale_fill_manual(values=c("white", "lightgray", "#707070")) +
  geom_hline(y=0, linetype="dashed") + 
  ylim(-13, 58) +
  ylab(expression(paste("%",Delta, R[0]))) +
  scale_x_discrete(labels = c(" ", "Negative", " ", " ", " ", " ", " ", " ", " ", "Positive", " ")) +
  xlab("Relationship between host competence and extirpation risk")
box2 <- ggplot(d4, aes(x=inv.fact, y=delta.Ro, fill=case)) + theme_classic() + 
  geom_boxplot(outlier.size=0) + facet_grid(~scenario) + 
  theme(legend.position="none") +
  scale_fill_manual(values=c("white", "lightgray", "#707070")) +
  geom_hline(y=0, linetype="dashed") + 
  ylim(-.3, 1) +
  ylab(expression(paste(Delta, R[0]))) +
  scale_x_discrete(labels = c(" ", "Negative", " ", " ", " ", " ", " ", " ", " ", "Positive", " ")) +
  xlab("Relationship between host competence and extirpation risk")
#pdf("test.pdf", width=14, height=10) # or jpeg, etc.
grid.arrange(arrangeGrob(dist1, dist2, dist3, ncol=3), 
             arrangeGrob(box2, box1, ncol=1), 
             ncol=1, heights=c(.3,1))
#dev.off() 
#
