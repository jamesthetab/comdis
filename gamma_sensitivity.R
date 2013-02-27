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
  sub1 <- subset(test1$data, !is.na(delta.Ro) & inversions >= median(inversions))
  # calculate percent change in R0
  sub1$percent.dR0 <- with(sub1, Ro / (Ro - delta.Ro))
  d1 <- rbind(d1, sub1[, c(9, 10, 7, 11)])
  #### dens/free
  test2 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="free", cij=.1)
  sub2 <- subset(test2$data, !is.na(delta.Ro) & inversions >= median(inversions))
  sub2$percent.dR0 <- with(sub2, Ro / (Ro - delta.Ro))
  d1 <- rbind(d1, sub2[, c(9, 10, 7, 11)])
  #### freq/fixed
  test3 <- ComPerm2(pool1, mode="freq", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub3 <- subset(test3$data, !is.na(delta.Ro) & inversions >= median(inversions))
  sub3$percent.dR0 <- with(sub3, Ro / (Ro - delta.Ro))
  d1 <- rbind(d1, sub3[, c(9, 10, 7, 11)])
  #### dens/fixed
  test4 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub4 <- subset(test4$data, !is.na(delta.Ro) & inversions >= median(inversions))
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
  sub1 <- subset(test1$data, !is.na(delta.Ro) & inversions >= median(inversions))
  # calculate percent change in R0
  sub1$percent.dR0 <- with(sub1, Ro / (Ro - delta.Ro))
  d2 <- rbind(d2, sub1[, c(9, 10, 7, 11)])
  #### dens/free
  test2 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="free", cij=.1)
  sub2 <- subset(test2$data, !is.na(delta.Ro) & inversions >= median(inversions))
  sub2$percent.dR0 <- with(sub2, Ro / (Ro - delta.Ro))
  d2 <- rbind(d2, sub2[, c(9, 10, 7, 11)])
  #### freq/fixed
  test3 <- ComPerm2(pool1, mode="freq", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub3 <- subset(test3$data, !is.na(delta.Ro) & inversions >= median(inversions))
  sub3$percent.dR0 <- with(sub3, Ro / (Ro - delta.Ro))
  d2 <- rbind(d2, sub3[, c(9, 10, 7, 11)])
  #### dens/fixed
  test4 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub4 <- subset(test4$data, !is.na(delta.Ro) & inversions >= median(inversions))
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
  sub1 <- subset(test1$data, !is.na(delta.Ro) & inversions >= median(inversions))
  # calculate percent change in R0
  sub1$percent.dR0 <- with(sub1, Ro / (Ro - delta.Ro))
  d3 <- rbind(d3, sub1[, c(9, 10, 7, 11)])
  #### dens/free
  test2 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="free", cij=.1)
  sub2 <- subset(test2$data, !is.na(delta.Ro) & inversions >= median(inversions))
  sub2$percent.dR0 <- with(sub2, Ro / (Ro - delta.Ro))
  d3 <- rbind(d3, sub2[, c(9, 10, 7, 11)])
  #### freq/fixed
  test3 <- ComPerm2(pool1, mode="freq", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub3 <- subset(test3$data, !is.na(delta.Ro) & inversions >= median(inversions))
  sub3$percent.dR0 <- with(sub3, Ro / (Ro - delta.Ro))
  d3 <- rbind(d3, sub3[, c(9, 10, 7, 11)])
  #### dens/fixed
  test4 <- ComPerm2(pool1, mode="dens", exmeth="deterministic",
                    iter=1, Kmeth="fixed", cij=.1)
  sub4 <- subset(test4$data, !is.na(delta.Ro) & inversions >= median(inversions))
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
g1 <- grid.arrange(dist1, dist2, dist3, ncol=3)



# Combine all data 
d4 <- rbind(d1, d2, d3)
d4$inv.fact <- NA
d4$inv.fact[which(d4$inversions==10)] <- "A"
d4$inv.fact[which(d4$inversions==9)] <- "B"
d4$inv.fact[which(d4$inversions==8)] <- "C"
d4$inv.fact[which(d4$inversions==7)] <- "D"
d4$inv.fact[which(d4$inversions==6)] <- "E"
d4$inv.fact[which(d4$inversions==5)] <- "F"
d4$percent.dR0 <- d4$percent.dR0*100
d4$scenario <- NA
d4$scenario[which(d4$scen == "dens/fixed")] <- "DD transmission, compensatory extirpations"
d4$scenario[which(d4$scen == "dens/free")] <- "DD transmission, subtractive extirpations"
d4$scenario[which(d4$scen == "freq/fixed")] <- "FD transmission, compensatory extirpations"
d4$scenario[which(d4$scen == "freq/free")] <- "FD transmission, subtractive extirpations"
box1 <- ggplot(d4, aes(x=inv.fact, y=percent.dR0 - 100, fill=case)) + theme_classic() + 
  geom_boxplot(outlier.size=0) + facet_grid(~scenario) + 
  theme(legend.position="none") + 
  scale_fill_manual(values=c("white", "lightgray", "#707070")) +
  geom_hline(y=0, linetype="dashed") + 
  ylim(-50, 70) +
  ylab(expression(paste("%",Delta, R[0]))) +
  scale_x_discrete(labels = c("Negative", " ", " ", " ", " ", "Random")) +
  xlab("Relationship between host competence and extirpation risk")
box2 <- ggplot(d4, aes(x=inv.fact, y=delta.Ro, fill=case)) + theme_classic() + 
  geom_boxplot(outlier.size=0) + facet_grid(~scenario) + 
  theme(legend.position="none") +
  scale_fill_manual(values=c("white", "lightgray", "#707070")) +
  geom_hline(y=0, linetype="dashed") + 
  ylim(-.3, 1) +
  ylab(expression(paste(Delta, R[0]))) +
  scale_x_discrete(labels = c("Negative", " ", " ", " ", " ", "Random")) +
  xlab("Relationship between host competence and extirpation risk")
g2 <- grid.arrange(box1, box2, ncol=1)
plot5<-grid.arrange(arrangeGrob(dist1, dist2, dist3, ncol=1), 
                    arrangeGrob(box1, box2, ncol=1), 
                    ncol=2, widths=c(.5,1))
#pdf("test.pdf", width=14, height=10) # or jpeg, etc.
grid.arrange(arrangeGrob(dist1, dist2, dist3, ncol=3), 
             arrangeGrob(box2, box1, ncol=1), 
             ncol=1, heights=c(.3,1))
#dev.off() 
#


























params <- expand.grid(k, w)
colnames(params) <- c("k", "w")
params$mean <- rep(NA, nrow(params))
params$median <- rep(NA, nrow(params))
params$sd <- rep(NA, nrow(params))

Nglobal=5

for (i in 1:nrow(params)){
  TG <- Truncate(Gammad(scale=params$k[i], shape=params$w[i]), 
                 lower=0, upper=10)
  R0ii <- r(TG)(Nglobal)
  params$mean[i] <- mean(R0ii)
  params$median[i] <- median(R0ii)
  params$sd[i] <- sd(R0ii)
}

plot(params$mean, params$sd)

scen <- c("f/s", "d/s", "f/c", "d/c")

sens.dat <- rbind(params, params, params, params)
sens.dat$scen <- rep(scen, each=nrow(params))
sens.dat$p.dilut <- NA

for (i in 1:nrow(params)){
  # create global pool for each k/w combination
  pool <- GPool(Nglobal=5, globmeth="allom",
                a=1.15, m=1.5, R0ii="TG", min.weight=1, kRecovery=10,
                k = params$k[i], w = params$w[i])
  # then under the 4 scenarios of interest, calculate 
  # 1) the proportion of permutations with only dilution
  #    - one value per scenario
  #      - 4 values per global pool
  # 2) transition points from net dilution to net amplification
  #    - one value per scenario (either NA if there is no transition, 
  #      or the value beyond which there is net amplification)
  #      - 4 values per global pool
  
  # freq/free
  out <- ComPerm2(pool, mode="freq", exmeth="deterministic",
                    iter=1, Kmeth="free", cij=.1)
  count.dilution <- sum(out$data$delta.Ro > 0, na.rm=T)
  count.total <- sum(!is.na(out$data$delta.Ro))
  p.dilut <- count.dilution / count.total
  # then plug that value into the output data frame
  indx2 <- which(sens.dat$k == params$k[i] & sens.dat$w == params$w[i] & sens.dat$scen == "f/s")
  sens.dat$p.dilut[indx2] <- p.dilut
  
  # dens/free
  out <- ComPerm2(pool, mode="dens", exmeth="deterministic",
                  iter=1, Kmeth="free", cij=.1)
  count.dilution <- sum(out$data$delta.Ro > 0, na.rm=T)
  count.total <- sum(!is.na(out$data$delta.Ro))
  p.dilut <- count.dilution / count.total
  # then plug that value into the output data frame
  indx2 <- which(sens.dat$k == params$k[i] & sens.dat$w == params$w[i] & sens.dat$scen == "d/s")
  sens.dat$p.dilut[indx2] <- p.dilut
  
  # freq/fixed
  out <- ComPerm2(pool, mode="freq", exmeth="deterministic",
                  iter=1, Kmeth="fixed", cij=.1)
  count.dilution <- sum(out$data$delta.Ro > 0, na.rm=T)
  count.total <- sum(!is.na(out$data$delta.Ro))
  p.dilut <- count.dilution / count.total
  # then plug that value into the output data frame
  indx2 <- which(sens.dat$k == params$k[i] & sens.dat$w == params$w[i] & sens.dat$scen == "f/c")
  sens.dat$p.dilut[indx2] <- p.dilut
  
  # dens/fixed
  out <- ComPerm2(pool, mode="dens", exmeth="deterministic",
                  iter=1, Kmeth="fixed", cij=.1)
  count.dilution <- sum(out$data$delta.Ro > 0, na.rm=T)
  count.total <- sum(!is.na(out$data$delta.Ro))
  p.dilut <- count.dilution / count.total
  # then plug that value into the output data frame
  indx2 <- which(sens.dat$k == params$k[i] & sens.dat$w == params$w[i] & sens.dat$scen == "d/c")
  sens.dat$p.dilut[indx2] <- p.dilut
}
require(rgl)
ggplot(subset(sens.dat, scen=="d/c"), aes(x=sd, y=p.dilut)) + 
  geom_point() + stat_smooth()
scatter3d(p.dilut ~ mean + sd, fit="smooth",
          data=subset(sens.dat, scen=="f/c"))

# 4 scenarios
# freq/free
# dens/free
# freq/fixed
# dens/fixed
require(lattice)
wireframe(p.dilut ~ mean*sd|scen, data=sens.dat, region=T, allow.multiple=T)
# each combination of shape and scale parameters specifies some distribution of R0ii values
# used to construct a global pool

# so, we want to construct a global pool for each k and w combination

# then under the 4 scenarios of interest, calculate 
# 1) the proportion of permutations with only dilution
#    - one value per scenario
#      - 4 values per global pool
# 2) transition points from net dilution to net amplification
#    - one value per scenario (either NA if there is no transition, 
#      or the value beyond which there is net amplification)
#      - 4 values per global pool
