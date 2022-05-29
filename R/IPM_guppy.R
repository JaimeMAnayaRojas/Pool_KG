setwd("~/Dropbox/Projects_JM/FSU/Pool_manipulation/Pool_R/")
pos <- readRDS("Posteriors_KG.RDS")

library(rethinking)
## Get guppy parameters

m.par.GR_G <- data.frame(
  ## survival
  surv<-   cbind(pos$Intercept_survG, pos$b_z_survG),
  ## growth 
  grow    <-  cbind(pos$Intercept_growG, pos$b_z_growG), 
  
  grow.sd   <-  sqrt(pos$sigma_growG), #summary(grow.mod)$sigma,
  ## reproduce or not
  #  repr      <-  cbind(pos$Intercept_repG, pos$b_z_repG), # coef(repr.mod),
  ## recruit or not
  recr      <-  cbind(pos$Intercept_recrG, pos$b_z_recrG) #coef(recr.mod),
  
)

head(m.par.GR_G)
m.par.GR_G <- as.data.frame(m.par.GR_G)
names(m.par.GR_G) <- c("surv.int", "surv.z","grow.int", "grow.z", "grow.sd", 
                       "recr.int", "recr.z")



m.par.NR_G <- list(
  ## survival
  surv<-   cbind(pos$Intercept_survG + pos$b_NK_survG, pos$b_z_survG + pos$b_zNK_survG),
  ## growth 
  grow    <-   cbind(pos$Intercept_growG  + pos$b_NK_growG, pos$b_z_growG + pos$b_zNK_growG), #coef(grow.mod),
  
  grow.sd   <-  sqrt(pos$sigma_growG), #summary(grow.mod)$sigma,
  ## reproduce or not
  #  repr      <-  cbind(pos$Intercept_repG + pos$b_NK_repG, pos$b_z_repG + pos$b_zNK_repG), # coef(repr.mod),
  ## recruit or not
  recr      <-  cbind(pos$Intercept_recrG + pos$b_NK_recrG, pos$b_z_recrG + pos$b_zNK_recrG) #coef(recr.mod),
  
)

m.par.NR_G <- as.data.frame(m.par.NR_G)
names(m.par.NR_G) <- c("surv.int", "surv.z","grow.int", "grow.z", "grow.sd", 
                       "recr.int", "recr.z")







#------------------------ Guppy IPM
#------------------------


nBigMatrix <- 100
# Make the meshpoints
min.size <- (2) #min(IPMdata$z1, na.rm = T) - sd(IPMdata$z1, na.rm = T)*0.5
max.size <- (40) #max(IPMdata$z1, na.rm = T) + sd(IPMdata$z1, na.rm = T)*0.5

U=max.size
L=min.size
m = nBigMatrix
h <- (U - L)/m
meshpts <- z1 <- z <- L + ((1:m) - 1/2) * h
size.cen <- (18)

### Functions
## Growth function, given you are size z now returns the pdf of size z1 next time

# m.par = round(m.par.GR_G[1,],3)

g_z1z <- function(z1, z, m.par){
  
  p.den.grow <- array(NA,c(nBigMatrix,nBigMatrix))
  mean <- as.numeric(m.par[,"grow.int"]) + as.numeric(m.par[,"grow.z"]) * (z -size.cen )           # mean size next year
  sd <- as.numeric(m.par[,"grow.sd"] )                               # sd about mean
  i=1
  for (i in 1:nBigMatrix){
    p.den.grow[,i] <- h*dnorm(z1, mean = mean[i], sd = sd) 
  }
  # pdf that you are size z1 given you were size z
  return((p.den.grow))
}



round(colSums(g_z1z(z1 = meshpts, z = meshpts, m.par = m.par.NR_G[1,])),3)


## Survival function, logistic regression
s_z <- function(z, m.par){
  
  #  linear.p <- logit(0.85) + 0.03 * (z - size.cen)
  linear.p <- as.numeric(m.par[,"surv.int"]) + as.numeric(m.par[,"surv.z"]) * (z - size.cen)       # linear predictor
  p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
  p <- diag(p)
  return(p)
}



plot(colSums(s_z(z = meshpts, m.par = m.par.GR_G[1,]))~ meshpts, type="l", ylim=c(0,1))
lines(y=colSums(s_z(z = meshpts, m.par = m.par.NR_G[1,])), x=  meshpts, type="l", ylim=c(0,1), col="red")


## Reproduction function, logistic regression

## Recruitment function (N.B - from birth in spring to first summer), logistic regression
pr_z <- function(z, m.par) {
  linear.p <- as.numeric(m.par[,"recr.int"]) + as.numeric(m.par[,'recr.z']) * (z - size.cen)                            # linear predictor
  p <- exp(linear.p)                                 # logistic transformation to probability
  p <- diag(p) * (1/2) 
  return(p)
}


plot(colSums(pr_z(z = meshpts, m.par = m.par.GR_G[1,]))~ meshpts, type="l", ylim=c(0,35))
lines(y=colSums(pr_z(z = meshpts, m.par = m.par.NR_G[1,])), x= meshpts, type="l", col="red")

## Recruit size function
c_z1z <- function(z1, z, m.par){
  p.den.rcsz <- array(NA,c(nBigMatrix,nBigMatrix))
  mean <- (7) + 0*z
  #mean <- as.numeric(m.par[,"rcsz.int"]) + as.numeric(m.par[,"rcsz.z"]) * (z - size.cen)           # mean size next year
  
  #sd <- m.par[,"rcsz.sd"]  
  sd <- 0.4# sd about mean
  
  for (i in 1:nBigMatrix){
    p.den.rcsz[,i] <- h*dnorm(z, mean = mean[i], sd = sd) 
  }
  # pdf that offspring are size z1 given you were size z
  return(p.den.rcsz)
}


colSums(c_z1z(z1,z, m.par.GR_G))
#


################ Invididual functions


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
P_z1z <- function (m.par) {
  return( g_z1z(z1, z, m.par) %*% s_z(z, m.par) )
}

## Define the reproduction kernel
F_z1z <- function (m.par) {
  return( c_z1z(z1, z, m.par) %*% pr_z(z,m.par) %*% s_z(z, m.par) )
}

## Build the discretized kernel
mk_K <- function(m, m.par) {
  F <- as.matrix(F_z1z(m.par=m.par))
  P <- as.matrix(P_z1z(m.par=m.par))
  K <- as.matrix(P + F)
  return(list(K = K, meshpts = meshpts, P = P, F = F))
}


###########  IPM For Guppies
#################

G_lamda.est <- as.data.frame(matrix(nrow = dim(m.par.GR_G)[1], ncol = 3))
head(G_lamda.est)

names(G_lamda.est) <- c("lamda.GR","lamda.NR","lamda.dif")
G_lamda.est$Sum.lam.effect<- NA
G_lamda.est$lam.growth<-  NA 
G_lamda.est$lam.fec<- NA 
G_lamda.est$lam.rcz<- NA  
G_lamda.est$lam.surv<- NA 
G_lamda.est$sum.comp.lamda <- NA
G_lamda.est$surv.p <- NA #[k] <- sum(lam.surv) / (sum(lam.gro) + sum(lam.rep) + sum(lam.fec) + sum(lam.rcz) + sum(lam.surv))
G_lamda.est$grow.p <- NA #[k] <- sum(lam.gro)  / (sum(lam.gro) + sum(lam.rep) + sum(lam.fec) + sum(lam.rcz) + sum(lam.surv))
G_lamda.est$fec.p <- NA #[k] <- sum(lam.fec)  / (sum(lam.gro) + sum(lam.rep) + sum(lam.fec) + sum(lam.rcz) + sum(lam.surv)) 
G_lamda.est$rcz.p <- NA #[k] <- sum(lam.rcz)  / (sum(lam.gro) + sum(lam.rep) + sum(lam.fec) + sum(lam.rcz) + sum(lam.surv))




Ggrow.mat <- array(NA,c(dim(G_lamda.est)[1],nBigMatrix))
Gfec.mat <- array(NA,c(dim(G_lamda.est)[1],nBigMatrix))
Grcz.mat <- array(NA,c(dim(G_lamda.est)[1],nBigMatrix))
Gsurv.mat <- array(NA,c(dim(G_lamda.est)[1],nBigMatrix))


k=1

head(m.par.GR_G)
# IPM RUN

start_time = Sys.time()
for (k in 1:dim(m.par.GR_G)[1]){
  
  ## make our projection kernels
  IPM.GR  <- (mk_K(nBigMatrix, m.par = m.par.GR_G[k,]))
  
  ## calculate the population growth rate
  lam.GR  <- Re(eigen( IPM.GR$K)$values[1])
  lam.GR
  G_lamda.est$lamda.GR[k] <- lam.GR
  
  ## make our projection kernels
  IPM.NR  <- mk_K(nBigMatrix, m.par.NR_G[k,])
  ## calculate the population growth rate
  lam.NR  <- Re(eigen( IPM.NR$K)$values[1])
  lam.NR
  G_lamda.est$lamda.NR[k] <- lam.NR
  G_lamda.est$lamda.dif[k] <- lam.NR - lam.GR
  
  # Average K matrix
  K.avg <- ( IPM.GR$K + IPM.NR$K) / 2
  
  ## normalised stable size distribution
  W.GR <- Re(eigen(IPM.GR$K)$vectors[,1]) # dominant right eigenvector
  W.GR <- W.GR / sum(W.GR)
  
  W.NR <- Re(eigen(IPM.NR$K)$vectors[,1]) # dominant right eigenvector
  W.NR <- W.NR / sum(W.NR)
  
  # plot(W.GR ~ meshpts, type="l")
  # lines(y=W.NR, x= meshpts, type="l", col="red")
  # 
  W.avg <- Re(eigen(K.avg)$vectors[,1]) # dominant right eigenvector
  W.avg <- W.avg / sum(W.avg)
  
  
  # reproductive value
  V.avg <- Re(eigen(t(K.avg))$vectors[,1]) # dominant left eigenvector
  V.avg <- V.avg/ c((t(V.avg)%*% W.avg))
  
  # sensitivity matrix
  
  sens.avg <- outer(V.avg,W.avg) 
  K.diff <-  IPM.NR$K - IPM.GR$K
  
  lam.effect <- array(NA,c(nBigMatrix,nBigMatrix))
  
  
  # 
  for (i in 1:nBigMatrix){
    for (j in 1:nBigMatrix){
      lam.effect[i,j] <- K.diff[i,j] * sens.avg[i,j]
    }
  }
  
  
  (G_lamda.est$Sum.lam.effect[k] <-  sum(lam.effect)) # effect matrix
  G_lamda.est$lamda.dif[k]
  
  
  # LIFE RESPONSE TABLE EXPERIMENT
  # Decompose the effects of each demographic parameter on lamda
  #---------------- Perturbations
  a <- 10E-8
  
  #
  one.matrix <- array(1, c(nBigMatrix,nBigMatrix))
  # functions differences
  growth.diff <- (g_z1z(z1, z, m.par = m.par.NR_G[k,]) - g_z1z(z1, z, m.par = m.par.GR_G[k,]))
  fec.diff <- one.matrix%*%(pr_z(z,m.par = m.par.NR_G[k,] ) - pr_z(z,m.par = m.par.GR_G[k,] ) )
  rcz.diff <- c_z1z(z1, z, m.par.NR_G[k,]) - c_z1z(z1, z, m.par.GR_G[k,])
  surv.diff <- one.matrix%*%(s_z(z,m.par = m.par.NR_G[k,]) - s_z(z,m.par = m.par.GR_G[k,] ))
  
  # functions averages
  growth.avg <- (g_z1z(z1, z, m.par = m.par.NR_G[k,]) + g_z1z(z1, z, m.par = m.par.GR_G[k,]))/2
  fec.avg <- one.matrix%*% ((pr_z(z,m.par = m.par.NR_G[k,] ) + pr_z(z,m.par = m.par.GR_G[k,] ) )/2)
  rcz.avg <- (c_z1z(z1, z, m.par.NR_G[k,]) + c_z1z(z1, z, m.par.GR_G[k,]))/2
  surv.avg <- one.matrix %*%(s_z(z,m.par = m.par.NR_G[k,]) + s_z(z,m.par = m.par.GR_G[k,] ) )/2
  
  
  # derivatives
  
  growth.deriv <- surv.avg
  fec.deriv <- rcz.avg * surv.avg
  rcz.deriv <-  fec.avg * surv.avg
  surv.deriv <- growth.avg +  fec.avg * rcz.avg
  
  lam.gro <- array(NA,c(nBigMatrix,nBigMatrix))
  lam.fec <- array(NA,c(nBigMatrix,nBigMatrix))
  lam.rcz <- array(NA,c(nBigMatrix,nBigMatrix))
  lam.surv <- array(NA,c(nBigMatrix,nBigMatrix))
  
  # Ron's original code
  for (i in 1:nBigMatrix){
    for (j in 1:nBigMatrix){
      
      lam.gro[i,j] <- growth.diff[i,j] * sens.avg[i,j] * growth.deriv[i,j] # growth influence
      lam.fec[i,j] <- fec.diff[i,j] * sens.avg[i,j] * fec.deriv[i,j] # fecundity influence
      lam.rcz[i,j] <- rcz.diff[i,j] * sens.avg[i,j] * rcz.deriv[i,j] # offspring size influence
      lam.surv[i,j] <- surv.diff[i,j] * sens.avg[i,j] * surv.deriv[i,j] # surv influence
    }
  }
  
  lam.gro[1:10,1:10]
  
  # Ron's original code
  
  #here is where i have to do the change for the size effects, following alden's advise
  #str(lam.gro)
  G_lamda.est$lam.growth[k]<-  sum(lam.gro)
  G_lamda.est$lam.fec[k]<- sum(lam.fec)
  G_lamda.est$lam.rcz[k]<- sum(lam.rcz)
  G_lamda.est$lam.surv[k]<- sum(lam.surv)
  
  (G_lamda.est$sum.comp.lamda[k] <-  (sum(lam.gro) + sum(lam.fec)
                                      + sum(lam.rcz) + sum(lam.surv))
  )
  
  G_lamda.est$Sum.lam.effect[k]
  
  G_lamda.est$surv.p[k] <- sum(lam.surv) / G_lamda.est$sum.comp.lamda[k] 
  G_lamda.est$grow.p[k] <- sum(lam.gro)  / G_lamda.est$sum.comp.lamda[k] 
  G_lamda.est$fec.p[k] <- sum(lam.fec)  / G_lamda.est$sum.comp.lamda[k] 
  G_lamda.est$rcz.p[k] <- sum(lam.rcz)  / G_lamda.est$sum.comp.lamda[k] 
  
  # Matrices for size effect
  Gsurv.mat[k,] <- (apply(lam.surv, 2, sum)) 
  Ggrow.mat[k,] <- (apply(lam.gro, 2, sum)) 
  Gfec.mat[k,] <- (apply(lam.fec, 2, sum)) 
  Grcz.mat[k,] <- (apply(lam.rcz, 2, sum)) 
  
}

end_time = Sys.time()

end_time - start_time


write.csv(G_lamda.est, "G_lamda.est.csv")
write.csv(Gsurv.mat, "Gsurv.mat.csv")
write.csv(Ggrow.mat, "Ggrow.mat.csv")
write.csv(Gfec.mat, "Gfec.matc.csv")
write.csv(Grcz.mat, "Grcz.mat.csv")