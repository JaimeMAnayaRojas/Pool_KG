setwd("~/Dropbox/Projects_JM/FSU/Pool_manipulation/KG_git/")
pos <- read.csv("Posteriors.csv")

pos[1:5, ]
getwd()
library(rethinking)
## Get parameters for killifish

m.par.GR_K <- data.frame(
  ## survival
  surv<-   cbind(pos$Intercept_survK, pos$b_z_survK),
  ## growth 
  grow    <-  cbind(pos$Intercept_growK, pos$b_z_growK), 
  
  grow.sd   <-  sqrt(pos$sigma_growK), #summary(grow.mod)$sigma,
  ## reproduce or not
  #  repr      <-  cbind(pos$Intercept_repK, pos$b_z_repK), # coef(repr.mod),
  ## recruit or not
  recr      <-  cbind(pos$Intercept_recrK, pos$b_z_recrK) #coef(recr.mod),
  
)

head(m.par.GR_K)
m.par.GR_K <- as.data.frame(m.par.GR_K)
names(m.par.GR_K) <- c("surv.int", "surv.z","grow.int", "grow.z", "grow.sd", 
                       "recr.int", "recr.z")

#--



m.par.NG_K <- list(
  ## survival
  surv<-   cbind(pos$Intercept_survK + pos$b_NG_survK, pos$b_z_survK + pos$b_zNG_survK),
  ## growth 
  grow    <-   cbind(pos$Intercept_growK  + pos$b_NG_growK, pos$b_z_growK + pos$b_zNG_growK), #coef(grow.mod),
  
  grow.sd   <-  sqrt(pos$sigma_growK), #summary(grow.mod)$sigma,
  ## reproduce or not
  #  repr      <-  cbind(pos$Intercept_repK + pos$b_NG_repK, pos$b_z_repK + pos$b_zNG_repK), # coef(repr.mod),
  ## recruit or not
  recr      <-  cbind(pos$Intercept_recrK, pos$b_z_recrK ) #coef(recr.mod),
  
)
m.par.NG_K <- as.data.frame(m.par.NG_K)
names(m.par.NG_K) <- c("surv.int", "surv.z","grow.int", "grow.z", "grow.sd", 
                       "recr.int", "recr.z")

m.par.NG_K[1:5,]
# IPM For Killifish
nBigMatrix <- 100

# Make the meshpoints
min.size <- (2) 
max.size <- (110) 

U=max.size
L=min.size
m = nBigMatrix
h <- (U - L)/m
meshpts <- z1 <- z <- L + ((1:m) - 1/2) * h
size.cen <- (18)

(1:m)
### Functions
## Growth function

g_z1z <- function(z1, z, m.par){
  
  p.den.grow <- array(NA,c(nBigMatrix ,nBigMatrix))
  mean <- ((as.numeric(m.par[,"grow.int"]) + as.numeric(m.par[,"grow.z"]) * (z -size.cen )) - z)/2 + z  # mean size for two weeks
  sd <- as.numeric(m.par[,"grow.sd"] )                    # sd about mean
  i=1
  for (i in 1:nBigMatrix){
    p.den.grow[,i] <- h*dnorm(z1, mean = mean[i], sd = sd) 
  }
  
  matex <- array(0,c(nBigMatrix+1,nBigMatrix+1))
  matex[-1,-1] <- p.den.grow
  
  # pdf that you are size z1 given you were size z
  return((matex))
}

z

round(colSums(g_z1z(z1 = meshpts, z = meshpts, m.par = m.par.NG_K[1,])),6)
m.par = m.par.GR_K[1,]

round(m.par$grow.int, 10)
za = 18
((as.numeric(m.par[,"grow.int"]) + as.numeric(m.par[,"grow.z"]) * (za -size.cen )) - za) + za


## Survival function, logistic regression
s_z <- function(z, m.par){
  
  #  linear.p <- logit(0.85) + 0.03 * (z - size.cen)
  linear.p <- as.numeric(m.par[,"surv.int"]) + as.numeric(m.par[,"surv.z"]) * (z - size.cen)       # linear predictor
  p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
  p <- diag(sqrt(p))
  
  matex <- array(0,c(nBigMatrix+1,nBigMatrix+1))
  matex[-1,-1] <- p
  
  return(matex)
}


s_z(z = meshpts, m.par = m.par.GR_K[1,])[1:5,1:5]



plot(colSums(s_z(z = meshpts, m.par = m.par.GR_K[1,])[-1,-1])~ meshpts, type="l", ylim=c(0,1))
lines(colSums(s_z(z = meshpts, m.par = m.par.NG_K[1,])[-1,-1])~ meshpts, col="red")




## Recruitment function

pr_z <- function(z, m.par) {
  linear.p <- as.numeric(m.par[,"recr.int"]) + as.numeric(m.par[,'recr.z']) * (z - size.cen)                            # linear predictor
  p <- exp(linear.p)* (1/2)                                 # logistic transformation to probability
  #  p <- diag(p) * (1/2) 
  matex <- as.data.frame(array(0,c(nBigMatrix+1,nBigMatrix+1)))
  matex[1,2:dim(matex)[2]] <- p
  
  return(as.matrix(matex))
}

m.par =m.par.GR_K[1,]

## Recruit size function
c_z1z <- function(z1, z, m.par){
  p.den.rcsz <- array(0,c(nBigMatrix+1,nBigMatrix+1))
  z0 <- c(0,z)
  
  mean <- (4.4) + 0*(z-18)
  #mean <- as.numeric(m.par[,"rcsz.int"]) + as.numeric(m.par[,"rcsz.z"]) * (z - size.cen)           # mean size next year
  
  #sd <- m.par[,"rcsz.sd"]  
  sd <- 0.6# sd about mean
  
  
  #  for (i in 1:(nBigMatrix)){
  
  p.den.rcsz[c(2:dim(p.den.rcsz)[1]),1] <- h*dnorm(z0[-1], mean = mean, sd = sd) 
  
  #   }
  # pdf that offspring are size z1 given you were size z
  return(p.den.rcsz)
}


colSums(c_z1z(z1,z, m.par.GR_K[1,]))

################ Invididual functions

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Build the kernel

#### Question for Ron: 
# Here is where things might be wrong. 
# Following Ron's paper, I should project A across two weeks.
# The measurements that we have correspond to 4 weeks, so logic dictates that I should divide the vital functions by 2.
# but I am not sure about it... 
# If my logic is correct, should I divide by 2 the vital functions, or the  kernel?

mk_K <- function(m, m.par, v=0.7) {
  F <- as.matrix(pr_z(z,m.par)%*%s_z(z, m.par)  )
  P <- as.matrix(g_z1z(z,z,m.par)%*%s_z(z, m.par))
  A <- (P + F) + v * c_z1z(z,z, m.par)
  K <-  A%*%A # project the population across the full census
  return(list(K = (K), meshpts = meshpts, P = P, F = F, A = A))
}



K_lamda.est <- as.data.frame(matrix(nrow = dim(m.par.GR_K)[1], ncol = 3))
head(K_lamda.est)

names(K_lamda.est) <- c("lamda.GR","lamda.NG","lamda.dif")
K_lamda.est$Sum.lam.effect<- NA
K_lamda.est$lam.growth<-  NA

K_lamda.est$lam.fec<- NA
K_lamda.est$lam.rcz<- NA
K_lamda.est$lam.surv<- NA
K_lamda.est$sum.comp.lamda <- NA
K_lamda.est$surv.p <- NA #[k] <- sum(lam.surv) / (sum(lam.gro) + sum(lam.rep) + sum(lam.fec) + sum(lam.rcz) + sum(lam.surv))
K_lamda.est$grow.p <- NA #[k] <- sum(lam.gro)  / (sum(lam.gro) + sum(lam.rep) + sum(lam.fec) + sum(lam.rcz) + sum(lam.surv))
K_lamda.est$fec.p <- NA #[k] <- sum(lam.fec)  / (sum(lam.gro) + sum(lam.rep) + sum(lam.fec) + sum(lam.rcz) + sum(lam.surv))
K_lamda.est$rcz.p <- NA #[k] <- sum(lam.rcz)  / (sum(lam.gro) + sum(lam.rep) + sum(lam.fec) + sum(lam.rcz) + sum(lam.surv))




Kgrow.mat <- array(NA,c(dim(K_lamda.est)[1],(nBigMatrix+1)))
Kfec.mat <- array(NA,c(dim(K_lamda.est)[1],nBigMatrix+1))
Krcz.mat <- array(NA,c(dim(K_lamda.est)[1],nBigMatrix)+1)
Ksurv.mat <- array(NA,c(dim(K_lamda.est)[1],nBigMatrix+1))


precis(m.par.GR_K)
precis(m.par.NG_K)
# IPM RUN
# k=1
# 
# 
# max(m.par.GR_K$grow.int)
# max(m.par.GR_K$grow.z)
# 
# ?sample
# 
# rsam = sample(1:dim(m.par.GR_K)[1], size = 2000)
# m.par.GR_K  = m.par.GR_K[rsam,]
# m.par.NG_K  = m.par.NG_K[rsam,]

k = 1

for (k in 1:dim(m.par.GR_K)[1]){
  
  ## make our projection kernels
  IPM.GR  <- (mk_K(nBigMatrix, m.par = m.par.GR_K[k,]))
  
  ## calculate the population growth rate
  lam.GR  <- Re(eigen( IPM.GR$K)$values[1])
  (K_lamda.est$lamda.GR[k] <- lam.GR)
  
  ## make our projection kernels
  IPM.NG  <- mk_K(nBigMatrix, m.par.NG_K[k,])
  ## calculate the population growth rate
  lam.NG  <- Re(eigen( IPM.NG$K)$values[1])
  (K_lamda.est$lamda.NG[k] <- lam.NG)
  (K_lamda.est$lamda.dif[k] <- lam.NG - lam.GR)
  
  
  
  
  
  # Average K matrix
  K.avg <-( IPM.GR$K + IPM.NG$K) / 2
  
  ## normalised stable size distribution
  W.GR <- Re(eigen(IPM.GR$K)$vectors[,1]) # dominant right eigenvector
  W.GR <- W.GR / sum(W.GR)
  
  W.NG <- Re(eigen(IPM.NG$K)$vectors[,1]) # dominant right eigenvector
  W.NG <- W.NG / sum(W.NG)
  
  W.avg <- Re(eigen(K.avg)$vectors[,1]) # dominant right eigenvector
  W.avg <- W.avg / sum(W.avg)
  
  
  # reproductive value
  V.avg <- Re(eigen(t(K.avg))$vectors[,1]) # dominant left eigenvector
  V.avg <- V.avg/ c((t(V.avg)%*% W.avg))
  
  # sensitivity matrix
  
  
  sens.avg <- outer(V.avg,W.avg)
  K.diff <-  IPM.NG$K - IPM.GR$K
  
  lam.effect <-  K.diff * sens.avg 
  
  K_lamda.est$Sum.lam.effect[k] <-  sum(lam.effect) # effect matrix
  K_lamda.est$lamda.dif[k]
  
  
  # LIFE RESPONSE TABLE EXPERIMENT
  # Decompose the effects of each demographic parameter on lamda
  
  one.matrix <- array(1, c(nBigMatrix+1,nBigMatrix+1))
  # functions differences
  growth.diff <- (g_z1z(z1, z, m.par = m.par.NG_K[k,]) - g_z1z(z1, z, m.par = m.par.GR_K[k,]))
  fec.diff <- (pr_z(z,m.par = m.par.NG_K[k,] ) - pr_z(z,m.par = m.par.GR_K[k,] ) )
  rcz.diff <- c_z1z(z1, z, m.par.NG_K[k,]) - c_z1z(z1, z, m.par.GR_K[k,])
  surv.diff <- (s_z(z,m.par = m.par.NG_K[k,]) - s_z(z,m.par = m.par.GR_K[k,] ))
  
  # functions averages
  growth.avg <- (g_z1z(z1, z, m.par = m.par.NG_K[k,]) + g_z1z(z1, z, m.par = m.par.GR_K[k,]))/2
  fec.avg <-  ((pr_z(z,m.par = m.par.NG_K[k,] ) + pr_z(z,m.par = m.par.GR_K[k,] ) )/2)
  rcz.avg <- (c_z1z(z1, z, m.par.NG_K[k,]) + c_z1z(z1, z, m.par.GR_K[k,]))/2
  surv.avg <- (s_z(z,m.par = m.par.NG_K[k,]) + s_z(z,m.par = m.par.GR_K[k,] ) )/2
  
  
  # derivatives
  I <- diag(rep(1, (nBigMatrix+1)))
  # growth.deriv[1:5,1:5]
  # dim(growth.deriv)
  # 
  # W = G%*%S
  # dW = I%*%dG%*%S + G%*%dS%*%I
  # dvecW = (t(S) %x% I) %*% dvecG + (I %x% G) %*% dvecS
  # 
  # dvecW / dvecG = t(S) %x% I + (I %x% G) %*% dvecS / dvecG
  # dvecW / dvecG = t(S) %x% I 
  # 
  # Q = DFS
  # 
  # K = W + Q = G%*%S + D%*%F%*%S
  # dK / dS:
  #   dK/ dS =      I%*%dG%*%S +     G%*%dS%*%I      + I%*%dD%*%F%*%S + D%*%dF%*%S + D%*%F%*%dS%*%I
  # dvecK /dvecS=       0       + (t(I)%x%G) %*% dvecS +     0         +      0     + (t(I) %x% (D%*%F)) %*% dvecS
  # 
  # dvecK / dvecS = (t(I)%x%G) + (t(I) %x% (D%*%F))
  # 
  # Brute force method:
  #   K.unpert
  # 1) Take the vital rates and increase element 1,1 in G by 0.001. so itwas 0.3 it will now be 0.301
  # 2) Use that G to make a new K and estimate lambda.
  # 3) dlambda /  dG11  ~  (lambda of K.pert - lambda of K.unpert) / 0.001
  # 
  # 
  # derivatives
  I <- diag(rep(1, (nBigMatrix+1)))
  growth.deriv <- t(surv.avg)%x%I
  
  surv.deriv <- (t(I) %x% growth.avg) +  (t(I) %x% (fec.avg %*% rcz.avg))
  fec.deriv <- t(surv.avg) %x% rcz.avg 
  rcz.deriv <-  t(fec.avg %*% surv.avg) %x% I  
  
  dA <- (I %x% K.avg) + ( t(K.avg) %x% I)#
  
  lam.gro <- (growth.diff) * array((array(sens.avg,c(1,(nBigMatrix+1)^2)) %*% dA %*% growth.deriv) ,c(nBigMatrix+1,nBigMatrix+1)) # growth influence
  lam.fec <- fec.diff *array( (array(sens.avg,c(1,(nBigMatrix+1)^2)) %*% dA %*% fec.deriv) ,c(nBigMatrix+1,nBigMatrix+1)) # fecundity influence
  lam.rcz <- rcz.diff * array( (array(sens.avg,c(1,(nBigMatrix+1)^2)) %*% dA %*% rcz.deriv) ,c(nBigMatrix+1,nBigMatrix+1))   # offspring size influence
  lam.surv <- surv.diff * array( (array(sens.avg,c(1,(nBigMatrix+1)^2)) %*% dA %*% surv.deriv) ,c(nBigMatrix+1,nBigMatrix+1))   # surv influence
  
  #here is where i have to do the change for the size effects, following alden's advise
  #str(lam.gro)
  K_lamda.est$lam.growth[k]<-  sum(lam.gro, na.rm = T)
  #K_lamda.est$lam.rep[k]<- sum(lam.rep, na.rm = T)
  K_lamda.est$lam.fec[k]<- sum(lam.fec, na.rm = T)
  K_lamda.est$lam.rcz[k]<- sum(lam.rcz, na.rm = T)
  K_lamda.est$lam.surv[k]<- sum(lam.surv, na.rm = T)
  
  (K_lamda.est$sum.comp.lamda[k] <-  (sum(lam.gro, na.rm = T) + sum(lam.fec, na.rm = T)
                                      + sum(lam.rcz, na.rm = T) + sum(lam.surv, na.rm = T)))
  
  K_lamda.est$lamda.dif[k]
  
  K_lamda.est$surv.p[k] <- sum(lam.surv, na.rm = T) / K_lamda.est$sum.comp.lamda[k]
  K_lamda.est$grow.p[k] <- sum(lam.gro, na.rm = T)  / K_lamda.est$sum.comp.lamda[k]
  K_lamda.est$fec.p[k] <- sum(lam.fec, na.rm = T)  / K_lamda.est$sum.comp.lamda[k]
  K_lamda.est$rcz.p[k] <- sum(lam.rcz, na.rm = T)  / K_lamda.est$sum.comp.lamda[k]
  
  
  Ksurv.mat[k,] <- (apply(lam.surv, c(2), sum))
  Kgrow.mat[k,] <- (apply(lam.gro, c(2), sum))
  Kfec.mat[k,] <- (apply(lam.fec, c(2), sum))
  Krcz.mat[k,] <- (apply(lam.rcz, c(2), sum))
  
  
}

# End of IPM loop

LOS <- function(x=NULL){
  
  out =   (length(which(x > 0)) / length(x))*100
  
  return(round(out, 3))
  
}

write.csv(K_lamda.est, "K_lamda.est.csv" )
write.csv(Ksurv.mat, "Ksurv.mat.csv")
write.csv(Kgrow.mat, "Kgrow.mat.csv")
write.csv(Kfec.mat, "Kfec.matc.csv")
write.csv(Krcz.mat, "Krcz.mat.csv")