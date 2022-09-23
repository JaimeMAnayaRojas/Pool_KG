####  
library(data.table)
library(brms)
library(ggplot2)
library(HDInterval)

rm(list=ls(all=TRUE))



###########################################################################################################
# First get the data
getwd()
#setwd("~/Dropbox/Jaime M/Projects_JM/FSU/Pool_manipulation/Pool_KG/")
source("R/Functions.R")

center = 18
# Guppy data cleaning -----------------------------------------------------

Gdata <- read.csv("data/GuppyIPM.csv")
Kdata <- read.csv("data/KillifishIPM.csv")

Gdata$z <- Gdata$SL1_mm - center
Gdata$z2 <- Gdata$SL1_mm^2 - center^2
Gdata$z1 <- Gdata$SL2_mm 

# Make sure the non reproductive have zeros
Gdata$Recr[which(Gdata$Repr==0 & Gdata$surv ==1)] = 0



# Killifish data cleaning -------------------------------------------------

Kdata$z <- Kdata$SL1_mm - center
Kdata$z1 <- Kdata$SL2_mm 
Kdata$z2 <- Kdata$SL1_mm^2 - center^2


# Make sure the non reproductive have zeros
Kdata$Recr[which(Kdata$Repr==0 & Kdata$surv ==1)] = 0



# Normalize the covariates within streams -----------------

Gdata$FishBiom <- Gdata$BiomassG1 + Gdata$BiomassK1
Kdata$FishBiom <- Kdata$BiomassG1 + Kdata$BiomassK1

Gdata$Density <-( Gdata$BiomassG1 + Gdata$BiomassK1) / Gdata$area
Kdata$Density <- (Kdata$BiomassG1 + Kdata$BiomassK1) / Kdata$area


library("plyr")
names(Gdata)
vars = c("site", "Location", "Pool_1", "area", "FishBiom", "canopy", "Density")
df = rbind(Gdata[,vars], Kdata[,vars])
df = ddply(df,c('site','Location','Pool_1'),summarise, area=mean(area), Biomass=mean(FishBiom), canopy=mean(canopy), Density = mean(Density))

means = ddply(df,c('Location'),summarise, mean_area=mean(area), mean_Biomass=mean(Biomass), mean_canopy=mean(canopy), mean_Density= mean(Density))
sd = ddply(df,c('Location'),summarise, sd_area=sd(area), sd_Biomass=sd(Biomass), sd_canopy=sd(canopy), sd_Density= sd(Density))

Gdata = merge(Gdata, means, by.x = "Location")
Gdata = merge(Gdata, sd, by.x = "Location")

Kdata = merge(Kdata, means, by.x = "Location")
Kdata = merge(Kdata, sd, by.x = "Location")




Gdata$FishBiom <- (Gdata$FishBiom - Gdata$mean_Biomass) /  Gdata$sd_Biomass
Kdata$FishBiom <- (Kdata$FishBiom - Kdata$mean_Biomass) /  Kdata$sd_Biomass

Gdata$area <- (Gdata$area - Gdata$mean_area) /  Gdata$sd_area
Kdata$area <- (Kdata$area - Kdata$mean_area) /  Kdata$sd_area


Gdata$canopy <- (Gdata$canopy - Gdata$mean_canopy) /  Gdata$sd_canopy
Kdata$canopy <- (Kdata$canopy - Kdata$mean_canopy) /  Kdata$sd_canopy


Gdata$Density <- (Gdata$Density - Gdata$mean_Density) /  Gdata$sd_Density
Kdata$Density <- (Kdata$Density - Kdata$mean_Density) /  Kdata$sd_Density



# Re-code the random effects (drainage id) --------------------------------


Kdata$stream <- factor(Kdata$Location)
levels(Kdata$stream) <- 1:length(levels(Kdata$stream))
Kdata$stream <- as.numeric(Kdata$stream)

Gdata$stream <- factor(Gdata$Location)
levels(Gdata$stream) <- 1:length(levels(Gdata$stream))
Gdata$stream <- as.numeric(Gdata$stream)

Gdata$z2 = Gdata$SL1_mm^2 - center^2
Kdata$z2 = Kdata$SL1_mm^2 - center^2


Gdata$growth = log(Gdata$SL2_mm/Gdata$SL1_mm)
Kdata$growth = log(Kdata$SL2_mm/Kdata$SL1_mm)



# Model selection for Guppies ---------------------------------------------

# Survival model ----------------------------------------------------------
# Guppies
S1 <- brm(surv ~ z + z2 + NK + z:NK + z2:NK  + Density + canopy + (1|stream), family = bernoulli(), Gdata,
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 10))
S1 <- add_criterion(S1, c("loo", "waic"))
## remove effects

S2 <- update(S1, formula. = ~ . -  z2:NK)
S2 <- add_criterion(S2, c("loo", "waic"))


S3 <- update(S2, formula. = ~ . -  z:NK)
S3 <- add_criterion(S3, c("loo", "waic"))


S4 <- update(S3, formula. = ~ . -  z2)
S4 <- add_criterion(S4, c("loo", "waic"))

(loo_compare(S1,S2,S3,S4, criterion = "loo"))
loo_model_weights(S1,S2,S3,S4, criterion = "loo")

summary(S4)
postS = data.frame(posterior_samples(S4))


Gdata$Treatment = ifelse(Gdata$NK == 1, "NK", "KG" ) 
KG = inv_logit_scaled(post$b_Intercept) *100
NK = inv_logit_scaled(post$b_Intercept + post$b_NK) *100


hdi(NK, credMass = .95)[1]

postS2 = posterior_samples(S2)


P_link= function(post, NK, size){
  z = size - 18
  z2 = (size^2) - 18^2
  
  mu = post$b_Intercept + post$b_NK * NK + (post$b_z + post$`b_z:NK`*NK) * z + post$b_z2*z2
  
  return(mu)
  
}

P_link(postS, NK = 0,size = 5)

M4_link= function(post, NK, size){
  z = size - 18
  z2 = (size^2) - 18^2
  
  mu = post$b_Intercept + post$b_NK * NK + post$b_z * z
  
  return(mu)
  
}

summary(S2)
postS = posterior_samples(S2)


size = seq(from=5, to=32, length= 100)
NGs = 100* inv_logit_scaled(sapply(1:length(size), function(i) P_link(post = postS, NK = 0, size = size[i])))
NGs = cbind(apply(NGs, 2, mean),
            t(apply(NGs, 2, hdi))[,1],
            t(apply(NGs, 2, hdi))[,2])

NKs = 100* inv_logit_scaled(sapply(1:length(size), function(i) P_link(post = postS, NK = 1, size = size[i])))
NKs = cbind(apply(NKs, 2, mean),
            t(apply(NKs, 2, hdi))[,1],
            t(apply(NKs, 2, hdi))[,2])



dfs = as.data.frame(rbind(NGs,NKs ))
names(dfs) <- c("surv", "l95", "u95")
dfs$Treatment = c(rep("KG", dim(NKs)[1]), rep("NK", dim(NKs)[1]))
dfs$z = c(size, size)

names(dfs)

Gdata$Treatment <- ifelse(Gdata$NK==1, "NK", "KG")

splot = Gdata[,c("SL1_mm", "Treatment", "surv")]
names(splot)[1] = "z"
splot$surv = splot$surv * 100
splot

library(latex2exp)
(plotA = ggplot(dfs, aes(x = z, y = surv, colour = Treatment)) + 
  geom_line(size = 1.5) + 
  ylab("Survival (%)") +
  xlab("Initial size (mm)") +
  geom_point(data=splot, aes(x=z, y=surv, colour= Treatment), alpha = 0.25, size = 0.75) +
  ylim(0, 100) + theme_classic() + theme(panel.grid.major = element_blank(), legend.position = "none") +
  geom_ribbon( aes(ymin = l95, ymax = u95, fill= Treatment), 
              alpha = 0.2, size = 0.01, show.legend = F) + 
  scale_color_manual(values=c("black", "orange")) +
  scale_fill_manual(values=c("black", "orange"))  + geom_vline(xintercept=18, lty=2, size=0.2)+
  xlim(10,30)
)


## Guppy growth

G1 <- brm(growth ~ z + z2 + NK + z:NK + z2:NK  + Density + canopy + (1|stream), family = gaussian(), subset(Gdata, surv == 1), 
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 11))
G1 <- add_criterion(G1, c("loo", "waic"))

summary(G1)
## remove effects

G2 <- update(G1, formula. = ~ . -  z2:NK)
G2 <- add_criterion(G2, c("loo", "waic"))


G3 <- update(G2, formula. = ~ . -  z:NK)
G3 <- add_criterion(G3, c("loo", "waic"))


G4 <- update(G3, formula. = ~ . -  z2)
G4 <- add_criterion(G4, c("loo", "waic"))



(loo_compare(G1,G2,G3,G4, criterion = "loo"))
loo_model_weights(G1,G2,G3,G4, criterion = "loo")

summary(G2)

M3_link= function(post, NK, size){
  z = size - 18
  z2 = (size^2 - 18^2)

 mu = post$b_Intercept + post$b_NK * NK + (post$b_z +  postG$b_z  ) * z + post$b_z2 * z2

  return(mu)
  
}


postG3 = posterior_samples(G3)


size = seq(from=5, to=32, length= 100)
KGg = sapply(1:length(size), function(i) M3_link(post = postG3, NK = 0, size = size[i]))
KGg = cbind(apply(KGg, 2, mean),
t(apply(KGg, 2, hdi))[,1],
t(apply(KGg, 2, hdi))[,2])

NKg = sapply(1:length(size), function(i) M3_link(post = postG3, NK = 1, size = size[i]))
NKg = cbind(apply(NKg, 2, mean),
            t(apply(NKg, 2, hdi))[,1],
            t(apply(NKg, 2, hdi))[,2])


dfG = as.data.frame(rbind(KGg,NKg ))
names(dfG) <- c("growth", "l95", "u95")
dfG$Treatment = c(rep("KG", dim(KGg)[1]), rep("NK", dim(NKg)[1]))
dfG$z = c(size, size)

dGpoints = subset(Gdata, surv == 1)
dGpoints$Treatment = ifelse(dGpoints$NK==1, "NK", "KG")
dGpoints = dGpoints[,c("SL1_mm", "growth", "Treatment")]

names(dGpoints)[1] = "z" 

#GrowthGup_plot
(plotB = ggplot(dfG, aes(x = z, y = growth, colour = Treatment)) + 
  geom_line(size = 1.5, show.legend = T) + 
  ylab("Growth ln(z1/z0)") +
  xlab("Initial size (mm)") +
  geom_point(data = dGpoints, aes(x=z, y=growth, colour= Treatment), alpha = 0.75, size = 1) +
  geom_ribbon( aes(ymin = l95, ymax = u95, fill= Treatment), 
                 alpha = 0.2, size = 0.01) + 
      theme_classic() + theme(panel.grid.major = element_blank(), legend.position = c(0.8,0.8)) +

  scale_color_manual(values=c("black", "orange")) +
  scale_fill_manual(values=c("black", "orange")) + geom_vline(xintercept=18, lty=2, size = 0.2) +
  xlim(10,30)
)

summary(G3)







## Guppy reproduction

R1 <- brm(Recr ~ z + z2 + NK + z:NK + z2:NK  + Density + canopy + (1|stream), family = negbinomial(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 11))
R1 <- add_criterion(R1, c("loo", "waic"))
## remove effects

R2 <- update(R1, formula. = ~ . -  z2:NK)
R2 <- add_criterion(R2, c("loo", "waic"))


R3 <- update(R2, formula. = ~ . -  z:NK)
R3 <- add_criterion(R3, c("loo", "waic"))


R4 <- update(R3, formula. = ~ . -  z2)
R4 <- add_criterion(R4, c("loo", "waic"))

R5 <- update(R4, formula. = ~ . -  z)
R5 <- add_criterion(R5, c("loo", "waic"))


as.data.frame(loo_compare(R1,R2,R3,R4,R5, criterion = "loo"))
loo_model_weights(R1,R2,R3,R4,R5, criterion = "loo")


postR2 = posterior_samples(R2)

summary(R2)


M2_link= function(post, NK, size){
  z = size - 18
  z2 = (size^2 - 18^2)
  
  mu = post$b_Intercept + post$b_NK * NK + (post$b_z  + post$`b_z:NK`*NK) * z + post$b_z2 * z2
  
  return(mu)
  
}



size = seq(from=5, to=32, length= 100)
KGr = sapply(1:length(size), function(i) M3_link(post = postR2, NK = 0, size = size[i]))
KGr = exp(KGr)
KGr = cbind(apply(KGr, 2, mean),
            t(apply(KGr, 2, hdi))[,1],
            t(apply(KGr, 2, hdi))[,2])

NKr = sapply(1:length(size), function(i) M3_link(post = postR2, NK = 1, size = size[i]))
NKr = exp(NKr)
NKr = cbind(apply(NKr, 2, mean),
            t(apply(NKr, 2, hdi))[,1],
            t(apply(NKr, 2, hdi))[,2])



dfR = as.data.frame(rbind(KGr,NKr))
names(dfR) <- c("Recr", "l95", "u95")
dfR$Treatment = c(rep("KG", dim(KGr)[1]), rep("NK", dim(NKr)[1]))
dfR$z = c(size, size)

dRpoints = subset(Gdata, surv == 1)
dRpoints$Treatment = ifelse(dRpoints$NK==1, "NK", "KG")
dRpoints = dRpoints[,c("SL1_mm", "Recr", "Treatment")]

names(dRpoints)[1] = "z" 

#GrowthGup_plot
(plotC = ggplot(dfR, aes(x = z, y = Recr, colour = Treatment)) + 
    geom_line(size = 1.5) + 
    ylab("Offspring (N)") +
    xlab("Initial size (mm)") +
    geom_point(data = dRpoints, aes(x=z, y=Recr, colour= Treatment), alpha = 0.75, size = 1) +
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = "none") +
    geom_ribbon( aes(ymin = l95, ymax = u95, fill= Treatment), 
                 alpha = 0.2, size = 0.01, show.legend = F) + 
    scale_color_manual(values=c("black", "orange")) +
    scale_fill_manual(values=c("black", "orange")) + geom_vline(xintercept=18, lty=2, size = 0.2) +
    xlim(10,30) + ylim(0,32)
)

summary(R2)

### Killifish 
# Model selection for Killifish ---------------------------------------------

# Survival model ----------------------------------------------------------

KS1 <- brm(surv ~ z + z2 + NG + z:NG + z2:NG  + Density + canopy + (1|stream), family = bernoulli(), Kdata,
           iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 10))
KS1 <- add_criterion(KS1, c("loo", "waic"))
summary(KS1)

## remove effects

KS2 <- update(KS1, formula. = ~ . -  z2:NG)
KS2 <- add_criterion(KS2, c("loo", "waic"))


KS3 <- update(KS2, formula. = ~ . -  z:NG)
KS3 <- add_criterion(KS3, c("loo", "waic"))


KS4 <- update(KS3, formula. = ~ . -  z2)
KS4 <- add_criterion(KS4, c("loo", "waic"))

(loo_compare(KS1,KS2,KS3,KS4, criterion = "loo"))
loo_model_weights(KS1,KS2,KS3,KS4, criterion = "loo")

postKS = data.frame(posterior_samples(KS3))


Kdata$Treatment = ifelse(Kdata$NG == 1, "NG", "KG" ) 

summary(KS3)

KM3_link= function(post, NG, size){
  z = size - 18
  z2 = (size^2) - 18^2
  
  mu = post$b_Intercept + post$b_NG * NG + post$b_z * z + post$b_z2 * z2
  
  return(mu)
  
}

postKS= posterior_samples(KS3)

summary(KS3)

size = seq(from=10, to=90, by= 0.1)

KKGs = 100* inv_logit_scaled(sapply(1:length(size), function(i) KM3_link(post = postKS, NG = 0, size = size[i])))
KKGs = cbind(apply(KKGs, 2, mean),
             t(apply(KKGs, 2, hdi))[,1],
             t(apply(KKGs, 2, hdi))[,2])

KNGs = 100* inv_logit_scaled(sapply(1:length(size), function(i) KM3_link(post = postKS, NG = 1, size = size[i])))
KNGs = cbind(apply(KNGs, 2, mean),
             t(apply(KNGs, 2, hdi))[,1],
             t(apply(KNGs, 2, hdi))[,2])



dfs = as.data.frame(rbind(KKGs,KNGs ))
names(dfs) <- c("surv", "l95", "u95")
dfs$Treatment = c(rep("KG", dim(KKGs)[1]), rep("NG", dim(KNGs)[1]))
dfs$z = c(size, size)

names(dfs)

Kdata$Treatment <- ifelse(Kdata$NG==1, "NG", "KG")

splot = Kdata[,c("SL1_mm", "Treatment", "surv")]
names(splot)[1] = "z"
splot$surv = splot$surv * 100
splot

library(latex2exp)

(plotKA = ggplot(dfs, aes(x = z, y = surv, colour = Treatment)) + 
    geom_line(size = 1.5) + 
    ylab("Survival (%)") +
    xlab(TeX("Initial size, $z_{0}$ (mm)")) +
    geom_point(data=splot, aes(x=z, y=surv, colour= Treatment), alpha = 0.25, size = 0.75) +
    ylim(0, 100) + theme_classic() + theme(panel.grid.major = element_blank(), legend.position = "none") +
    geom_ribbon( aes(ymin = l95, ymax = u95, fill= Treatment), 
                 alpha = 0.2, size = 0.01, show.legend = F) + 
    scale_color_manual(values=c("black", "orange")) +
    scale_fill_manual(values=c("black", "orange"))  + geom_vline(xintercept=18, lty=2, size=0.2)+
    xlim(10,90)
)


## Guppy growth
# Growth model ------------------------------------------------------------

#----------- picewise regression


bprior <-   prior(normal(0, 1), nlpar = "bz1") +
  prior(normal(0, 1), nlpar = "bz2") +
  prior(normal(20, 10), nlpar = "omega") + 
  prior(normal(0, 5), nlpar = "b0")

Kdata$zr = Kdata$SL1_mm 
range(Kdata$zr)

bform <- bf(
  growth ~ b0 + bz1 * (zr - omega) * step(omega - zr) + 
    bz2 * (zr - omega) * step(zr - omega),
  b0 ~  NG  + Density + canopy +  (1|stream),
  bz1 + bz2 ~ NG,
  omega ~ 1,
  nl = TRUE
)

#make_stancode(bform2, data = subset(Kdata, surv == 1), prior = bprior)
fit1 <- brm(bform, data = subset(Kdata, surv == 1), prior = bprior, chains = 4, save_pars = save_pars(all=TRUE))

fit1 <- add_criterion(fit1, c("loo", "waic"))
summary(fit1)
## remove effects




postG = posterior_samples(fit1)


size = seq(from=5, to=90, length = 200)
G2_linkK = function(post, z, NG){
  
  bz1= post$b_bz1_Intercept
  bz2= post$b_bz2_Intercept
  a= post$b_b0_Intercept
  bNG= post$"b_b0_NG"
  bz1NG= post$"b_bz1_NG"
  bz2NG= post$"b_bz2_NG"
  w = post$b_omega_Intercept
  
  
  #  for(j in 1:length(u)){
  if(z < mean(w)){
    u = a + bNG * NG + (bz1 + bz1NG * NG) * (z - w) # linear predictor
  }else{
    u = a + bNG* NG + (bz2 + bz2NG * NG) * (z - w) # linear predictor
    #   }
  }
  
  return(u)
  
}


size = seq(from=5, to=90, by = 0.2)
KG = sapply(1:length(size), function(i) G2_linkK(post = postG, z = size[i], NG = 0))
KGg = cbind(apply(KG, 2, mean),
            t(apply(KG, 2, hdi))[,1],
            t(apply(KG, 2, hdi))[,2])

NG =sapply(1:length(size), function(i) G2_linkK(post = postG, z = size[i], NG = 1))
NGg = cbind(apply(NG, 2, mean),
            t(apply(NG, 2, hdi))[,1],
            t(apply(NG, 2, hdi))[,2])



dfG = as.data.frame(rbind(KGg,NGg ))

head(dfG)
names(dfG) <- c("growth", "l95", "u95")
dfG$Treatment = c(rep("KG", dim(KGg)[1]), rep("NG", dim(NGg)[1]))
dfG$z = c(size, size)

dGpoints = subset(Kdata, surv == 1)
dGpoints$Treatment = ifelse(dGpoints$NG==1, "NG", "KG")
dGpoints = dGpoints[,c("SL1_mm", "growth", "Treatment")]

names(dGpoints)[1] = "z" 

#GrowthGup_plot
(plotKB = ggplot(dfG, aes(x = z, y = growth, colour = Treatment)) + 
    geom_line(size = 1.5, show.legend = T) + 
    ylab("Growth ln(z1/z0)") +
    xlab("Initial size (mm)") +
    geom_point(data = dGpoints, aes(x=z, y=growth, colour= Treatment), alpha = 0.75, size = 1) +
    geom_ribbon( aes(ymin = l95, ymax = u95, fill= Treatment), 
                 alpha = 0.2, size = 0.01) + 
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = c(0.8,0.8)) +
    
    scale_color_manual(values=c("black", "orange")) +
    scale_fill_manual(values=c("black", "orange")) + geom_vline(xintercept=18, lty=2, size = 0.2) +
    xlim(10,90)
)

summary(fit1)


## Killifish reproduction

KR1 <- brm(Recr ~ z + z2 + NG + z:NG + z2:NG  + Density + canopy + (1|stream), family = negbinomial(), subset(Kdata, surv == 1),
           iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 11))
KR1 <- add_criterion(KR1, c("loo", "waic"))
## remove effects

KR2 <- update(KR1, formula. = ~ . -  z2:NG)
KR2 <- add_criterion(KR2, c("loo", "waic"))


KR3 <- update(KR2, formula. = ~ . -  z:NG)
KR3 <- add_criterion(KR3, c("loo", "waic"))


KR4 <- update(KR3, formula. = ~ . -  z2)
KR4 <- add_criterion(KR4, c("loo", "waic"))



as.data.frame(loo_compare(KR1,KR2,KR3,KR4, criterion = "loo"))
loo_model_weights(KR1,KR2,KR3,KR4, criterion = "loo")

summary(KR4)

KM4_link= function(post, NG, size){
  z = size - 18
  
  mu = post$b_Intercept + post$b_NG * NG + post$b_z * z
  
  return(exp(mu))
  
}


postKR4 = posterior_samples(KR4)

size = seq(from=5, to=90, by= 0.5)
KKGr = sapply(1:length(size), function(i)  KM4_link(post = postKR4, NG = 0, size = size[i]))
KKGr = cbind(apply(KKGr, 2, mean),
             t(apply(KKGr, 2, hdi))[,1],
             t(apply(KKGr, 2, hdi))[,2])

KNGr = sapply(1:length(size), function(i) KM4_link(post = postR4, NG = 1, size = size[i]))
KNGr = cbind(apply(KNGr, 2, mean),
             t(apply(KNGr, 2, hdi))[,1],
             t(apply(KNGr, 2, hdi))[,2])


dfR = as.data.frame(rbind(KKGr,KNGr))
dim(dfR)
names(dfR) <- c("Recr", "l95", "u95")
dfR$Treatment = c(rep("KG", dim(KNGr)[1]), rep("NG", dim(KNGr)[1]))
dfR$z = c(size, size)

dRpoints = subset(Kdata, surv == 1 & Recr >= 0 )
dRpoints$Treatment = ifelse(dRpoints$NG==1, "NG", "KG")
dRpoints = dRpoints[,c("SL1_mm", "Recr", "Treatment")]


names(dRpoints)[1] = "z" 

#GrowthGup_plot
(plotKC = ggplot(dRpoints, aes(x = z, y = Recr, colour = Treatment)) + 
    geom_point(alpha = 0.75, size = 1) +
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = "none") +
    #geom_ribbon( aes(ymin = l95, ymax = u95, fill= Treatment), 
    #            alpha = 0.2, size = 0.01, show.legend = F) + 
    ylab("Offspring (N)") +
    xlab("Initial size (mm)") +
    scale_color_manual(values=c("black", "orange")) +
    scale_fill_manual(values=c("black", "orange")) + geom_vline(xintercept=18, lty=2, size = 0.2) +
    xlim(10,90) + ylim(0,60)
)








library("ggpubr")
theme_set(
  theme_classic() 
  
)



figure1 <- ggarrange(NULL, plotA, NULL, plotB, NULL, plotC,
                    labels = c("A)", "", "B)","", "C)"),
                    ncol = 1, nrow = 6, heights = c(0.1, 1, 0.1, 1, 0.1, 1),
                    font.label = list(size = 10, color = "black", face = "bold", family = NULL)
)
figure1
svg("plots/Figure1.svg", width = 4, height = 7 )
figure1
graphics.off()


###-------------

theme_set(
  theme_classic() 
  
)


figure2 <- ggarrange(NULL, plotKA, NULL, plotKB, NULL, plotKC,
                     labels = c("A)", "", "B)","", "C)"),
                     ncol = 1, nrow = 6, heights = c(0.1, 1, 0.1, 1, 0.1, 1),
                     font.label = list(size = 10, color = "black", face = "bold", family = NULL)
)




svg("plots/Figure2.svg", width = 4, height = 7 )
figure2
graphics.off()