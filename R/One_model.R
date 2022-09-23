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


Gdata$KR = Gdata$NK
Kdata$GR = Kdata$NG

# Model selection for Guppies ---------------------------------------------

# Survival model ----------------------------------------------------------
# Guppies
S1 <- brm(surv ~ z + z2 + KR + z:KR   + Density + canopy + (1|stream), family = bernoulli(), Gdata,
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 10))

postS = data.frame(posterior_samples(S1))

predS = posterior_predict(S1, newdata = Gdata)

posterior

ppS = pp_check(S1, ndraws = 200)

svg("plots/PPSurvival.svg", width = 4, height = 4)
ppS
graphics.off()

Gdata$Treatment = ifelse(Gdata$NK == 1, "KR", "KG" ) 

P_link= function(post, KR, size){
  z = size - 18
  z2 = (size^2) - 18^2
  
  mu = post$b_Intercept + post$ b_KR  * KR + (post$b_z + post$`b_z:KR`*KR) * z + post$b_z2*z2
  
  return(mu)
  
}

P_link(postS, KR = 0,size = 5)

summary(S1)
postS = posterior_samples(S1)


round(cbind((apply(postS, 2, mean)), t(apply(postS, 2, hdi)),(apply(postS, 2, LOS))),3)

#
mean(inv_logit_scaled(P_link(post = postS, KR = 0, size = 10)))*100
hdi(inv_logit_scaled(P_link(post = postS, KR = 0, size = 10)))*100

mean(inv_logit_scaled(P_link(post = postS, KR = 1, size = 10)))*100
hdi(inv_logit_scaled(P_link(post = postS, KR = 1, size = 10)))*100




#

size = seq(from=5, to=32, length= 100)
NGs = 100* inv_logit_scaled(sapply(1:length(size), function(i) P_link(post = postS, KR = 0, size = size[i])))
NKs = 100* inv_logit_scaled(sapply(1:length(size), function(i) P_link(post = postS, KR = 1, size = size[i])))

pp = paste(round(LOS(NKs-NGs),1),"%", sep="")


NGs = cbind(apply(NGs, 2, mean),
            t(apply(NGs, 2, hdi))[,1],
            t(apply(NGs, 2, hdi))[,2])

NKs = cbind(apply(NKs, 2, mean),
            t(apply(NKs, 2, hdi))[,1],
            t(apply(NKs, 2, hdi))[,2])



dfs = as.data.frame(rbind(NGs,NKs ))
names(dfs) <- c("surv", "l95", "u95")
dfs$Treatment = c(rep("KG", dim(NKs)[1]), rep("KR", dim(NKs)[1]))
dfs$z = c(size, size)

names(dfs)

Gdata$Treatment <- ifelse(Gdata$NK==1, "KR", "KG")

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
    xlim(10,30) + 
    annotate("text", x=28, y=10, label= pp)
  
)

plotA
## Guppy growth

G1 <- brm(growth ~ z + z2 + KR + z:KR + Density + canopy + (1|stream), family = gaussian(), subset(Gdata, surv == 1), 
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 11))
G1 <- add_criterion(G1, c("loo", "waic"))

summary(G1)

postG = posterior_samples(G1)

round(cbind((apply(postG, 2, mean)), t(apply(postG, 2, hdi)),(apply(postG, 2, LOS))),3)
apply(postG, 2, hdi)
##
exp(mean((P_link(post = postG, KR = 0, size = 10))))
exp(hdi((P_link(post = postG, KR = 0, size = 10))))


exp(mean((P_link(post = postG, KR = 1, size = 10))))
exp(hdi((P_link(post = postG, KR = 1, size = 10))))




size = seq(from=5, to=32, length= 100)
KGg = sapply(1:length(size), function(i) P_link(post = postG, KR = 0, size = size[i]))
NKg = sapply(1:length(size), function(i) P_link(post = postG, KR = 1, size = size[i]))


pp = paste(round(LOS(NKg-KGg),1),"%", sep="")


KGg = cbind(apply(KGg, 2, mean),
            t(apply(KGg, 2, hdi))[,1],
            t(apply(KGg, 2, hdi))[,2])

NKg = cbind(apply(NKg, 2, mean),
            t(apply(NKg, 2, hdi))[,1],
            t(apply(NKg, 2, hdi))[,2])


dfG = as.data.frame(rbind(KGg,NKg ))
names(dfG) <- c("growth", "l95", "u95")
dfG$Treatment = c(rep("KG", dim(KGg)[1]), rep("KR", dim(NKg)[1]))
dfG$z = c(size, size)

dGpoints = subset(Gdata, surv == 1)
dGpoints$Treatment = ifelse(dGpoints$NK==1, "KR", "KG")
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
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = c(0.8,0.75)) +
    
    scale_color_manual(values=c("black", "orange")) +
    scale_fill_manual(values=c("black", "orange")) + geom_vline(xintercept=18, lty=2, size = 0.2) +
    xlim(10,30) + annotate("text", x = 13, y=.75, label = pp )
)

summary(G1)

plotB

ppG = pp_check(G1, ndraws = 200)

svg("plots/PPGrowthG.svg", width = 4, height = 4)
ppG
graphics.off()






## Guppy reproduction

R1 <- brm(Recr ~ z + z2 + KR + z:KR + Density + canopy + (1|stream), family = negbinomial(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 11))
postR = posterior_samples(R1)

apply(postR, 2, LOS)
apply(postR, 2, mean)
summary(R1)


round(cbind((apply(postR, 2, mean)), t(apply(postR, 2, hdi)),(apply(postR, 2, LOS))),3)

size = seq(from=5, to=32, length= 100)
KGr = sapply(1:length(size), function(i) P_link(post = postR, KR = 0, size = size[i]))
NKr = sapply(1:length(size), function(i) P_link(post = postR, KR = 1, size = size[i]))
pp = paste(round(LOS(NKr-KGr),1),"%", sep="")

KGr = exp(KGr)
KGr = cbind(apply(KGr, 2, mean),
            t(apply(KGr, 2, hdi))[,1],
            t(apply(KGr, 2, hdi))[,2])




NKr = exp(NKr)
NKr = cbind(apply(NKr, 2, mean),
            t(apply(NKr, 2, hdi))[,1],
            t(apply(NKr, 2, hdi))[,2])


dfR = as.data.frame(rbind(KGr,NKr))
names(dfR) <- c("Recr", "l95", "u95")
dfR$Treatment = c(rep("KG", dim(KGr)[1]), rep("KR", dim(NKr)[1]))
dfR$z = c(size, size)

dRpoints = subset(Gdata, surv == 1)
dRpoints$Treatment = ifelse(dRpoints$NK==1, "KR", "KG")
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
    xlim(10,30) + ylim(0,32) + annotate("text", x = 28, y=30, label = pp )
)
summary(R1)
plotC


(ppR = pp_check(R1, ndraws = 200))

svg("plots/PPRecproG.svg", width = 4, height = 4)
ppR
graphics.off()

### Killifish 
# Model Killifish ---------------------------------------------

KP_link= function(post,GR, size){
  z = size - 18
  z2 = (size^2) - 18^2
  
  mu = post$b_Intercept + post$b_NG *GR + (post$b_z + post$`b_z:GR` * GR) * z + post$b_z2*z2
  
  return(mu)
  
}

# Survival model ----------------------------------------------------------

KS1 <- brm(surv ~ z + z2 +GR + z:GR  + Density + canopy + (1|stream), family = bernoulli(), Kdata,
           iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 10))
summary(KS1)

(ppSK = pp_check(KS1, ndraws = 200))

svg("plots/PPSurvK.svg", width = 4, height = 4)
ppSK
graphics.off()


postKS = data.frame(posterior_samples(KS1))
Kdata$Treatment = ifelse(Kdata$NG == 1,  "GR", "KG" ) 


postKS= posterior_samples(KS1)

round(cbind((apply(postKS, 2, mean)), t(apply(postKS, 2, hdi)),(apply(postKS, 2, LOS))),3)

size = seq(from=10, to=90, by= 0.1)



KKGs = 100* inv_logit_scaled(sapply(1:length(size), function(i) KP_link(post = postKS,GR = 0, size = size[i])))
KNGs = 100* inv_logit_scaled(sapply(1:length(size), function(i) KP_link(post = postKS,GR = 1, size = size[i])))

pp = paste(round(LOS(KNGs-KKGs),1),"%", sep="")

KKGs = cbind(apply(KKGs, 2, mean),
             t(apply(KKGs, 2, hdi))[,1],
             t(apply(KKGs, 2, hdi))[,2])

KNGs = cbind(apply(KNGs, 2, mean),
             t(apply(KNGs, 2, hdi))[,1],
             t(apply(KNGs, 2, hdi))[,2])



dfs = as.data.frame(rbind(KKGs,KNGs ))
names(dfs) <- c("surv", "l95", "u95")
dfs$Treatment = c(rep("KG", dim(KKGs)[1]), rep( "GR", dim(KNGs)[1]))
dfs$z = c(size, size)

names(dfs)

Kdata$Treatment <- ifelse(Kdata$NG==1,  "GR", "KG")

splot = Kdata[,c("SL1_mm", "Treatment", "surv")]
names(splot)[1] = "z"
splot$surv = splot$surv * 100
splot

(plotKA = ggplot(dfs, aes(x = z, y = surv, colour = Treatment)) + 
    geom_line(size = 1.5) + 
    ylab("Survival (%)") +
    xlab("Initial size (mm)") +
    geom_point(data=splot, aes(x=z, y=surv, colour= Treatment), alpha = 0.25, size = 0.75) +
    ylim(0, 100) + theme_classic() + theme(panel.grid.major = element_blank(), legend.position = "none") +
    geom_ribbon( aes(ymin = l95, ymax = u95, fill= Treatment), 
                 alpha = 0.2, size = 0.01, show.legend = F) + 
    scale_color_manual(values=c("black", "orange")) +
    scale_fill_manual(values=c("black", "orange"))  + geom_vline(xintercept=18, lty=2, size=0.2)+
    xlim(10,90) + annotate("text", x = 88, y=90, label = pp )
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
  b0 ~ GR  + Density + canopy +  (1|stream),
  bz1 + bz2 ~GR,
  omega ~ 1,
  nl = TRUE
)

#make_stancode(bform2, data = subset(Kdata, surv == 1), prior = bprior)
fit1 <- brm(bform, data = subset(Kdata, surv == 1), prior = bprior, chains = 4, save_pars = save_pars(all=TRUE))

fit1 <- add_criterion(fit1, c("loo", "waic"))
summary(fit1)
## remove effects


(ppGK = pp_check(fit1, ndraws = 200))

svg("plots/PPSurvK.svg", width = 4, height = 4)
ppSK
graphics.off()


postG = posterior_samples(fit1)

round(cbind((apply(postG, 2, mean)), t(apply(postG, 2, hdi)),(apply(postG, 2, LOS))),3)


size = seq(from=5, to=90, length = 200)
G2_linkK = function(post, z,GR){
  
  bz1= post$b_bz1_Intercept
  bz2= post$b_bz2_Intercept
  a= post$b_b0_Intercept
  bNG= post$"b_b0_NG"
  bz1NG= post$"b_bz1_NG"
  bz2NG= post$"b_bz2_NG"
  w = post$b_omega_Intercept
  
  
  #  for(j in 1:length(u)){
  if(z < mean(w)){
    u = a + bNG *GR + (bz1 + bz1NG *GR) * (z - w) # linear predictor
  }else{
    u = a + bNG*GR + (bz2 + bz2NG *GR) * (z - w) # linear predictor
    #   }
  }
  
  return(u)
  
}


size = seq(from=5, to=90, by = 0.2)
KG = sapply(1:length(size), function(i) G2_linkK(post = postG, z = size[i],GR = 0))

NG =sapply(1:length(size), function(i) G2_linkK(post = postG, z = size[i],GR = 1))
pp = paste(round(LOS(KG-NG),1),"%", sep="")

KGg = cbind(apply(KG, 2, mean),
            t(apply(KG, 2, hdi))[,1],
            t(apply(KG, 2, hdi))[,2])



NGg = cbind(apply(NG, 2, mean),
            t(apply(NG, 2, hdi))[,1],
            t(apply(NG, 2, hdi))[,2])



dfG = as.data.frame(rbind(KGg,NGg ))

head(dfG)
names(dfG) <- c("growth", "l95", "u95")
dfG$Treatment = c(rep("KG", dim(KGg)[1]), rep( "GR", dim(NGg)[1]))
dfG$z = c(size, size)

dGpoints = subset(Kdata, surv == 1)
dGpoints$Treatment = ifelse(dGpoints$NG==1,  "GR", "KG")
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
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = c(0.25,0.85)) +
    
    scale_color_manual(values=c("black", "orange")) +
    scale_fill_manual(values=c("black", "orange")) + geom_vline(xintercept=29.67, lty=2, size = 0.2) +
    xlim(10,90) + annotate("text", x = 88, y=0.65, label = pp )
)

summary(fit1)


## Killifish reproduction

KR1 <- brm(Recr ~ z + z2 +GR + z:GR + Density + canopy + (1|stream), family = negbinomial(), subset(Kdata, surv == 1),
           iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 11))

postKR = posterior_samples(KR1)

round(cbind((apply(postKR, 2, mean)), t(apply(postKR, 2, hdi)),(apply(postKR, 2, LOS))),3)


size = seq(from=5, to=90, by= 0.5)
KKGr = sapply(1:length(size), function(i)  KP_link(post = postKR,GR = 0, size = size[i]))
KKGr = cbind(apply(KKGr, 2, mean),
             t(apply(KKGr, 2, hdi))[,1],
             t(apply(KKGr, 2, hdi))[,2])

KNGr = sapply(1:length(size), function(i)  KP_link(post = postKR,GR = 1, size = size[i]))
KNGr = cbind(apply(KNGr, 2, mean),
             t(apply(KNGr, 2, hdi))[,1],
             t(apply(KNGr, 2, hdi))[,2])


dfR = as.data.frame(rbind(KKGr,KNGr))
dim(dfR)
names(dfR) <- c("Recr", "l95", "u95")
dfR$Treatment = c(rep("KG", dim(KNGr)[1]), rep( "GR", dim(KNGr)[1]))
dfR$z = c(size, size)

dRpoints = subset(Kdata, surv == 1 & Recr >= 0 )
dRpoints$Treatment = ifelse(dRpoints$NG==1,  "GR", "KG")
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


## Recruitment


RecrData <- read.csv("data/Recruitment_data.csv")

head(RecrData) 

RecrData$Pool = factor(paste(RecrData$Removal,RecrData$Sp, sep = "-"))

levels(RecrData$Pool) <- c("KG", "KG", "KR",  "GR")
RecrData$Pool = factor(RecrData$Pool, levels = c("KG", "KR",  "GR"))

# get_prior(Recrt ~ 0 + Pool + (1|Location), family = negbinomial(), Data)


priors2 = prior(normal(0,1), class = b, coef = PoolNK) + # +
  prior(normal(3,0.2), class = Intercept) #+
# prior(normal(3,0.5), class = b, coef = PoolNG) +
# prior(normal(3,0.5), class = b, coef = PoolNK) 

RectG = subset(RecrData, Sp=="Guppy")

get_prior(Recrt ~ Pool + Density + Canopy , family = negbinomial(), RectG)
mG <- brm( Recrt ~ Pool + Density + Canopy , family = negbinomial(), RectG, prior = priors2,
           iter = 2000, warmup = 1000, control = list(adapt_delta = 0.98, max_treedepth = 13))
post_Recr = posterior_samples(mG)
summary(mG)

t(rbind(apply(post_Recr, 2, mean),apply(post_Recr, 2, hdi),apply(post_Recr, 2, LOS)))

conditional_effects(mG, effects = "Pool")

KGrec = exp(post_Recr$b_Intercept)
NKrec = exp(post_Recr$b_Intercept + post_Recr$b_PoolNK)



RecrDF = as.data.frame(cbind(apply(cbind(KGrec,NKrec), 2, mean),
                             t(apply(cbind(KGrec,NKrec), 2, hdi, credMass=.95))[,1],
                             t(apply(cbind(KGrec,NKrec), 2, hdi, credMass=.95))[,2]))

RecrDF$Pool <- c("KG", "KR")
names(RecrDF)[1:3] <- c("mean", "l95", "u95" )

RecrDF

pp = LOS(NKrec - KGrec)
pp = paste(round(pp,1),"%", sep = "")






RecrData$Treatment =  RecrData$Pool
levels(RecrData$Treatment) <- c("KG", "KR", "KG",  "GR")

(plotD <- ggplot(RecrDF, aes(x=Pool, y=mean, colour = Pool)) + 
    geom_point(size = 4) + scale_color_manual(values=c("black", "orange")) +
    geom_errorbar(aes(ymin=l95, ymax=u95), width=.1) +
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = "none") +
    geom_point(data = subset(RecrData, Sp == "Guppy"), 
               mapping = aes(x = Treatment, y = Recrt, color = Treatment),
               position = position_jitterdodge(dodge.width=1.5),
               size = 1, alpha = 0.7) +
    ylab("Recruits (N)") +
    xlab("Treatment") + annotate("text", x = 2.5, y=85, label = pp )
  
)



#### Killifish


RectK = subset(RecrData, Sp=="Killifish")

RectK$Treatment = factor(RectK$Pool)

get_prior(Recrt ~ Pool + Density + Canopy , family = negbinomial(), RectK)

priors3 = prior(normal(0,0.5), class = b, coef = PoolNG) + # +
  prior(normal(2.5,0.5), class = Intercept) #+
# prior(normal(3,0.5), class = b, coef = PoolNG) +
# prior(normal(3,0.5), class = b, coef = PoolNK) 

log(mean(RectK$Recrt))

mK <- brm( Recrt ~ Pool + Density + Canopy , family = negbinomial(), RectK, prior = priors3,
           iter = 2000, warmup = 1000, control = list(adapt_delta = 0.98, max_treedepth = 13))
post_RecrK = posterior_samples(mK)
summary(mK)

t(rbind(apply(post_RecrK, 2, mean),apply(post_RecrK, 2, hdi),apply(post_RecrK, 2, LOS)))

conditional_effects(mK, effects = "Pool")

KKGrec = exp(post_RecrK$b_Intercept)
KNGrec = exp(post_RecrK$b_Intercept + post_RecrK$b_PoolNG)



RecrDF = as.data.frame(cbind(apply(cbind(KKGrec,KNGrec), 2, mean),
                             t(apply(cbind(KKGrec,KNGrec), 2, hdi, credMass=.95))[,1],
                             t(apply(cbind(KKGrec,KNGrec), 2, hdi, credMass=.95))[,2]))

RecrDF$Pool <- c("KG",  "GR")
names(RecrDF)[1:3] <- c("mean", "l95", "u95" )

RecrDF

ppk = LOS(KNGrec - KKGrec)
(ppk = paste(round(ppk,1),"%", sep = ""))

(plotDK <- ggplot(RecrDF, aes(x=Pool, y=mean, colour = Pool)) + 
    geom_point(size = 4) + scale_color_manual(values=c("black", "orange")) +
    geom_errorbar(aes(ymin=l95, ymax=u95), width=.1) +
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = "none") +
    geom_point(data = RectK, 
               mapping = aes(x = Treatment, y = Recrt, color = Treatment),
               position = position_jitterdodge(dodge.width=1.5),
               size = 1, alpha = 0.7) +
    ylab("Recruits (N)") +
    xlab("Treatment") + annotate("text", x = 2.5, y=60, label = ppk )
)

###
library("ggpubr")
theme_set(
  theme_classic() 
  
)



figure1 <- ggarrange(plotA, plotB, plotC, plotD,
                     labels = c("A)", "B)", "C)", "D)"),
                     ncol = 2, nrow = 2, 
                     #heights = c(0.9, 0.9, 0.9, 0.9),
                     #withds =  c(0.9, 0.9, 0.9, 0.9),
                     font.label = list(size = 10, color = "black", face = "bold", family = NULL)
)

figure1

svg("plots/Figure1.svg", width = 7, height = 7 )
figure1
graphics.off()


###-------------

theme_set(
  theme_classic() 
  
)


figure2 <- ggarrange(plotKA, plotKB, plotKC, plotDK,
                     labels = c("A)", "B)", "C)", "D)"),
                     ncol = 2, nrow = 2, 
                     #heights = c(0.9, 0.9, 0.9, 0.9),
                     #withds =  c(0.9, 0.9, 0.9, 0.9),
                     font.label = list(size = 10, color = "black", face = "bold", family = NULL)
)
figure2

svg("plots/Figure2.svg", width = 7, height = 7 )
figure2
graphics.off()




