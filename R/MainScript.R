####  
library(data.table)
library(brms)
library(ggplot2)
library(HDInterval)
library("plyr")
library(dplyr)
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

### Use fish density as biomass per area
Gdata$Density <-( Gdata$BiomassG1 + Gdata$BiomassK1) / Gdata$area
Kdata$Density <- (Kdata$BiomassG1 + Kdata$BiomassK1) / Kdata$area


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


##
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



Gdata$Male = ifelse(Gdata$Sex1 == "M", 1, 0)

Kdata[-which(Kdata$Sex1 == Kdata$Sex2), c("Sex1", "Sex2")]
Kdata$Sex2 = factor(Kdata$Sex2)
levels(Kdata$Sex2)<- c("F", "J", "M")
levels(Kdata$Sex2)

Kdata[, c("Sex1", "Sex2")]

Kdata$Sex = Kdata$Sex1 
Kdata$Sex[-which(is.na(Kdata$Sex2) )] = as.character(Kdata$Sex2[-which(is.na(Kdata$Sex2))]) 
Kdata$Male = ifelse(Kdata$Sex == "M", 1, 0)

# Model selection for Guppies ---------------------------------------------

# Survival model ----------------------------------------------------------
# Guppies
S1 <- brm(surv ~ z*Male*KR + z2    + Density + canopy + (1|stream), family = bernoulli(), Gdata,
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 10))

postS = as.data.frame(as_draws_df(S1))

names(postS)

summary(S1)

(ppS = pp_check(S1, ndraws = 200))

svg("plots/PPSurvival.svg", width = 4, height = 4)
ppS
graphics.off()

Gdata$Treatment = ifelse(Gdata$NK == 1, "KR", "KG" ) 



post= postS
KR = 0; size = 18; male = 0

names(post)
P_link= function(post, KR = 0, size = 18, male = 0, center = center){

  z = size - center
  z2 = (size^2) - center^2
  

  mu = post$b_Intercept + post$b_KR  * KR + (post$b_z + post$`b_z:KR`*KR) * z + post$b_z2*z2 +
      post$b_Male * male + post$`b_Male:KR` *(male*KR) + post$`b_z:Male` *(male*z) + post$`b_z:Male:KR` * (z*male*KR)    
  return(mu)
  
}

P_link(post = postS,KR = 0, size = 10, male = 0, center)

summary(S1)


round(cbind((apply(postS, 2, mean)), t(apply(postS, 2, hdi)),(apply(postS, 2, LOS))),3)

#
mean(inv_logit_scaled(P_link(post = postS, KR = 0, size = 10, male = 0, center = center)))*100
hdi(inv_logit_scaled(P_link(post = postS, KR = 0, size = 10, center = center)))*100

mean(inv_logit_scaled(P_link(post = postS, KR = 1, size = 10, center= center )))*100
hdi(inv_logit_scaled(P_link(post = postS, KR = 1, size = 10, center= center)))*100




#
HPDI68 <- function(x){
  hdi(x, credMass = .68)
}


sizeF = seq(from=5, to=32, length= 100)


NGsF = 100* inv_logit_scaled(sapply(1:length(sizeF), function(i) P_link(post = postS, KR = 0, size = sizeF[i], center= center)))
NKsF = 100* inv_logit_scaled(sapply(1:length(sizeF), function(i) P_link(post = postS, KR = 1, size = sizeF[i], center= center)))

sizeM = seq(from=5, to=max(Gdata$SL1_mm[which(Gdata$Male == 1)], na.rm = T), length= 100)
NGsM = 100* inv_logit_scaled(sapply(1:length(sizeM), function(i) P_link(post = postS, KR = 0, size = sizeM[i], male = 1,center= center)))
NKsM = 100* inv_logit_scaled(sapply(1:length(sizeM), function(i) P_link(post = postS, KR = 1, size = sizeM[i], male = 1, center= center)))


(ppF = paste(round(LOS(NKsF-NGsF),1),"%", sep=""))
(ppM = paste(round(LOS(NKsM-NGsM),1),"%", sep=""))

(ppFM = paste(round(LOS((NKsF +NGsF)/2 - (NKsM +NGsM)/2),1),"%", sep=""))




NGsF = cbind(apply(NGsF, 2, mean),
            t(apply(NGsF, 2, HPDI68))[,1],
            t(apply(NGsF, 2, HPDI68))[,2])

NKsF = cbind(apply(NKsF, 2, mean),
            t(apply(NKsF, 2, HPDI68))[,1],
            t(apply(NKsF, 2, HPDI68))[,2])


NGsM = cbind(apply(NGsM, 2, mean),
             t(apply(NGsM, 2, HPDI68))[,1],
             t(apply(NGsM, 2, HPDI68))[,2])

NKsM = cbind(apply(NKsM, 2, mean),
             t(apply(NKsM, 2, HPDI68))[,1],
             t(apply(NKsM, 2, HPDI68))[,2])



dfs = as.data.frame(rbind(NGsF,NKsF,NGsM,NKsM))
names(dfs) <- c("surv", "LC", "UC")
dfs$Pool = c(rep("KG", dim(NKsF)[1]), rep("KR", dim(NKsF)[1]), rep("KG", dim(NKsF)[1]), rep("KR", dim(NKsF)[1]))
dfs$Sex = c(rep("F", dim(NKsF)[1]), rep("F", dim(NKsF)[1]), rep("M", dim(NKsF)[1]), rep("M", dim(NKsF)[1]))
dim(dfs)
dfs$z = c(sizeF, sizeF,sizeM, sizeM)
dfs$Treatment = factor(paste(dfs$Sex, dfs$Pool, sep="-"), levels= c("M-KG","F-KG", "M-KR", "F-KR"))

names(dfs)

Gdata <- Gdata %>% mutate(Pool = factor(ifelse(NK ==1, "KR", "KG"), levels = c("KG", "KR")),
                 Sex = factor(ifelse(Sex1=="M", "M", "F"), levels = c("F", "M")),
                 Treatment =  factor(paste(Sex, Pool, sep="-"), levels= c("M-KG","F-KG", "M-KR", "F-KR"))) 
  


splot = Gdata %>% select(SL1_mm, Treatment, surv, Sex) %>%
     rename(z = SL1_mm) %>% mutate(surv = surv*100)


library(latex2exp)
(plotA = ggplot(dfs, aes(x = z, y = surv, group = Treatment)) + 
    geom_line(size = 0.5,  aes(linetype=Treatment, color=Treatment), show.legend = F) + 
    ylab("Survival (%)") +
    xlab("Initial size (mm)") +
    geom_point(data = splot, aes(x=z, y=surv, fill= Treatment, shape = Treatment),  show.legend = F, alpha = 0.75, size = 1.75) +
    geom_ribbon( aes(ymin = LC, ymax = UC, fill= Treatment),  show.legend = F, alpha = 0.2, size = 0.01) + 
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = c(0.8,0.75), legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(values=c(2,1,2,1)) +
    scale_shape_manual(values = c(21,24,21,24))+
    scale_color_manual(values=c("black", "black","orange", "orange")) +
    scale_fill_manual(values=c("black","black", "orange", "orange"))  + geom_vline(xintercept=18, lty=2, size = 0.2) +
    xlim(10,30) + ylim(0,101)  #scale_y_discrete(limits= (seq(0, 100, by= 25))) + 
#    annotate("text", x = c(13,20, 25), y=c(106, 106, 106), label = c(TeX("$\\F_{KR>KG}$: 50.3%"), TeX("$\\M_{KR>KG}$: 33.4%"),  TeX("F-M: 80.1%")) ) 
)
ppF
ppM
ppFM
## Guppy growth

G1 <- brm(growth ~ z*Male*KR + z2  + Density + canopy + (1|stream), family = gaussian(), subset(Gdata, surv == 1), 
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 11))
G1 <- add_criterion(G1, c("loo", "waic"))

summary(G1)

postG = as.data.frame(as_draws_df(G1))

round(cbind((apply(postG, 2, mean)), t(apply(postG, 2, hdi)),(apply(postG, 2, LOS))),3)
apply(postG, 2, hdi)
##
exp(mean((P_link(post = postG, KR = 0, size = 10, male = 0, center = center))))
exp(hdi((P_link(post = postG, KR = 0, size = 10))))


exp(mean((P_link(post = postG, KR = 1, size = 10))))
exp(hdi((P_link(post = postG, KR = 1, size = 10))))


# Plot growth---------------
dim(NGgF)

KGgF =  (sapply(1:length(sizeF), function(i) P_link(post = postG, KR = 0, size = sizeF[i], center= center, male = 0)))
NKgF =  (sapply(1:length(sizeF), function(i) P_link(post = postG, KR = 1, size = sizeF[i], center= center, male = 0)))
KGgM =  (sapply(1:length(sizeM), function(i) P_link(post = postG, KR = 0, size = sizeM[i], center= center, male = 1)))
NKgM =  (sapply(1:length(sizeM), function(i) P_link(post = postG, KR = 1, size = sizeM[i], center= center, male = 1)))


(ppF = paste(round(LOS(NKgF-KGgF),1),"%", sep=""))
(ppM = paste(round(LOS(NKgM-KGgM),1),"%", sep=""))
(ppM = paste(round(LOS((NKgF+KGgF) - (NKgM+KGgM) ),1),"%", sep=""))


KGgF = cbind(apply(KGgF, 2, mean),
             t(apply(KGgF, 2, HPDI68))[,1],
             t(apply(KGgF, 2, HPDI68))[,2])

NKgF = cbind(apply(NKgF, 2, mean),
             t(apply(NKgF, 2, HPDI68))[,1],
             t(apply(NKgF, 2, HPDI68))[,2])


KGgM = cbind(apply(KGgM, 2, mean),
             t(apply(KGgM, 2, HPDI68))[,1],
             t(apply(KGgM, 2, HPDI68))[,2])

NKgM = cbind(apply(NKgM, 2, mean),
             t(apply(NKgM, 2, HPDI68))[,1],
             t(apply(NKgM, 2, HPDI68))[,2])



dfg = as.data.frame(rbind(KGgF,NKgF,KGgM,NKgM))
names(dfg) <- c("growth", "LC", "UC")

dfg$Pool = c(rep("KG", dim(KGgF)[1]), rep("KR", dim(NKgF)[1]), rep("KG", dim(KGgF)[1]), rep("KR", dim(NKgF)[1]))
dfg$Sex = c(rep("F", dim(NKgF)[1]), rep("F", dim(KGgF)[1]), rep("M", dim(KGgM)[1]), rep("M", dim(NKgM)[1]))
dim(dfg)
dfg$z = c(sizeF, sizeF,sizeM, sizeM)
dfg$Treatment = factor(paste(dfg$Sex, dfg$Pool, sep="-"), levels= c("M-KG","F-KG", "M-KR", "F-KR"))
names(dfg)


dGpoints = Gdata %>% filter(surv == 1) %>%
  select(c(growth, SL1_mm, Treatment, Sex)) %>% 
  rename(z = SL1_mm) 


#GrowthGup_plot
(plotB = ggplot(dfg, aes(x = z, y = growth, group = Treatment)) + 
    geom_line(size = 0.5,  aes(linetype=Treatment, color=Treatment), show.legend = T) + 
    ylab("Growth ln(final size/initial size)") +
    xlab("Initial size (mm)") +
    geom_point(data = dGpoints, aes(x=z, y=growth, fill= Treatment, shape = Treatment),  show.legend = T, alpha = 0.75, size = 1.75) +
    geom_ribbon( aes(ymin = LC, ymax = UC, fill= Treatment),  show.legend = F, alpha = 0.2, size = 0.01) + 
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = c(0.8,0.75), legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(values=c(2,1,2,1)) +
    scale_shape_manual(values = c(21,24,21,24))+
    scale_color_manual(values=c("black", "black","orange", "orange")) +
    scale_fill_manual(values=c("black","black", "orange", "orange"))  + geom_vline(xintercept=18, lty=2, size = 0.2) +
    xlim(10,30) #+ 
  #  annotate("text", x = c(13,20, 25), y=c(.9, .9, .9), label = c(TeX("$\\F_{KR>KG}$: 59%"), TeX("$\\M_{KR>KG}$: 27.3%"),  TeX("F>M: 80.1%")) ) 
)
ppF
ppM
ppFM
summary(G1)

(ppG = pp_check(G1, ndraws = 200))

svg("plots/PPGrowthG.svg", width = 4, height = 4)
ppG
graphics.off()






## Guppy reproduction

R1 <- brm(Recr ~ z*KR  + Density + canopy + (1|stream), family = negbinomial(), subset(Gdata, surv == 1 & Sex != "M" ),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 11))


postR = as.data.frame(as_draws_df(R1))

apply(postR, 2, LOS)
apply(postR, 2, mean)
summary(R1)


round(cbind((apply(postR, 2, mean)), t(apply(postR, 2, hdi)),(apply(postR, 2, LOS))),3)


PR_link= function(post, KR = 0, size = 18, center = center){
  
  z = size - center
#  z2 = (size^2) - center^2
  
  
  mu = post$b_Intercept + post$b_KR  * KR + (post$b_z + post$`b_z:KR`*KR) * z #+ post$b_z2*z2
  
  return(mu)
  
}

size = seq(from=5, to=32, length= 100)
KGr = sapply(1:length(size), function(i) PR_link(post = postR, KR = 0, size = size[i], center = center ))
NKr = sapply(1:length(size), function(i) PR_link(post = postR, KR = 1, size = size[i], center = center ))
pp = paste(round(LOS(NKr-KGr),1),"%", sep="")

KGr = exp(KGr)
KGr = cbind(apply(KGr, 2, mean),
            t(apply(KGr, 2, HPDI68))[,1],
            t(apply(KGr, 2, HPDI68))[,2])




NKr = exp(NKr)
NKr = cbind(apply(NKr, 2, mean),
            t(apply(NKr, 2, HPDI68))[,1],
            t(apply(NKr, 2, HPDI68))[,2])


dfR = as.data.frame(rbind(KGr,NKr))
names(dfR) <- c("Recr", "LC", "UC")
dfR$Treatment = c(rep("KG", dim(KGr)[1]), rep("KR", dim(NKr)[1]))
dfR$z = c(size, size)

# dRpoints = subset(Gdata, surv == 1 )
# dRpoints$Treatment = ifelse(dRpoints$NK==1, "KR", "KG")
# dRpoints = dRpoints[,c("SL1_mm", "Recr", "Treatment")]

dRpoints = Gdata %>% filter(surv == 1, Sex != "M") %>%
  mutate(Treatment = factor(ifelse(Treatment == "F-KG", "KG", "KR"))) %>%
  select(c(Recr, SL1_mm, Treatment, Sex)) %>% 
  rename(z = SL1_mm) 

(plotC = ggplot(dfR, aes(x = z, y = Recr, group = Treatment)) + 
    geom_line(size = 1, aes(x=z, y=Recr, colour= Treatment )) + 
    ylab("Offspring (N)") +
    xlab("Initial size (mm)") +
    geom_point(data = dRpoints, aes(x=z, y=Recr,  fill= Treatment ), shape = 21, alpha = 0.75, size = 1.75) +
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = "none") +
    geom_ribbon( aes(ymin = LC, ymax = UC, fill= Treatment), 
                 alpha = 0.2, size = 0.01, show.legend = F) + 
    scale_color_manual(values=c("black", "orange")) +
    #scale_shape_manual(values = c(21,21))+
    scale_fill_manual(values=c("black", "orange")) + geom_vline(xintercept=18, lty=2, size = 0.2) +
    xlim(10,30) + ylim(0,32) #+ annotate("text", x = 27, y=30, label = pp ) + 
#    geom_tile()
  
)
summary(R1)
plotC


# plot(apply(posterior_predict(R1), 2, mean) ~ dRpoints$Recr)
# abline(a = 0, b = 1)
# 
# svg("plots/PPRecproG.svg", width = 4, height = 4)
# ppR
# graphics.off()



my_summary <- function(model, variable, species){
  
  df_sum = data.frame(Mean =apply(as_draws_df(model), 2, mean),
                      lC=t(apply(as_draws_df(model), 2, hdi))[,1],
                      uC=t(apply(as_draws_df(model), 2, hdi))[,2],
                      PP =apply(as_draws_df(model), 2, LOS)
                      )
  
  
  df_sum  <-round(df_sum, 3)
  df_sum <- df_sum[-((dim(df_sum)[1]-5):dim(df_sum)[1]),]
  df_sum$var <- variable 
  df_sum$sp <- species
  
  return(df_sum)
  
  
  
}

SumTabG <- rbind(my_summary(S1, variable = "Survival", species = "Guppy"), 
                 my_summary(G1,variable = "Growth", species = "Guppy"), 
                 my_summary(R1,variable = "Fecundity", species = "Guppy"))




SumTabG$Sig <- ifelse(SumTabG$PP < 2.5 | SumTabG$PP > 97.5, "*", "")

SumTabG


# make summary table for guppies here:

write.csv(SumTabG, "outputs/Table-S3.csv")





### Killifish 
# Model Killifish ---------------------------------------------

names(postKS)
KP_link= function(post,GR, size, male, center){
  
  z = size - center
  z2 = (size^2) - center^2
  
  if(which(names(post) == "b_Male:GR") >2 ){
    mu = post$b_Intercept + post$b_GR *GR + (post$b_z + post$`b_z:GR` * GR) * z + post$b_z2*z2 + 
      post$b_Male * male + post$`b_Male:GR` *(male*GR) + post$`b_z:Male` *(male*z) + post$`b_z:Male:GR` * (z*male*GR)  
  }else{
    mu = post$b_Intercept + post$b_GR *GR + (post$b_z + post$`b_z.GR` * GR) * z + post$b_z2*z2 + 
      post$b_Male * male + post$`b_Male.GR` *(male*GR) + post$`b_z.Male` *(male*z) + post$`b_z.Male.GR` * (z*male*GR)  
    
    
    }
  
  return(mu)
  
}

# Survival model ----------------------------------------------------------

KS1 <- brm(surv ~ z*Male*GR + z2   + Density + canopy + (1|stream), family = bernoulli(), Kdata,
           iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 10))


(ppSK = pp_check(KS1, ndraws = 200))

svg("plots/PPSurvK.svg", width = 4, height = 4)
ppSK
graphics.off()


postKS = data.frame(as_draws_df(KS1))
names(postKS)
post_names = paste("b_",rownames(summary(KS1)$fixed), sep = "")
names(postKS)[1:length(post_names)] <- post_names


Kdata$Pool = ifelse(Kdata$NG == 1,  "GR", "KG" ) 

my_summary(KS1, variable = "surv", species = "kill")
size = seq(from=10, to=90, by= 0.1)



KKGsF = 100* inv_logit_scaled(sapply(1:length(size), function(i) KP_link(post = postKS,GR = 0, size = size[i], male = 0, center = center )))
KNGsF = 100* inv_logit_scaled(sapply(1:length(size), function(i) KP_link(post = postKS,GR = 1, size = size[i], male = 0, center = center)))


KKGsM = 100* inv_logit_scaled(sapply(1:length(size), function(i) KP_link(post = postKS,GR = 0, size = size[i], male = 1, center = center )))
KNGsM = 100* inv_logit_scaled(sapply(1:length(size), function(i) KP_link(post = postKS,GR = 1, size = size[i], male = 1, center = center)))

(ppF = paste(round(LOS(KNGsF-KKGsF),1),"%", sep=""))
(pM = paste(round(LOS(KNGsM-KKGsM),1),"%", sep=""))
(pFM = paste(round(LOS( (KNGsF+KKGsF)/2 - (KNGsM+KKGsM)/2    ),1),"%", sep=""))



KKGsF = cbind(apply(KKGsF, 2, median),
             t(apply(KKGsF, 2, HPDI68))[,1],
             t(apply(KKGsF, 2, HPDI68))[,2])

KNGsF = cbind(apply(KNGsF, 2, median),
             t(apply(KNGsF, 2, HPDI68))[,1],
             t(apply(KNGsF, 2, HPDI68))[,2])

KKGsM= cbind(apply(KKGsM, 2, median),
              t(apply(KKGsM, 2, HPDI68))[,1],
              t(apply(KKGsM, 2, HPDI68))[,2])

KNGsM = cbind(apply(KNGsM, 2, median),
              t(apply(KNGsM, 2, HPDI68))[,1],
              t(apply(KNGsM, 2, HPDI68))[,2])



dfs = as.data.frame(rbind(KKGsF,KNGsF,KKGsM,KNGsM ))
names(dfs) <- c("surv", "LC", "UC")
dfs$Pool = factor(c(rep("KG", dim(KKGsF)[1]), rep( "GR", dim(KNGsF)[1]),rep("KG", dim(KKGsM)[1]), rep( "GR", dim(KNGsM)[1])), levels = c("KG", "GR"))
dfs$Sex = factor(c(rep("F", dim(KKGsF)[1]), rep( "F", dim(KNGsF)[1]),rep("M", dim(KKGsM)[1]), rep( "M", dim(KNGsM)[1])), levels = c("F", "M"))
dfs$Treatment = factor(paste(dfs$Sex, dfs$Pool, sep = "-"), levels = c("F-KG", "M-KG", "F-GR", "M-GR"))
dfs$z = c(size, size,size, size)


Kdata$Treatment <- factor(paste(ifelse(Kdata$Sex== "M", "M", "F"), Kdata$Pool, sep = "-"), levels = c("F-KG", "M-KG", "F-GR", "M-GR"))

splot = Kdata[,c("SL1_mm", "Treatment", "surv")]
names(splot)[1] = "z"
splot$surv = splot$surv * 100
splot



(plotKA = ggplot(dfs, aes(x = z, y = surv, group = Treatment)) + 
    geom_line(size = 0.5,  aes(linetype=Treatment, color=Treatment), show.legend = F) + 
    ylab("Survival (%)") +
    xlab("Initial size (mm)") +
    geom_point(data = splot, aes(x=z, y=surv, fill= Treatment, shape = Treatment),  show.legend = F, alpha = 0.75, size = 1.75) +
    geom_ribbon( aes(ymin = LC, ymax = UC, fill= Treatment),  show.legend = F, alpha = 0.2, size = 0.01) + 
    theme_classic() + theme(panel.grid.major = element_blank())+
  scale_linetype_manual(values=c(2,1,2,1)) +
  scale_shape_manual(values = c(21,24,21,24))+
    scale_color_manual(values=c("black", "black","orange", "orange")) +
    scale_fill_manual(values=c("black","black", "orange", "orange"))  + geom_vline(xintercept=18, lty=2, size = 0.2) +
    xlim(10,90) + ylim(0,101)# + #scale_y_discrete(limits= (seq(0, 100, by= 25))) + 
   # annotate("text", x = c(25, 50, 75), y=c(106, 106, 106), label = c(TeX("$\\F_{KR>KG}$: 38.9%"), TeX("$\\M_{KR>KG}$: 27.3%"),  TeX("F>M: 80.1%")) )  
)

ppF
ppM
ppFM
#------------ Killifish growth
#----------- picewise regression


bprior <-   prior(normal(0, 5), nlpar = "bz1") +
  prior(normal(0, 5), nlpar = "bz2") +
  prior(normal(20, 10), nlpar = "omega") + 
  prior(normal(0, 10), nlpar = "b0")

Kdata$zr = Kdata$SL1_mm 
range(Kdata$zr)

bform <- bf(
  growth ~ b0 + bz1 * (zr - omega) * step(omega - zr) + 
    bz2 * (zr - omega) * step(zr - omega),
  b0 ~ GR*Male  + Density + canopy +  (1|stream),
  bz1 + bz2 ~ GR,
  omega ~ 1,
  nl = TRUE
)


#make_stancode(bform2, data = subset(Kdata, surv == 1), prior = bprior)
fit1 <- brm(bform, data = subset(Kdata, surv == 1), prior = bprior, chains = 4, save_pars = save_pars(all=TRUE),
            iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 12))

summary(fit1)

(ppGK = pp_check(fit1, ndraws = 200))

svg("plots/PPSurvK.svg", width = 4, height = 4)
ppSK
graphics.off()


postKG = as.data.frame(as_draws_df(fit1))



(a = summary(fit1))



ids = which(names(postKG) %in% c(".chain",".draw",".iteration") )


postKG = postKG[,  -ids]
postKG = postKG[,  sort(names(postKG))]

names(postKG)
fixedNames = sort(rownames(a$fixed))


names(postKG)[1:length(fixedNames)] <- fixedNames

my_summary(model = fit1, variable = "growth", species = "kill")


size = seq(from=5, to=90, length = 200)

G2_linkK = function(post, z,GR, male){
  
  
  a= post$b0_Intercept

  bGR= post$b0_GR 
  bMale = post$b0_Male
  bGRMale = post$`b0_GR:Male`
  bz1= post$bz1_Intercept
  bz1GR= post$bz1_GR
  # bz1Male= post$`bz1_GR:Male`
  # bz1GRMale= post$`bz1_GR:Male`

  bz2= post$bz2_Intercept
  bz2GR= post$bz2_GR
  # bz2Male= post$bz2_Male
  # bz2GRMale= post$`bz2_GR:Male` 
  
  w = post$omega_Intercept
  
  
  #  for(j in 1:length(u)){
  if(z < mean(w)){
    u = a + bGR *GR + bMale * male + bGRMale * (GR*male) + (bz1  + bz1GR *GR ) * (z - w) # linear predictor
  }else{
    u = a + bGR*GR + bMale * male + bGRMale * (GR*male) + (bz2 + bz2GR *GR) * (z - w) # linear predictor
  
  }
  
  return(u)
  
}


size = seq(from=5, to=90, by = 0.2)


KGf = sapply(1:length(size), function(i) G2_linkK(post = postKG, z = size[i],GR = 0, male = 0))
NGf =sapply(1:length(size), function(i) G2_linkK(post = postKG, z = size[i],GR = 1, male = 0))

KGm = sapply(1:length(size), function(i) G2_linkK(post = postKG, z = size[i],GR = 0, male = 1))
NGm =sapply(1:length(size), function(i) G2_linkK(post = postKG, z = size[i],GR = 1, male = 1))



(ppF = paste(round(LOS(KGf-NGf),1),"%", sep=""))
(ppM = paste(round(LOS(KGm-NGm),1),"%", sep=""))

(ppFM = paste(round(LOS((KGf + NGf)/2 -(KGf + NGf)/2),1),"%", sep=""))




KGgkf = cbind(apply(KGf, 2, mean),
            t(apply(KGf, 2, HPDI68))[,1],
            t(apply(KGf, 2, HPDI68))[,2])



NGgkf = cbind(apply(NGf, 2, mean),
            t(apply(NGf, 2, HPDI68))[,1],
            t(apply(NGf, 2, HPDI68))[,2])


KGgkm = cbind(apply(KGm, 2, mean),
              t(apply(KGm, 2, HPDI68))[,1],
              t(apply(KGm, 2, HPDI68))[,2])



NGgkm = cbind(apply(NGm, 2, mean),
             t(apply(NGm, 2, HPDI68))[,1],
             t(apply(NGm, 2, HPDI68))[,2])



dfG = as.data.frame(rbind(KGgkf,NGgkf, KGgkm,NGgkm ))

head(dfG)
names(dfG) <- c("growth", "LC", "UC")
dfG$Pool = factor(c(rep("KG", dim(KGgkf)[1]), rep( "GR", dim(KGgkf)[1]),
                         rep("KG", dim(KGgkf)[1]), rep( "GR", dim(KGgkf)[1])
                         ), levels = c("KG", "GR"))

dfG$Sex = factor(c(rep("F", dim(KGgkf)[1]), rep( "F", dim(KGgkf)[1]),
                    rep("M", dim(KGgkf)[1]), rep( "M", dim(KGgkf)[1])
), levels = c("F", "M"))



dfG$Treatment <- factor(paste(dfG$Sex, dfG$Pool, sep = "-"), levels=c("F-KG", "M-KG", "F-GR", "M-GR" ))
levels(dfG$Treatment)

dfG$z = c(size, size,size, size)

dGpoints = subset(Kdata, surv == 1)

dGpoints$Pool = factor(ifelse(dGpoints$NG==1,  "GR", "KG"), levels = c("KG", "GR"))

dGpoints$Sex <- ifelse(dGpoints$Sex == "M", "M", "F")
dGpoints$Treatment = factor(paste(dGpoints$Sex, dGpoints$Pool, sep = "-"), levels=c("F-KG", "M-KG", "F-GR", "M-GR" ))
levels(dGpoints$Treatment)

dGpoints = dGpoints[,c("SL1_mm", "growth", "Treatment")]


names(dGpoints)[1] = "z" 


#Growth


#GrowthGup_plot
(plotKB = ggplot(dfG, aes(x = z, y = growth, group = Treatment)) + 
    geom_line(size = 0.5,  aes(linetype=Treatment, color=Treatment), show.legend = T) + 
    ylab("Growth ln(final size/initial size)") +
    xlab("Initial size (mm)") +
    geom_point(data = dGpoints, aes(x=z, y=growth, fill= Treatment, shape = Treatment),  show.legend = T, alpha = 0.75, size = 1.75) +
    geom_ribbon( aes(ymin = LC, ymax = UC, fill= Treatment),  show.legend = F, alpha = 0.2, size = 0.01) + 
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = c(0.8,0.75), legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(values=c(2,1,2,1)) +
    scale_shape_manual(values = c(21,24,21,24))+
    scale_color_manual(values=c("black", "black","orange", "orange")) +
    scale_fill_manual(values=c("black","black", "orange", "orange"))  + geom_vline(xintercept=mean(postKG$omega_Intercept), lty=2, size = 0.2) +
    xlim(10,90) #+ 
  #  annotate("text", x = c(25, 50, 75), y=c(106, 106, 106), label = c(TeX("$\\F_{KR>KG}$: 38.9%"), TeX("$\\M_{KR>KG}$: 27.3%"),  TeX("F>M: 80.1%")) )  
)

ppF
ppM
ppFM


summary(fit1)


## Killifish reproduction

KR1 <- brm(Recr ~ z + z2 +GR + z:GR + Density + canopy + (1|stream), family = negbinomial(), subset(Kdata, surv == 1 & Sex != "M"),
           iter = 2000, warmup = 1000, control = list(adapt_delta = 0.90, max_treedepth = 11))

postKR = as.data.frame(as_draws_df(KR1))

names(postKR)


KR_link= function(post,GR, size, center){
  
  z = size - center
  z2 = (size^2) - center^2
  
  if(which(names(post) == "b_z:GR") >2 ){
    mu = post$b_Intercept + post$b_GR *GR + (post$b_z + post$`b_z:GR` * GR) * z + post$b_z2*z2  
  }else{
    mu = post$b_Intercept + post$b_GR *GR + (post$b_z + post$`b_z.GR` * GR) * z + post$b_z2*z2   
    
  }
  
  return(exp(mu))
  
}


size = seq(from=5, to=90, by= 0.5)
KKGr = sapply(1:length(size), function(i)  KR_link(post = postKR,GR = 0, size = size[i], center = center))
KKGr = cbind(apply(KKGr, 2, mean),
             t(apply(KKGr, 2, hdi))[,1],
             t(apply(KKGr, 2, hdi))[,2])

KNGr = sapply(1:length(size), function(i)  KR_link(post = postKR,GR = 1, size = size[i], center = center))
KNGr = cbind(apply(KNGr, 2, mean),
             t(apply(KNGr, 2, hdi))[,1],
             t(apply(KNGr, 2, hdi))[,2])


dfR = as.data.frame(rbind(KKGr,KNGr))
dim(dfR)
names(dfR) <- c("Recr", "LC", "UC")
dfR$Treatment = factor(c(rep("KG", dim(KNGr)[1]), rep( "GR", dim(KNGr)[1])), levels = c("KG", "GR"))
dfR$z = c(size, size)

dRpoints = subset(Kdata, surv == 1 & Recr >= 0 )
dRpoints$Treatment = factor(ifelse(dRpoints$NG==1,  "GR", "KG"), levels = c("KG", "GR"))
dRpoints = dRpoints[,c("SL1_mm", "Recr", "Treatment")]


names(dRpoints)[1] = "z" 

#GrowthGup_plot
(plotKC = ggplot(dRpoints, aes(x = z, y = Recr, colour = Treatment)) + 
    geom_point(alpha = 0.75, size = 1) +
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = "none") +
    #geom_ribbon( aes(ymin = LC, ymax = UC, fill= Treatment), 
    #            alpha = 0.2, size = 0.01, show.legend = F) + 
    ylab("Offspring (N)") +
    xlab("Initial size (mm)") +
    scale_color_manual(values=c("black", "orange")) +
    scale_fill_manual(values=c("black", "orange")) + geom_vline(xintercept=18, lty=2, size = 0.2) +
    xlim(10,90) + ylim(0,60)
)



##


SumTabK <- rbind(my_summary(KS1, variable = "Survival", species = "Killifish"), 
                 my_summary(fit1,variable = "Growth", species = "Killifish"), 
                 my_summary(KR1,variable = "Fecundity", species = "Killifish"))

write.csv(SumTabK, "outputs/Table-S4.csv")

## Recruitment


RecrData <- read.csv("data/Recruitment_data.csv")

head(RecrData) 

RecrData$Pool = factor(paste(RecrData$Removal,RecrData$Sp, sep = "-"))

levels(RecrData$Pool) <- c("KG", "KG", "KR",  "GR")
RecrData$Pool = factor(RecrData$Pool, levels = c("KG", "KR",  "GR"))

# get_prior(Recrt ~ 0 + Pool + (1|Location), family = negbinomial(), Data)


priors2 = prior(normal(0,1), class = b, coef = PoolKR) + # +
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
NKrec = exp(post_Recr$b_Intercept + post_Recr$b_PoolKR)



RecrDF = as.data.frame(cbind(apply(cbind(KGrec,NKrec), 2, mean),
                             t(apply(cbind(KGrec,NKrec), 2, hdi, credMass=.95))[,1],
                             t(apply(cbind(KGrec,NKrec), 2, hdi, credMass=.95))[,2]))

RecrDF$Pool <- factor(c("KG", "KR"), levels = c("KG", "KR"))
names(RecrDF)[1:3] <- c("mean", "LC", "UC" )

RecrDF

pp = LOS(NKrec - KGrec)
pp = paste(round(pp,1),"%", sep = "")






RecrData$Treatment =  RecrData$Pool
levels(RecrData$Treatment) <- c("KG", "KR", "KG",  "GR")
RecrData$Treatment = factor(RecrData$Treatment, levels = c("KG", "KR", "GR"))

(plotD <- ggplot(RecrDF, aes(x=Pool, y=mean, colour = Pool)) + 
    geom_point(size = 4) + scale_color_manual(values=c("black", "orange")) +
    geom_errorbar(aes(ymin=LC, ymax=UC), width=.1) +
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = "none") +
    geom_point(data = subset(RecrData, Sp == "Guppy"), 
               mapping = aes(x = Treatment, y = Recrt, color = Treatment),
               position = position_jitterdodge(dodge.width=1.5),
               size = 1, alpha = 0.7) +
    ylab("Recruits (N)") +
    xlab("Treatment")# + annotate("text", x = 2.5, y=85, label = pp )
  
)



#### Killifish


RectK = subset(RecrData, Sp=="Killifish")

RectK$Treatment = factor(RectK$Pool)
levels(RectK$Treatment)

get_prior(Recrt ~ Pool + Density + Canopy , family = negbinomial(), RectK)

priors3 = prior(normal(0,0.5), class = b, coef = PoolGR) + # +
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
KNGrec = exp(post_RecrK$b_Intercept + post_RecrK$b_PoolGR)


SumTabRec <- rbind(summary(mG)$fixed, summary(mK)$fixed)

(SumTabRec <- round(SumTabRec, 3))

write.csv(SumTabRec, "outputs/Table-Recr.csv")



RecrDF = as.data.frame(cbind(apply(cbind(KKGrec,KNGrec), 2, mean),
                             t(apply(cbind(KKGrec,KNGrec), 2, hdi, credMass=.95))[,1],
                             t(apply(cbind(KKGrec,KNGrec), 2, hdi, credMass=.95))[,2]))

RecrDF$Pool <- factor(c("KG",  "GR"), levels = c("KG",  "GR"))

names(RecrDF)[1:3] <- c("mean", "LC", "UC" )

RecrDF

ppk = LOS(KNGrec - KKGrec)
(ppk = paste(round(ppk,1),"%", sep = ""))

(plotDK <- ggplot(RecrDF, aes(x=Pool, y=mean, colour = Pool)) + 
    geom_point(size = 4) + scale_color_manual(values=c("black", "orange")) +
    geom_errorbar(aes(ymin=LC, ymax=UC), width=.1) +
    theme_classic() + theme(panel.grid.major = element_blank(), legend.position = "none") +
    geom_point(data = RectK, 
               mapping = aes(x = Treatment, y = Recrt, color = Treatment),
               position = position_jitterdodge(dodge.width=1.5),
               size = 1, alpha = 0.7) +
    ylab("Recruits (N)") +
    xlab("Treatment") #+ annotate("text", x = 2.5, y=60, label = ppk )
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




