####  
library(data.table)
library(brms)
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

S1 <- brm(surv ~ z + z2 + NK + z:NK + z2:NK  + Density + canopy + (1|stream), family = bernoulli(), Gdata,
         iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

S2 <- brm(surv ~ z + z2 + NK + z:NK  + Density + canopy + (1|stream), family = bernoulli(), Gdata,
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

S3 <- brm(surv ~ z + NK + z:NK  + Density + canopy + (1|stream), family = bernoulli(), Gdata,
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

S4 <- brm(surv ~ z + NK + Density + canopy + (1|stream), family = bernoulli(), Gdata,
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

S5 <- brm(surv ~ NK + Density + canopy + (1|stream), family = bernoulli(), Gdata,
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))


S1 <- add_criterion(S1, criterion = c("loo", "waic"))
S2 <- add_criterion(S2, criterion = c("loo", "waic"))
S3 <- add_criterion(S3, criterion = c("loo", "waic"))
S4 <- add_criterion(S4, criterion = c("loo", "waic"))
S5 <- add_criterion(S5, criterion = c("loo", "waic"))


as.data.frame(loo_compare(S1,S2,S3,S4,S5, criterion = "loo"))
loo_model_weights(S1,S2,S3,S4,S5, criterion = "loo")

summary(S5)
conditional_effects(S5, effects = 'NK')
Survival_G = Model_selection(S5, name = "Survival", species = "Guppy")
summary(Survival_G)
post = data.frame(posterior_samples(Survival_G))

apply(post, 2, LOS)
write.csv(post, "outputs/Post_Survival_G.csv")

names(post)
# Growth model ------------------------------------------------------------
G1 <- brm(growth ~ z + z2 + NK + z:NK + z2:NK  + Density + canopy + (1|stream), family = gaussian(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

G2 <- brm(growth ~ z + z2 + NK + z:NK  + Density + canopy + (1|stream), family = gaussian(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

G3 <- brm(growth ~ z  + NK + z:NK   + Density + canopy + (1|stream), family = gaussian(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))


G1 <- add_criterion(G1, criterion = c("loo", "waic"))
G2 <- add_criterion(G2, criterion = c("loo", "waic"))
G3 <- add_criterion(G3, criterion = c("loo", "waic"))

as.data.frame(loo_compare(G1,G2,G3, criterion = "loo"))
loo_model_weights(G1,G2,G3, criterion = "loo")

summary(G2)

Growth_G = Model_selection(G2, name = "Growth", "Guppy")
summary(Growth_G)
post = data.frame(posterior_samples(Growth_G))

apply(post, 2, LOS)

write.csv(post, "outputs/Post_Growth_G.csv")


# Reproduction model ------------------------------------------------------------

R1 <- brm(Recr ~ z + z2 + NK + z:NK + z2:NK   + Density + canopy + (1|stream), family = negbinomial(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

R2 <- brm(Recr ~ z + z2 + NK + z:NK   + Density + canopy + (1|stream), family = negbinomial(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

R3 <- brm(Recr ~ z + NK + z:NK  + Density + canopy + (1|stream), family = negbinomial(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))


R1 <- add_criterion(R1, criterion = c("loo", "waic"))
R2 <- add_criterion(R2, criterion = c("loo", "waic"))
R3 <- add_criterion(R3, criterion = c("loo", "waic"))

as.data.frame(loo_compare(R1,R2,R3, criterion = "loo"))
loo_model_weights(R1,R2,R3, criterion = "loo")


Repro_G = Model_selection(R2, name = "Repr", species = "Guppy")
Repro_G
post = data.frame(posterior_samples(Repro_G))
apply(post, 2, LOS)


write.csv(post, "outputs/Post_Repr_G.csv")



# Model structure for guppies


Mod_Sec = rbind(as.data.frame(loo_compare(S1,S2,S3, criterion = "loo")), 
as.data.frame(loo_compare(G1,G2,G3, criterion = "loo")),
as.data.frame(loo_compare(R1,R2,R3, criterion = "loo")))

write.csv(Mod_Sec, "outputs/Model_structuresG.csv")

# ######## Models for killifish -------------------------------------------


# Survival model ----------------------------------------------------------

S1 <- brm(surv ~ z + z2 + NG + z:NG + z2:NG  + Density + canopy + (1|stream), family = bernoulli(), Kdata,
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13), chains = 1)

S2 <- brm(surv ~ z + z2 + NG + z:NG   + Density + canopy + (1|stream), family = bernoulli(), Kdata,
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13), chains = 1)

S3 <- brm(surv ~ z +  NG + z:NG   + Density + canopy + (1|stream), family = bernoulli(), Kdata,
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13), chains = 1)



S1 <- add_criterion(S1, criterion = c("loo", "waic"))
S2 <- add_criterion(S2, criterion = c("loo", "waic"))
S3 <- add_criterion(S3, criterion = c("loo", "waic"))

as.data.frame(loo_compare(S1,S2,S3, criterion = "loo"))
loo_model_weights(S1,S2,S3, criterion = "loo")



Survival_K = Model_selection(S2, name = "Survival", species = "Killifish")
summary(Survival_K)
post = data.frame(posterior_samples(Survival_K))
write.csv(post, "outputs/Post_Survival_K.csv")

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


bform2 <- bf(
  growth ~ b0 + bz1 * (zr - omega) * step(omega - zr) + 
    bz2 * (zr - omega) * step(zr - omega),
  b0 ~  NG + canopy +  (1|stream),
  bz1 + bz2 ~ NG,
  omega ~ 1,
  nl = TRUE
)

#make_stancode(bform2, data = subset(Kdata, surv == 1), prior = bprior)
fit2 <- brm(bform2, data = subset(Kdata, surv == 1), prior = bprior, chains = 4, save_pars = save_pars(all=TRUE))



bform3 <- bf(
  growth ~ b0 + bz1 * (zr - omega) * step(omega - zr) + 
    bz2 * (zr - omega) * step(zr - omega),
  b0 ~  NG   + Density +  (1|stream),
  bz1 + bz2 ~ NG,
  omega ~ 1,
  nl = TRUE
)

#make_stancode(bform2, data = subset(Kdata, surv == 1), prior = bprior)
fit3 <- brm(bform3, data = subset(Kdata, surv == 1), prior = bprior, chains = 4, save_pars = save_pars(all=TRUE))


bform4 <- bf(
  growth ~ b0 + bz1 * (zr - omega) * step(omega - zr) + 
    bz2 * (zr - omega) * step(zr - omega),
  b0 ~  NG    +  (1|stream),
  bz1 + bz2 ~ NG,
  omega ~ 1,
  nl = TRUE
)

#make_stancode(bform2, data = subset(Kdata, surv == 1), prior = bprior)
G4 <- brm(bform6, data = subset(Kdata, surv == 1), prior = bprior, chains = 4, save_pars = save_pars(all=TRUE))



G1 <- add_criterion(fit1, criterion = c("loo", "waic"))
G2 <- add_criterion(fit2, criterion = c("loo", "waic"))
G3 <- add_criterion(fit3, criterion = c("loo", "waic"))
G4 <- add_criterion(fit4, criterion = c("loo", "waic"))

as.data.frame(loo_compare(G1,G2, G3, G4,  criterion = "loo"))
loo_model_weights(G1,G2, G3, G4, criterion = "loo")

write.csv(as.data.frame(loo_compare(G1,G2, G3, G4,  criterion = "loo")), "outputs/Growth_Killifish.csv")
# you need the github version of brms for this to run
conditional_effects(G3, effects = "zr:NG")
summary(G3)


post = data.frame(posterior_samples(G3))

apply(post, 2, mean)
write.csv(post, "outputs/Post_Growth_K.csv")


# Reproduction model ------------------------------------------------------------

R1 <- brm(Recr ~ z + NG  + Density + canopy + (1|stream), family = negbinomial(), subset(Kdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

Repro_K = Model_selection(R1, name = "Reproduction", "Killifish")

post = data.frame(posterior_samples(Repro_K))
Repro_K

apply(post, 2, LOS)

write.csv(post, "outputs/Post_Repr_K.csv")

