####  
library(data.table)
library(brms)
rm(list=ls(all=TRUE))



###########################################################################################################
# First get the data
getwd()
setwd("~/Dropbox/Jaime M/Projects_JM/FSU/Pool_manipulation/KG/")
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

S1 <- brm(surv ~ z*NK + z2 +  FishBiom + Density + canopy + (1|stream), family = bernoulli(), Gdata,
         iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))
summary(S1)

conditional_effects(S1, effects = 'z:NK')
Survival_G = Model_selection(S1, name = "Survival", species = "Guppy")
summary(Survival_G)
post = data.frame(posterior_samples(Survival_G))
write.csv(post, "outputs/Post_Survival_G.csv")

# Growth model ------------------------------------------------------------

G1 <- brm(z1 ~ z*NK + + z2 + FishBiom + Density + canopy + (1|stream), family = gaussian(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))


Gdata$growth
Gdata$z2

G2 <- brm(growth ~ z*NK + z2 + FishBiom + Density + canopy + (1|stream), family = gaussian(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

summary(G2)


Growth_G = Model_selection(G2, name = "Growth", "Guppy")
post = data.frame(posterior_samples(Growth_G))
write.csv(post, "outputs/Post_Growth_G.csv")


# Reproduction model ------------------------------------------------------------

R1 <- brm(Recr ~ z*NK  + FishBiom + Density + canopy + (1|stream), family = negbinomial(), subset(Gdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))


Repro_G = Model_selection(R1, name = "Repr", species = "Guppy")
post = data.frame(posterior_samples(Repro_G))
write.csv(post, "outputs/Post_Repr_G.csv")



# ######## Models for killifish -------------------------------------------


# Survival model ----------------------------------------------------------

S1 <- brm(surv ~ z*NG + FishBiom + Density + canopy + (1|stream), family = bernoulli(), Kdata,
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13), chains = 1)
summary(S1)

Survival_K = Model_selection(S1, name = "Survival", species = "Killifish")
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
  b0 ~  NG + FishBiom + Density + canopy +  (1|stream),
  bz1 + bz2 ~ NG,
  omega ~ 1,
  nl = TRUE
)



#make_stancode(bform2, data = subset(Kdata, surv == 1), prior = bprior)
fit1 <- brm(bform, data = subset(Kdata, surv == 1), prior = bprior, chains = 4, save_pars = save_pars(all=TRUE))

# you need the github version of brms for this to run
conditional_effects(fit1, effects = "zr:NG")
summary(fit1)


post = data.frame(posterior_samples(fit1))
write.csv(post, "outputs/Post_Growth_K.csv")


# Reproduction model ------------------------------------------------------------

R1 <- brm(Recr ~ z + FishBiom + Density + canopy + (1|stream), family = negbinomial(), subset(Kdata, surv == 1),
          iter = 2000, warmup = 1000, control = list(adapt_delta = 0.92, max_treedepth = 13))

Repro_K = Model_selection(R1, name = "Reproduction", "Killifish")
post = data.frame(posterior_samples(Repro_K))
write.csv(post, "outputs/Post_Repr_K.csv")

