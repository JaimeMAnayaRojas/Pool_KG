####  
library(rethinking)
library(data.table)
library(brms)
rm(list=ls(all=TRUE))

# Functions
LOS <- function(x=NULL){
  
  out =   (length(which(x > 0)) / length(x))*100
  
  return(round(out, 3))
  
}


###########################################################################################################
# First get the data
getwd()
#setwd("~/Dropbox/Jaime M/Projects_JM/FSU/Pool_manipulation/KG_git/")


center = 18
# Guppy data cleaning -----------------------------------------------------

Gdata <- read.csv("data/GuppyIPM.csv")
Kdata <- read.csv("data/KillifishIPM.csv")

Gdata$z <- Gdata$SL1_mm - center
Gdata$z1 <- Gdata$SL2_mm 
Gdata <- Gdata[-which(Gdata$Sex2=="M"),]

# Make sure the non reproductive have zeros
Gdata$Recr[which(Gdata$Repr==0 & Gdata$surv ==1)] = 0



# Killifish data cleaning -------------------------------------------------

Kdata$z <- Kdata$SL1_mm - center
Kdata$z1 <- Kdata$SL2_mm 


Kdata <- Kdata[-which(Kdata$Sex2=="M"),]

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

Gdata$growth = log(Gdata$SL2_mm/Gdata$SL1_mm)
Kdata$growth = log(Kdata$SL2_mm/Kdata$SL1_mm)





# Collect the data for the stan model -------------------------------------

data_stan = list(
  
  
  # size of the variables (n) -----------------------------------------------
  
  N_survG = length(Gdata$surv),
  N_growG = length(which(Gdata$surv ==1)),
  N_repG = length(which(Gdata$surv ==1)),
  N_recrG = length(which(Gdata$surv ==1)),
  
  N_survK = length(Kdata$surv),
  N_growK = length(which(Kdata$surv ==1)),
  N_repK = length(which(Kdata$surv ==1)),
  N_recrK = length(which(Kdata$surv ==1)),
  
  N_stream =  length(unique(Gdata$stream)),
  
  
  # Guppy: survival data ----------------------------------------------------
  
  Surv_G = Gdata$surv,
  z_survG = Gdata$z,
  NK_survG = Gdata$NK,
  stream_survG = Gdata$stream,
  canopy_survG = Gdata$canopy,
  BK1_survG = Gdata$K1s,
  BK2_survG = Gdata$K2s,
  BG1_survG = Gdata$G1s,
  BG2_survG = Gdata$G2s,
  FB_survG = Gdata$FishBiom,
  

  # Guppy: growth -----------------------------------------------------------
  z1_G   = Gdata$z1[which(Gdata$surv ==1)],
  z_growG = Gdata$z[which(Gdata$surv ==1)],
  NK_growG = Gdata$NK[which(Gdata$surv ==1)],
  stream_growG = Gdata$stream[which(Gdata$surv ==1)],
  canopy_growG = Gdata$canopy[which(Gdata$surv ==1)],
  BK1_growG = Gdata$K1s[which(Gdata$surv ==1)],
  BK2_growG = Gdata$K2s[which(Gdata$surv ==1)],
  BG1_growG = Gdata$G1s[which(Gdata$surv ==1)],
  BG2_growG = Gdata$G2s[which(Gdata$surv ==1)],
  FB_growG = Gdata$FishBiom[which(Gdata$surv ==1)],
  
  
  # Guppy: Reproduction data ------------------------------------------------
  
  Recr_G = Gdata$Recr[which(Gdata$surv ==1)],
  z_recrG = Gdata$z[which(Gdata$surv ==1)],
  NK_recrG = Gdata$NK[which(Gdata$surv ==1)],
  stream_recrG = Gdata$stream[which(Gdata$surv ==1)],
  canopy_recrG = Gdata$canopy[which(Gdata$surv ==1)],
  BK1_recrG = Gdata$K1s[which(Gdata$surv ==1)],
  BK2_recrG = Gdata$K2s[which(Gdata$surv ==1)],
  BG1_recrG = Gdata$G1s[which(Gdata$surv ==1)],
  BG2_recrG = Gdata$G2s[which(Gdata$surv ==1)],
  FB_recrG = Gdata$FishBiom[which(Gdata$surv ==1)],
  
  # Killifish: survival data ------------------------------------------------
  
  Surv_K = Kdata$surv,
  z_survK = Kdata$z,
  NG_survK = Kdata$NG,
  stream_survK = Kdata$stream,
  canopy_survK = Kdata$canopy,
  BK1_survK = Kdata$K1s,
  BK2_survK = Kdata$K2s,
  BG1_survK = Kdata$G1s,
  BG2_survK = Kdata$G2s,
  FB_survK = Kdata$FishBiom,
  
  
  
  # Killifish: growth data --------------------------------------------------
  
  
  z1_K   = Kdata$z1[which(Kdata$surv ==1)],
  z_growK = Kdata$z[which(Kdata$surv ==1)],
  NG_growK = Kdata$NG[which(Kdata$surv ==1)],
  stream_growK = Kdata$stream[which(Kdata$surv ==1)],
  canopy_growK = Kdata$canopy[which(Kdata$surv ==1)],
  BK1_growK = Kdata$K1s[which(Kdata$surv ==1)],
  BK2_growK = Kdata$K2s[which(Kdata$surv ==1)],
  BG1_growK = Kdata$G1s[which(Kdata$surv ==1)],
  BG2_growK = Kdata$G2s[which(Kdata$surv ==1)],
  FB_growK = Kdata$FishBiom[which(Kdata$surv ==1)],
  
  
  # Killifish: Reproduction data --------------------------------------------
  
  
  Recr_K = Kdata$Recr[which(Kdata$surv ==1)],
  z_recrK = Kdata$z[which(Kdata$surv ==1)],
  NG_recrK = Kdata$NG[which(Kdata$surv ==1)],
  stream_recrK = Kdata$stream[which(Kdata$surv ==1)],
  canopy_recrK = Kdata$canopy[which(Kdata$surv ==1)],
  BK1_recrK = Kdata$K1s[which(Kdata$surv ==1)],
  BK2_recrK = Kdata$K2s[which(Kdata$surv ==1)],
  BG1_recrK = Kdata$G1s[which(Kdata$surv ==1)],
  BG2_recrK = Kdata$G2s[which(Kdata$surv ==1)],
  FB_recrK = Kdata$FishBiom[which(Kdata$surv ==1)]
  

  
)
# 

modG = stan("R/models/pool_mod.stan", data = data_stan, cores = 4, chains = 4, 
            iter = 6000, warmup = 4500, control = list(adapt_delta = 0.92, max_treedepth = 12))
saveRDS(modG, "outputs/Model_KG.RDS")


modG = readRDS("outputs/Model_KG.RDS")
sum = as.data.frame(precis(modG, digits = 5, prob = .95, depth = 2))
sum$Pars = rownames(sum)


precis(modG, digits = 5, prob = .95, depth = 2)




#traceplot_ulam(modG)

#

post <- as.data.frame(extract(modG))
write.csv(as.data.frame(post), "outputs/Posteriors.csv")

PP = as.vector(apply(as.data.frame(post), 2, LOS))

PP

names(post)
model.summary <- as.data.frame(precis(modG, prob = .95, digits = 3, depth = 3))

model.summary$LOS_l = PP[-which(names(post) == "lp__")]
model.summary$LOS_u = 100 - model.summary$LOS_l
model.summary[model.summary$`2.5%` > 0,]
model.summary[model.summary$`97.5%` < 0,] 


model.summary[model.summary$LOS_l > 95,] 
model.summary[model.summary$LOS_l < 5,] 

write.csv(model.summary, "outputs/Model_sum.csv")
