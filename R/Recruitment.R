library(brms)

Data <- read.csv("data/Recruitment_data.csv")

Data$Pool = factor(paste(Data$Removal,Data$Sp, sep = "-"))

levels(Data$Pool) <- c("KGG", "KGK", "NK", "NG")
Data$Pool = factor(Data$Pool, levels = c("KGG", "NK", "KGK", "NG"))

# get_prior(Recrt ~ 0 + Pool + (1|Location), family = negbinomial(), Data)


prios = prior(normal(3,1), class = b, coef = PoolKGG) +
        prior(normal(3,1), class = b, coef = PoolKGK) +
        prior(normal(3,1), class = b, coef = PoolNG) +
        prior(normal(3,1), class = b, coef = PoolNK) 
mG <- brm( Recrt ~ 0 + Pool + (1|Location), family = negbinomial(), Data, prior = prios,
           iter = 2000, warmup = 1000, control = list(adapt_delta = 0.98, max_treedepth = 13))
post_mG = posterior_samples(mG)
summary(mG)
conditional_effects(mG, effects = "Pool")


library(ggplot2)

Data$Treatment = factor(Data$Removal)
ggplot(Data, aes(x=Pool, y = Recrt, colour = Sp)) + geom_boxplot() 
hist(post_mG$b_PoolNG - post_mG$b_PoolKGK )
hist(post_mG$b_PoolNK - post_mG$b_PoolKGG )


