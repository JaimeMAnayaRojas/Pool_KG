Model_selection <- function(Full_model, name, species){
  M1 = Full_model
  M1 <- add_criterion(M1, c("loo", "waic"))
  
  ## remove effects
  M2 <- update(M1, formula. = ~ . - FishBiom)
  M2 <- add_criterion(M2, c("loo", "waic"))
  
  
  M3 <- update(M1, formula. = ~ . - Density)
  M3 <- add_criterion(M3, c("loo", "waic"))
  
  
  M4 <- update(M1, formula. = ~ . - canopy)
  M4 <- add_criterion(M4, c("loo", "waic"))
  
  M5 <- update(M1, formula. = ~ . -FishBiom - canopy)
  M5 <- add_criterion(M5, c("loo", "waic"))
  
  
  M6 <- update(M1, formula. = ~ . -Density - canopy)
  M6 <- add_criterion(M6, c("loo", "waic"))
  
  
  M7 <- update(M1, formula. = ~ . -Density - FishBiom - canopy)
  M7 <- add_criterion(M7, c("loo", "waic"))
  
  tab = data.frame(loo_compare(M1,M2, M3, M4, M5,M6, M7, criterion = "loo"))
  mod_lis = list(M1 = M1, M2= M2, M3= M3, M4= M4, M5= M5, M6= M6, M7= M7)

  
  
  
  tab$Formula = c(mod_lis[rownames(tab)[1]]$formula[1]$formula[3],
                  mod_lis[rownames(tab)[2]]$formula[1]$formula[3],
                  mod_lis[rownames(tab)[3]]$formula[1]$formula[3],
                  mod_lis[rownames(tab)[4]]$formula[1]$formula[3],
                  mod_lis[rownames(tab)[5]]$formula[1]$formula[3],
                  mod_lis[rownames(tab)[6]]$formula[1]$formula[3],
                  mod_lis[rownames(tab)[7]]$formula[1]$formula[3]
  )
  
  write.csv(tab, paste("outputs/",name, "_" ,species,".csv", sep = ""))
  

  return(get(rownames(tab)[1]))
  
  post= posterior_samples(get(rownames(tab)[1]))
  
  write.csv(post, paste("post_", name, species, ".csv", sep = ""))
}

LOS = function(x){
  
  p =  100 * length(which(x > 0))/length(x)
  
  return(p)
  
}
