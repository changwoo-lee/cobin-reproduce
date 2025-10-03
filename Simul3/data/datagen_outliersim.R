rm(list = ls())
library(cobin)
library(betareg)
nsim = 110 # 10 more in case for failed simul 

n = 500
dataall = list()
for(isim in 1:nsim){
  df_outlier_y10m2 = data.frame(y = 0.01, x1 = 2, x2 = 2)
  df_outlier_y10m3 = data.frame(y = 0.001, x1 = 2, x2 = 2)
  
  set.seed(isim*3333)
  
  x1 = rnorm(n); x2 = rnorm(n);
  # beta
  y = betareg::rbetar(n, cobin::bftprime(-5 + x1), phi = rep(17,n))
  df_beta = data.frame(y = y, x1 = x1, x2 = x2)
  
  # cobin
  y = cobin::rcobin(n, -5 + x1, rep(6,n))
  df_cobin = data.frame(y = y, x1 = x1, x2 = x2)
  # first observation is outlier
  dataall[[isim]] = list("beta_y10m2" = rbind(df_outlier_y10m2, df_beta), 
                         "beta_y10m3" = rbind(df_outlier_y10m3, df_beta),
                         "cobin_y10m2" = rbind(df_outlier_y10m2, df_cobin),  
                         "cobin_y10m3" = rbind(df_outlier_y10m3, df_cobin))
}


saveRDS(dataall, file = paste0("data_outliersim.rds"))


