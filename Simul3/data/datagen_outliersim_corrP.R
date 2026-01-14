rm(list = ls())
library(cobin)
library(betareg)
nsim = 110 # 10 more in case for failed simul
library(MASS)

n = 500
dataall = list()
for(isim in 1:nsim){
  df_outlier_y10m2 = data.frame(y = 0.01, x1 = 6, x2 = 6)
  df_outlier_y10m3 = data.frame(y = 0.001, x1 = 6, x2 = 6)

  set.seed(isim*33333)

  #x1 = rnorm(n); x2 = rnorm(n);
  # correlated with 0.9
  Sigma = 9*matrix(c(1,0.9,0.9,1),2,2)
  X = MASS::mvrnorm(n, mu = c(0,0), Sigma = Sigma)
  x1 = X[,1]; x2 = X[,2]
  # beta
  y = betareg::rbetar(n, cobin::bftprime(-6 + x1), phi = rep(17,n))
  df_beta = data.frame(y = y, x1 = x1, x2 = x2)

  # cobin
  y = cobin::rcobin(n, -6 + x1, rep(6,n))
  df_cobin = data.frame(y = y, x1 = x1, x2 = x2)
  # first observation is outlier
  dataall[[isim]] = list("beta_y10m2" = rbind(df_outlier_y10m2, df_beta),
                         "beta_y10m3" = rbind(df_outlier_y10m3, df_beta),
                         "cobin_y10m2" = rbind(df_outlier_y10m2, df_cobin),
                         "cobin_y10m3" = rbind(df_outlier_y10m3, df_cobin))
}


saveRDS(dataall, file = paste0("data_outliersim_corrP.rds"))


