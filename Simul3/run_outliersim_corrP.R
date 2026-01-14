rm(list = ls())

# slurm_id corresponds to 1, 2, ..., 100 for each simulation
#
# slurm <- TRUE
# if(slurm){
#   slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# } else{
#   slurm_id <- 0
# 
# }

# manual setting of slurm_id for this example (actual sim was run with slurm job array)
slurm_id = 1 


# load data
dataall = readRDS(file = "data/data_outliersim_corrP.rds")
data = dataall[[slurm_id]]
library(rstan)
library(brms)
library(cobin)

nsave = 5000
nburn = 1000
nthin = 1


# beta data, y=10^-2

df_results = data.frame()

for(jj in 1:4){
  ## jj = 1 : beta_y10m2
  ## jj = 2 : beta_y10m3
  ## jj = 3 : cobin_y10m2
  ## jj = 4 : cobin_y10m3
  df_outlier = data[[jj]]
  df = df_outlier[-1,]
  
  ######## beta regression ##########
  # whole dataset
  stan_data_o = brms::standata(brm(y ~ ., data = df_outlier, 
                                   family = Beta(), empty = TRUE) )
  fit_beta_o = rstan::stan(file = "betareg_fixedeffect.stan",
                           data = stan_data_o,
                           chains = 1, seed = jj*slurm_id*4,
                           iter = nburn + nsave, warmup = nburn, thin = nthin)  
  
  fit_beta_o_samp = extract(fit_beta_o, permute = FALSE)
  beta_mess = mcmcse::multiESS(fit_beta_o_samp[,1,c("b_Intercept","b[1]","b[2]")])
  runtime = sum(rstan::get_elapsed_time(fit_beta_o))
  
  # including intercept
  betahat_o = summary(fit_beta_o)$summary[c("b_Intercept","b[1]","b[2]"),"mean"]
  beta_o_lbd = summary(fit_beta_o)$summary[c("b_Intercept","b[1]","b[2]"),"2.5%"]
  beta_o_ubd = summary(fit_beta_o)$summary[c("b_Intercept","b[1]","b[2]"),"97.5%"]
  
  res_o = cbind(matrix(betahat_o, ncol = 3), 
                matrix(beta_o_lbd, ncol = 3), matrix(beta_o_ubd, ncol = 3))
  
  colnames(res_o) = c("betahat_o_x0", "betahat_o_x1", "betahat_o_x2",
                      "beta_o_lbd_x0", "beta_o_lbd_x1", "beta_o_lbd_x2",
                      "beta_o_ubd_x0", "beta_o_ubd_x1", "beta_o_ubd_x2")
  # linear predictor for non-outliers
  eta_o = cbind(1, df$x1,df$x2)%*%betahat_o
  
  
  # without outlier
  stan_data = brms::standata(brm(y ~ ., data = df, 
                                 family = Beta(), empty = TRUE) )
  fit_beta = rstan::stan("betareg_fixedeffect.stan",
                               data = stan_data,
                               chains = 1, seed = jj*slurm_id*4,
                               iter = nburn + nsave, warmup = nburn, thin = nthin)

  betahat = summary(fit_beta)$summary[c("b_Intercept","b[1]","b[2]"),"mean"]
  beta_lbd = summary(fit_beta)$summary[c("b_Intercept","b[1]","b[2]"),"2.5%"]
  beta_ubd = summary(fit_beta)$summary[c("b_Intercept","b[1]","b[2]"),"97.5%"]
  
  res = cbind(matrix(betahat, ncol = 3), 
                matrix(beta_lbd, ncol = 3), matrix(beta_ubd, ncol = 3))
  
  colnames(res) = c("betahat_x0", "betahat_x1", "betahat_x2",
                      "beta_lbd_x0", "beta_lbd_x1", "beta_lbd_x2",
                      "beta_ubd_x0", "beta_ubd_x1", "beta_ubd_x2")
  
  # linear predictor for non-outliers
  eta = cbind(1, df$x1,df$x2)%*%betahat
  
  res_diff = res_o[1,1:3,drop=F] - res[1,1:3,drop=F]
  colnames(res_diff) = c("d_betahat_x0", "d_betahat_x1", "d_betahat_x2")
  
  df_temp = data.frame("dataname" = names(data)[[jj]],
                       "model" = "betareg",
                       "isim" = slurm_id,
                       "d_etahat" = sqrt(sum((eta_o - eta)^2)), 
                       "beta_mess" = beta_mess,
                       "runtime" = runtime)
  
  df_temp = cbind(df_temp, res_diff,res_o, res)
  df_results = rbind(df_results, df_temp)
  
  ######## beta rect regression ##########
  # whole dataset
  stan_data_o = brms::standata(brm(y ~ ., data = df_outlier, 
                                   family = Beta(), empty = TRUE) )
  fit_betarec_o = rstan::stan("betarecreg_fixedeffect.stan",
                               data = stan_data_o,
                               chains = 1, seed = jj*slurm_id*4,
                               iter = nburn + nsave, warmup = nburn, thin = nthin)

  fit_betarec_o_samp = extract(fit_betarec_o, permute = FALSE)
  beta_mess = mcmcse::multiESS(fit_betarec_o_samp[,1,c("b_Intercept","b[1]","b[2]")])
  runtime = sum(rstan::get_elapsed_time(fit_betarec_o))
  
  # including intercept
  betahat_o = summary(fit_betarec_o)$summary[c("b_Intercept","b[1]","b[2]"),"mean"]
  beta_o_lbd = summary(fit_betarec_o)$summary[c("b_Intercept","b[1]","b[2]"),"2.5%"]
  beta_o_ubd = summary(fit_betarec_o)$summary[c("b_Intercept","b[1]","b[2]"),"97.5%"]
  
  res_o = cbind(matrix(betahat_o, ncol = 3), 
                matrix(beta_o_lbd, ncol = 3), matrix(beta_o_ubd, ncol = 3))
  
  colnames(res_o) = c("betahat_o_x0", "betahat_o_x1", "betahat_o_x2",
                      "beta_o_lbd_x0", "beta_o_lbd_x1", "beta_o_lbd_x2",
                      "beta_o_ubd_x0", "beta_o_ubd_x1", "beta_o_ubd_x2")
  # linear predictor for non-outliers
  eta_o = cbind(1, df$x1,df$x2)%*%betahat_o
  
  
  # without outlier
  stan_data = brms::standata(brm(y ~ ., data = df, 
                                 family = Beta(), empty = TRUE) )
  fit_betarec = rstan::stan("betarecreg_fixedeffect.stan",
                               data = stan_data,
                               chains = 1, seed = jj*slurm_id*4,
                               iter = nburn + nsave, warmup = nburn, thin = nthin)

  betahat = summary(fit_betarec)$summary[c("b_Intercept","b[1]","b[2]"),"mean"]
  beta_lbd = summary(fit_betarec)$summary[c("b_Intercept","b[1]","b[2]"),"2.5%"]
  beta_ubd = summary(fit_betarec)$summary[c("b_Intercept","b[1]","b[2]"),"97.5%"]
  
  res = cbind(matrix(betahat, ncol = 3), 
              matrix(beta_lbd, ncol = 3), matrix(beta_ubd, ncol = 3))
  
  colnames(res) = c("betahat_x0", "betahat_x1", "betahat_x2",
                    "beta_lbd_x0", "beta_lbd_x1", "beta_lbd_x2",
                    "beta_ubd_x0", "beta_ubd_x1", "beta_ubd_x2")
  
  # linear predictor for non-outliers
  eta = cbind(1, df$x1,df$x2)%*%betahat
  
  res_diff = res_o[1,1:3,drop=F] - res[1,1:3,drop=F]
  colnames(res_diff) = c("d_betahat_x0", "d_betahat_x1", "d_betahat_x2")
  
  df_temp = data.frame("dataname" = names(data)[[jj]],
                       "model" = "betarecreg",
                       "isim" = slurm_id,
                       "d_etahat" = sqrt(sum((eta_o - eta)^2)), 
                       "beta_mess" = beta_mess,
                       "runtime" = runtime)
  
  df_temp = cbind(df_temp, res_diff,res_o, res)
  df_results = rbind(df_results, df_temp)
  
  
  ######## cobin regression ##########
  set.seed(jj*slurm_id*5)
  fit_cobin_o = cobin::cobinreg(y ~ ., data = df_outlier, link = "cobit",
                                nburn = nburn, nsave = nsave, nthin = nthin)
  
  beta_mess = mcmcse::multiESS(fit_cobin_o$post_save[,1:3])
  runtime = as.numeric(fit_cobin_o$t_mcmc)
  
  betahat_o = colMeans(fit_cobin_o$post_save[,1:3])
  beta_o_lbd = apply(fit_cobin_o$post_save[,1:3], 2, function(x) quantile(x, probs = 0.025))
  beta_o_ubd = apply(fit_cobin_o$post_save[,1:3], 2, function(x) quantile(x, probs = 0.975))
  
  
  res_o = cbind(matrix(betahat_o, ncol = 3), 
                matrix(beta_o_lbd, ncol = 3), matrix(beta_o_ubd, ncol = 3))
  # linear predictor for non-outliers
  eta_o = cbind(1, df$x1,df$x2)%*%betahat_o
  
  colnames(res_o) = c("betahat_o_x0", "betahat_o_x1", "betahat_o_x2",
                      "beta_o_lbd_x0", "beta_o_lbd_x1", "beta_o_lbd_x2",
                      "beta_o_ubd_x0", "beta_o_ubd_x1", "beta_o_ubd_x2")
  
  set.seed(jj*slurm_id*55)
  fit_cobin = cobin::cobinreg(y ~ ., data = df, link = "cobit",
                              nburn = nburn, nsave = nsave, nthin = nthin)
  betahat = colMeans(fit_cobin$post_save[,1:3])
  beta_lbd = apply(fit_cobin$post_save[,1:3], 2, function(x) quantile(x, probs = 0.025))
  beta_ubd = apply(fit_cobin$post_save[,1:3], 2, function(x) quantile(x, probs = 0.975))
  
  res = cbind(matrix(betahat, ncol = 3), 
              matrix(beta_lbd, ncol = 3), matrix(beta_ubd, ncol = 3))
  # linear predictor for non-outliers
  eta = cbind(1, df$x1,df$x2)%*%betahat
  
  colnames(res) = c("betahat_x0", "betahat_x1", "betahat_x2",
                    "beta_lbd_x0", "beta_lbd_x1", "beta_lbd_x2",
                    "beta_ubd_x0", "beta_ubd_x1", "beta_ubd_x2")
  
  res_diff = res_o[1,1:3,drop=F] - res[1,1:3,drop=F]
  colnames(res_diff) = c("d_betahat_x0", "d_betahat_x1", "d_betahat_x2")
  
  df_temp = data.frame("dataname" = names(data)[[jj]],
                       "model" = "cobinreg",
                       "isim" = slurm_id,
                       "d_etahat" = sqrt(sum((eta_o - eta)^2)),
                       "beta_mess" = beta_mess,
                       "runtime" = runtime)
  df_temp = cbind(df_temp, res_diff,res_o, res)
  df_results = rbind(df_results, df_temp)
  
  ######## micobin regression ##########
  set.seed(jj*slurm_id*555)
  fit_micobin_o = cobin::micobinreg(y ~ ., data = df_outlier, link = "cobit",
                                    nburn = nburn, nsave = nsave, nthin = nthin)
  beta_mess = mcmcse::multiESS(fit_micobin_o$post_save[,1:3])
  runtime = as.numeric(fit_micobin_o$t_mcmc)
  
  betahat_o = colMeans(fit_micobin_o$post_save[,1:3])
  beta_o_lbd = apply(fit_micobin_o$post_save[,1:3], 2, function(x) quantile(x, probs = 0.025))
  beta_o_ubd = apply(fit_micobin_o$post_save[,1:3], 2, function(x) quantile(x, probs = 0.975))
  res_o = cbind(matrix(betahat_o, ncol = 3), 
                matrix(beta_o_lbd, ncol = 3), matrix(beta_o_ubd, ncol = 3))
  
  colnames(res_o) = c("betahat_o_x0", "betahat_o_x1", "betahat_o_x2",
                      "beta_o_lbd_x0", "beta_o_lbd_x1", "beta_o_lbd_x2",
                      "beta_o_ubd_x0", "beta_o_ubd_x1", "beta_o_ubd_x2")
  # linear predictor for non-outliers
  eta_o = cbind(1, df$x1,df$x2)%*%betahat_o
  
  set.seed(jj*slurm_id*5555)
  fit_micobin = cobin::micobinreg(y ~ ., data = df, link = "cobit",
                                  nburn = nburn, nsave = nsave, nthin = nthin)
  betahat = colMeans(fit_micobin$post_save[,1:3])
  beta_lbd = apply(fit_micobin$post_save[,1:3], 2, function(x) quantile(x, probs = 0.025))
  beta_ubd = apply(fit_micobin$post_save[,1:3], 2, function(x) quantile(x, probs = 0.975))
  res = cbind(matrix(betahat, ncol = 3), 
              matrix(beta_lbd, ncol = 3), matrix(beta_ubd, ncol = 3))
  
  colnames(res) = c("betahat_x0", "betahat_x1", "betahat_x2",
                    "beta_lbd_x0", "beta_lbd_x1", "beta_lbd_x2",
                    "beta_ubd_x0", "beta_ubd_x1", "beta_ubd_x2")
  # linear predictor for non-outliers
  eta = cbind(1, df$x1,df$x2)%*%betahat
  
  res_diff = res_o[1,1:3,drop=F] - res[1,1:3,drop=F]
  colnames(res_diff) = c("d_betahat_x0", "d_betahat_x1", "d_betahat_x2")
  
  df_temp = data.frame("dataname" = names(data)[[jj]],
                       "model" = "micobinreg",
                       "isim" = slurm_id,
                       "d_etahat" = sqrt(sum((eta_o - eta)^2)),
                       "beta_mess" = beta_mess,
                       "runtime" = runtime)
  df_temp = cbind(df_temp, res_diff,res_o, res)
  
  df_results = rbind(df_results, df_temp)
}
saveRDS(df_results, file = paste0("out_outliersim_corrP/outliersim_", slurm_id, ".rds"))

devtools::session_info()



