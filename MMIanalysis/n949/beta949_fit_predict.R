rm(list = ls())
# set path as current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df = read.csv("../mmi_lakecat.csv")
library(dplyr)
library(cobin)
packageVersion("cobin") # > 0.0.1.5
library(spam)
library(spNNGP)
library(coda)
library(mcmcse)
library(loo)
library(bayesplot)
library(betareg)

# step 1: remove one row with Y = 0
df = df[df$MMI_BENT != 0,]
dim(df) # nrow = 949

# step 2: transform covariates
df = df %>% mutate(agkffact_log = log2(1+agkffact),
                   bfi_log = log2(1+bfi),
                   cbnf_log = log2(1+cbnf),
                   conif_log = log2(1+conif),
                   crophay_log = log2(1+crophay),
                   fert_log = log2(1+fert),
                   manure_log = log2(1+manure),
                   pestic1997_log = log2(1+pestic1997),
                   urbmdhi_log = log2(1+urbmdhi)
)

# step 3: standardize covariates
df = df %>% mutate(agkffact_logstd = scale(agkffact_log),
                   bfi_logstd = scale(bfi_log),
                   cbnf_logstd = scale(cbnf_log),
                   conif_logstd = scale(conif_log),
                   crophay_logstd = scale(crophay_log),
                   fert_logstd = scale(fert_log),
                   manure_logstd = scale(manure_log),
                   pestic1997_logstd = scale(pestic1997_log),
                   urbmdhi_logstd = scale(urbmdhi_log)
)

# find mean and sd of each covariate
x_center = df %>% select(agkffact_log,
                         bfi_log,
                         cbnf_log,
                         conif_log,
                         crophay_log,
                         fert_log,
                         manure_log,
                         pestic1997_log,
                         urbmdhi_log) %>% summarise_all(mean)
x_sd = df %>% select(agkffact_log,
                     bfi_log,
                     cbnf_log,
                     conif_log,
                     crophay_log,
                     fert_log,
                     manure_log,
                     pestic1997_log,
                     urbmdhi_log) %>% summarise_all(sd)

scaleback <- function(x, centers, scales){
  as.numeric(c(x[1] - sum((x[-1] * centers) / scales),
               x[-1] / scales))
}

x_center = as.numeric(x_center)
x_sd = as.numeric(x_sd)


M = 15                 # Number of Nearest Neighbors
source("NNmatrix.R")
NN.matrix <- NNMatrix(coords = cbind(df$easting, df$northing), n.neighbors = M, n.omp.threads = 1)
str(NN.matrix)
o =  order(df$easting) # default order for NNGP
all.equal(NN.matrix$ord, o) # should be same as previous order
library(brms)
library(rstan)
nburn = 1000
nsave = 5000
nthin = 1
stan_data = standata(brm(MMI_BENT ~ agkffact_logstd + bfi_logstd +
                           cbnf_logstd + conif_logstd +
                           crophay_logstd + fert_logstd +
                           manure_logstd + pestic1997_logstd +
                           urbmdhi_logstd,
                         data = df[o,], ## important! it requires ordered input
                         family = Beta(), empty = TRUE) )

stan_data$phi_spatial = 0.005
stan_data$M = M
stan_data$NN_ind = NN.matrix$NN_ind
stan_data$NN_dist = NN.matrix$NN_dist
stan_data$NN_distM = NN.matrix$NN_distM

RUN = F
if(RUN){
  fit_beta1 = rstan::stan(file = "betareg_spatial_NNGP.stan",
                          data = stan_data,
                          chains = 1, seed = 3,
                          iter = nburn + nsave, warmup = nburn, thin = nthin)
  saveRDS(fit_beta1, file = "res/fit_beta1_n949.rds")
  fit_beta2 = rstan::stan(file = "betareg_spatial_NNGP.stan",
                          data = stan_data,
                          chains = 1, seed = 33,
                          iter = nburn + nsave, warmup = nburn, thin = nthin)
  saveRDS(fit_beta2, file = "res/fit_beta2_n949.rds")
  fit_beta3 = rstan::stan(file = "betareg_spatial_NNGP.stan",
                          data = stan_data,
                          chains = 1, seed = 333,
                          iter = nburn + nsave, warmup = nburn, thin = nthin)
  saveRDS(fit_beta3, file = "res/fit_beta3_n949.rds")
}else{
  fit_beta1 = readRDS("res/fit_beta1_n949.rds")
  fit_beta2 = readRDS("res/fit_beta2_n949.rds")
  fit_beta3 = readRDS("res/fit_beta3_n949.rds")
}

out1 = rstan::extract(fit_beta1, permute = FALSE)
fit_beta1_betasave = cbind(out1[,1,"b_Intercept"] ,out1[,1,grep("^b\\[\\d+\\]$", dimnames(out1)$parameters)])
fit_beta1_usave = out1[,1,grep("^u\\[\\d+\\]$", dimnames(out1)$parameters)]
fit_beta1_phisave = out1[,1,"phi"]

out2 = rstan::extract(fit_beta2, permute = FALSE)
fit_beta2_betasave = cbind(out2[,1,"b_Intercept"] ,out2[,1,grep("^b\\[\\d+\\]$", dimnames(out2)$parameters)])
fit_beta2_usave = out2[,1,grep("^u\\[\\d+\\]$", dimnames(out2)$parameters)]
fit_beta2_phisave = out2[,1,"phi"]

out3 = rstan::extract(fit_beta3, permute = FALSE)
fit_beta3_betasave = cbind(out3[,1,"b_Intercept"] ,out3[,1,grep("^b\\[\\d+\\]$", dimnames(out3)$parameters)])
fit_beta3_usave = out3[,1,grep("^u\\[\\d+\\]$", dimnames(out3)$parameters)]
fit_beta3_phisave = out3[,1,"phi"]

mycolnames = c("Intercept", "agkffact", "bfi", "cbnf", "conif", "crophay", "fert", "manure", "pestic1997", "urbmdhi")
colnames(fit_beta1_betasave) = mycolnames
colnames(fit_beta2_betasave) = mycolnames
colnames(fit_beta3_betasave) = mycolnames

fit_beta_betasave = as.mcmc.list(list(as.mcmc(fit_beta1_betasave),
                                      as.mcmc(fit_beta2_betasave),
                                      as.mcmc(fit_beta3_betasave)))


X = as.matrix(stan_data$X[order(o),]) # stan NNGP implementation uses permuted data; permute back
# all.equal(X, fit_cobin1$X, check.attributes = FALSE)
# all.equal(X, fit_micobin1$X, check.attributes = FALSE)


fit_beta1_linpredsave = matrix(0, nsave, nrow(df))
fit_beta1_logliksave = matrix(0, nsave, nrow(df))
fit_beta2_linpredsave = matrix(0, nsave, nrow(df))
fit_beta2_logliksave = matrix(0, nsave, nrow(df))
fit_beta3_linpredsave = matrix(0, nsave, nrow(df))
fit_beta3_logliksave = matrix(0, nsave, nrow(df))
for(isave in 1:nsave){
  fit_beta1_linpredsave[isave,] = as.numeric(X %*% fit_beta1_betasave[isave,] + fit_beta1_usave[isave,order(o)]) # permute back random effect
  fit_beta1_logliksave[isave,] = betareg::dbetar(df$MMI_BENT,
                                                 cobin::bftprime(fit_beta1_linpredsave[isave,]),
                                                 fit_beta1_phisave[isave], log = TRUE)
  fit_beta2_linpredsave[isave,] = as.numeric(X %*% fit_beta2_betasave[isave,] + fit_beta2_usave[isave,order(o)]) # permute back random effect
  fit_beta2_logliksave[isave,] = betareg::dbetar(df$MMI_BENT,
                                                 cobin::bftprime(fit_beta2_linpredsave[isave,]),
                                                 fit_beta2_phisave[isave], log = TRUE)
  fit_beta3_linpredsave[isave,] = as.numeric(X %*% fit_beta3_betasave[isave,] + fit_beta3_usave[isave,order(o)]) # permute back random effect
  fit_beta3_logliksave[isave,] = betareg::dbetar(df$MMI_BENT,
                                                 cobin::bftprime(fit_beta3_linpredsave[isave,]),
                                                 fit_beta3_phisave[isave], log = TRUE)
}

mean(c(as.numeric(mcmcse::multiESS(fit_beta1_betasave)),
       as.numeric(mcmcse::multiESS(fit_beta2_betasave)),
       as.numeric(mcmcse::multiESS(fit_beta3_betasave))))

mean(c(sum(get_elapsed_time(fit_beta1)),
       sum(get_elapsed_time(fit_beta2)),
       sum(get_elapsed_time(fit_beta3))))/60

loo::loo(rbind(fit_beta1_logliksave,
               fit_beta2_logliksave,
               fit_beta3_logliksave))

loo::waic(rbind(fit_beta1_logliksave,
                fit_beta2_logliksave,
                fit_beta3_logliksave))

fit_beta_betasave_all = as.matrix(fit_beta_betasave)
fit_beta_betasave_all_orig = coda::mcmc(t(apply(fit_beta_betasave_all, 1, scaleback, centers = x_center,
                                                scales = x_sd)))

mycolnames = c("Intercept", "agkffact", "bfi", "cbnf", "conif", "crophay", "fert", "manure", "pestic1997", "urbmdhi")

colnames(fit_beta_betasave_all_orig) = mycolnames
summary(fit_beta_betasave_all_orig)

fit_beta_parsave <- mcmc.list(
  as.mcmc(cbind(
    as.matrix(fit_beta_betasave[[1]]),
    phi  = out1[,1,"phi"],
    sigmasq = out1[,1,"sigmasq"],
    loglik = rowSums(fit_beta1_logliksave)
  )),
  as.mcmc(cbind(
    as.matrix(fit_beta_betasave[[2]]),
    phi  = out2[,1,"phi"],
    sigmasq = out2[,1,"sigmasq"],
    loglik = rowSums(fit_beta2_logliksave)
  )),
  as.mcmc(cbind(
    as.matrix(fit_beta_betasave[[3]]),
    phi  = out3[,1,"phi"],
    sigmasq = out3[,1,"sigmasq"],
    loglik = rowSums(fit_beta3_logliksave)
  ))
)

coda::gelman.diag(fit_beta_parsave)









# average marginal effect, logit link
nsave = 15000
fit_beta_linpredsave_all = rbind(fit_beta1_linpredsave, fit_beta2_linpredsave, fit_beta3_linpredsave)
avgslope_beta_save = matrix(0, nsave, 9)
colnames(avgslope_beta_save) = mycolnames[-1]
for(isave in 1:nsave){
  for(j in 1:9){
    avgslope_beta_save[isave,j] = mean( cobin::bftprimeprime(fit_beta_linpredsave_all[isave, ]) * fit_beta_betasave_all_orig[isave,j+1])
  }
}
round(colMeans(avgslope_beta_save), 3)
round(apply(avgslope_beta_save, 2, function(x) quantile(x, 0.025)), 3)
round(apply(avgslope_beta_save, 2, function(x) quantile(x, 0.975)), 3)


#fields::quilt.plot(df$lon, df$lat, colMeans(fit_cobin_linpredsave_all))
tt = cbind(round(colMeans(avgslope_beta_save), 3), 
           round(apply(avgslope_beta_save, 2, function(x) quantile(x, 0.025)), 3),
           round(apply(avgslope_beta_save, 2, function(x) quantile(x, 0.975)), 3))

formatted <- apply(tt, 1, function(x) {
  sprintf("%.3f (%.3f, %.3f)", x[1], x[2], x[3])
})

formatted


fit_beta_linpred_mean = colMeans(fit_beta_linpredsave_all)
fit_beta_linpred_sd = apply(fit_beta_linpredsave_all, 2, sd)

fit_beta_mu_mean = colMeans(cobin::bftprime(fit_beta_linpredsave_all))
fit_beta_mu_sd = apply(cobin::bftprime(fit_beta_linpredsave_all), 2, sd)

saveRDS(list(
  coef_std_mean = colMeans(fit_beta_betasave_all),
  coef_std_sd = apply(fit_beta_betasave_all, 2, sd),
  coef_orig_mean = colMeans(fit_beta_betasave_all_orig),
  coef_orig_sd = apply(fit_beta_betasave_all_orig, 2, sd),
  coef_orig_95lo = apply(fit_beta_betasave_all_orig, 2, function(x) quantile(x, 0.025)),
  coef_orig_95up = apply(fit_beta_betasave_all_orig, 2, function(x) quantile(x, 0.975)),
  phi_mean = mean(c(fit_beta1_phisave, fit_beta2_phisave, fit_beta3_phisave)),
  linpred_mean = fit_beta_linpred_mean,
  linpred_sd = fit_beta_linpred_sd,
  mu_mean = fit_beta_mu_mean,
  mu_sd = fit_beta_mu_sd),
  file = "res/summary_beta_n949.rds")

###################################
######### predict #################
###################################


# load prediction data

df4 = read.csv("../lakecat-over40000m2.csv")

################################
# x_center and x_sd must be the one used in this script

df4 = df4 %>% mutate(agkffact_logstd = scale(log2(1+agkffact), center = x_center[1], scale = x_sd[1]),
                     bfi_logstd = scale(log2(1+bfi), center = x_center[2], scale = x_sd[2]),
                     cbnf_logstd = scale(log2(1+cbnf), center = x_center[3], scale = x_sd[3]),
                     conif_logstd = scale(log2(1+conif), center = x_center[4], scale = x_sd[4]),
                     crophay_logstd = scale(log2(1+crophay), center = x_center[5], scale = x_sd[5]),
                     fert_logstd = scale(log2(1+fert), center = x_center[6], scale = x_sd[6]),
                     manure_logstd = scale(log2(1+manure), center = x_center[7], scale = x_sd[7]),
                     pestic1997_logstd = scale(log2(1+pestic1997), center = x_center[8], scale = x_sd[8]),
                     urbmdhi_logstd = scale(log2(1+urbmdhi), center = x_center[9], scale = x_sd[9])
)
# select
Xnew = data.frame(
  df4$agkffact_logstd,
  df4$bfi_logstd,
  df4$cbnf_logstd,
  df4$conif_logstd,
  df4$crophay_logstd,
  df4$fert_logstd,
  df4$manure_logstd,
  df4$pestic1997_logstd,
  df4$urbmdhi_logstd
)
Xnew = as.matrix(cbind(1, Xnew))
str(Xnew)


fit_cobin1 = readRDS("res/fit_micobin1_n949.rds")

#o = order(df$easting) # default order for NNGP

fit_beta1_spNNGPfit = fit_cobin1$spNNGPfit # placeholder

# for(i in M:949){
#   print(all.equal(
#     fit_cobin1$spNNGPfit$neighbor.info$n.indx[[i+1]],
#     stan_data$NN_ind[i,]))
# }

# plug-in posterior samples
out1 = rstan::extract(fit_beta1, permute = FALSE)
fit_beta1_betasave = cbind(out1[,1,"b_Intercept"] ,out1[,1,grep("^b\\[\\d+\\]$", dimnames(out1)$parameters)])
fit_beta1_usave = out1[,1,grep("^u\\[\\d+\\]$", dimnames(out1)$parameters)]


fit_beta1_spNNGPfit$p.beta.samples = fit_beta1_betasave
# order important
fit_beta1_spNNGPfit$p.theta.samples = cbind(out1[,1,"sigmasq"], 0.005)
fit_beta1_spNNGPfit$p.w.samples = t(as.matrix(fit_beta1_usave[,order(o)])) # important to order back


# predict using spNNGP; takes time
pred_beta1 <- predict(fit_beta1_spNNGPfit, 
                      X.0 = Xnew, 
                      coords.0 = cbind(df4$easting,df4$northing),
                      sub.sample = list(start = 1, end = 5000, thin = 1))


transformed_coords4 <- sf::st_coordinates(usmap::usmap_transform(df4[,c("lon","lat")]))

###################################################
# replace NaN to NA: likely due to numerical issue with too close coordinates
pred_beta1$p.w.0[is.nan(pred_beta1$p.w.0)] = NA
pred_beta1$p.y.0[is.nan(pred_beta1$p.y.0)] = NA

wnew_hat_beta = rowMeans(pred_beta1$p.w.0, na.rm = TRUE)
wnew_sd_beta = apply(pred_beta1$p.w.0, 1, sd, na.rm = TRUE)
munew_hat_beta = rowMeans(cobin::bftprime(pred_beta1$p.y.0), na.rm = TRUE)
munew_sd_beta = apply(cobin::bftprime(pred_beta1$p.y.0), 1, sd, na.rm = TRUE)


drawdf4_beta = data.frame(
  x = transformed_coords4[,1],
  y = transformed_coords4[,2],
  urbmdhi = log2(1+df4$urbmdhi),
  wnew_hat = wnew_hat_beta,
  wnew_sd = wnew_sd_beta,
  munew_hat = munew_hat_beta,
  munew_sd = munew_sd_beta
)


saveRDS(drawdf4_beta, "res/predict_beta.rds")
















