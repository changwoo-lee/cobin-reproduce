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
# also remove two rows with comid = 9201925, 22845861
excludeidx = which(df$comid %in% c(9201925, 22845861))
df[excludeidx,]
df = df[-excludeidx,]

dim(df) # nrow = 947

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


# step 3 : fit the model
RUN = TRUE


o = order(df$easting) # default order for NNGP


# step 3.1. cobin

if(RUN){
  set.seed(1)
  fit_cobin1 = spcobinreg(MMI_BENT ~ agkffact_logstd + bfi_logstd +
                            cbnf_logstd + conif_logstd +
                            crophay_logstd + fert_logstd +
                            manure_logstd + pestic1997_logstd +
                            urbmdhi_logstd,
                          data = df, coords = cbind(df$easting, df$northing),
                          NNGP = T,
                          priors = list(beta_intercept_scale = 100, beta_scale = 100, beta_df = Inf,
                                        phi_lb = 0.005, phi_ub = 0.005),
                          nngp.control = list(n.neighbors = 15, ord = o),
                          nburn = 1000, nsave = 5000)
  #saveRDS(fit_cobin1, file = "fit_cobin1_n949.rds")
  set.seed(11)
  fit_cobin2 = spcobinreg(MMI_BENT ~ agkffact_logstd + bfi_logstd +
                            cbnf_logstd + conif_logstd +
                            crophay_logstd + fert_logstd +
                            manure_logstd + pestic1997_logstd +
                            urbmdhi_logstd,
                          data = df, coords = cbind(df$easting, df$northing),
                          NNGP = T,
                          priors = list(beta_intercept_scale = 100, beta_scale = 100, beta_df = Inf,
                                        phi_lb = 0.005, phi_ub = 0.005),
                          nngp.control = list(n.neighbors = 15, ord = o),
                          nburn = 1000, nsave = 5000)
  #saveRDS(fit_cobin2, file = "fit_cobin2_n950.rds")
  set.seed(111)
  fit_cobin3 = spcobinreg(MMI_BENT ~ agkffact_logstd + bfi_logstd +
                            cbnf_logstd + conif_logstd +
                            crophay_logstd + fert_logstd +
                            manure_logstd + pestic1997_logstd +
                            urbmdhi_logstd,
                          data = df, coords = cbind(df$easting, df$northing),
                          NNGP = T,
                          priors = list(beta_intercept_scale = 100, beta_scale = 100, beta_df = Inf,
                                        phi_lb = 0.005, phi_ub = 0.005),
                          nngp.control = list(n.neighbors = 15, ord = o),
                          nburn = 1000, nsave = 5000)
  #saveRDS(fit_cobin3, file = "fit_cobin3_n950.rds")
}else{
  fit_cobin1 = readRDS("fit_cobin1_n947.rds")
  fit_cobin2 = readRDS("fit_cobin2_n947.rds")
  fit_cobin3 = readRDS("fit_cobin3_n947.rds")
}


fit_cobin_betasave = as.mcmc.list(list(fit_cobin1$post_save[,1:ncol(fit_cobin1$X)],
                                       fit_cobin2$post_save[,1:ncol(fit_cobin2$X)],
                                       fit_cobin3$post_save[,1:ncol(fit_cobin3$X)]))
coda::gelman.diag(fit_cobin_betasave)
bayesplot::mcmc_trace(fit_cobin_betasave)# based on standardized covariates

mean(c(as.numeric(fit_cobin1$t_mcmc),
       as.numeric(fit_cobin2$t_mcmc),
       as.numeric(fit_cobin3$t_mcmc)))/60

mean(c(as.numeric(mcmcse::multiESS(fit_cobin_betasave[[1]])),
       as.numeric(mcmcse::multiESS(fit_cobin_betasave[[2]])),
       as.numeric(mcmcse::multiESS(fit_cobin_betasave[[3]]))))

loo::waic(rbind(fit_cobin1$loglik_save,
                fit_cobin2$loglik_save,
                fit_cobin3$loglik_save))

fit_cobin_betasave_all = as.matrix(fit_cobin_betasave)
fit_cobin_betasave_all_orig = coda::mcmc(t(apply(fit_cobin_betasave_all, 1, scaleback, centers = x_center,
                                                 scales = x_sd)))

mycolnames = c("Intercept", "agkffact", "bfi", "cbnf", "conif", "crophay", "fert", "manure", "pestic1997", "urbmdhi")

colnames(fit_cobin_betasave_all_orig) = mycolnames
summary(fit_cobin_betasave_all_orig)


# step 3.2. micobin

if(RUN){
  set.seed(2)
  fit_micobin1 = spmicobinreg(MMI_BENT ~ agkffact_logstd + bfi_logstd +
                                cbnf_logstd + conif_logstd +
                                crophay_logstd + fert_logstd +
                                manure_logstd + pestic1997_logstd +
                                urbmdhi_logstd,
                              data = df, coords = cbind(df$easting, df$northing),
                              NNGP = T,
                              priors = list(beta_intercept_scale = 100, beta_scale = 100, beta_df = Inf,
                                            phi_lb = 0.005, phi_ub = 0.005),
                              nngp.control = list(n.neighbors = 15, ord = o),
                              nburn = 1000, nsave = 5000)
  #saveRDS(fit_micobin1, file = "fit_micobin1_n949.rds")
  set.seed(22)
  fit_micobin2 = spmicobinreg(MMI_BENT ~ agkffact_logstd + bfi_logstd +
                                cbnf_logstd + conif_logstd +
                                crophay_logstd + fert_logstd +
                                manure_logstd + pestic1997_logstd +
                                urbmdhi_logstd,
                              data = df, coords = cbind(df$easting, df$northing),
                              NNGP = T,
                              priors = list(beta_intercept_scale = 100, beta_scale = 100, beta_df = Inf,
                                            phi_lb = 0.005, phi_ub = 0.005),
                              nngp.control = list(n.neighbors = 15, ord = o),
                              nburn = 1000, nsave = 5000)
  #saveRDS(fit_micobin2, file = "fit_micobin2_n949.rds")
  set.seed(222)
  fit_micobin3 = spmicobinreg(MMI_BENT ~ agkffact_logstd + bfi_logstd +
                                cbnf_logstd + conif_logstd +
                                crophay_logstd + fert_logstd +
                                manure_logstd + pestic1997_logstd +
                                urbmdhi_logstd,
                              data = df, coords = cbind(df$easting, df$northing),
                              NNGP = T,
                              priors = list(beta_intercept_scale = 100, beta_scale = 100, beta_df = Inf,
                                            phi_lb = 0.005, phi_ub = 0.005),
                              nngp.control = list(n.neighbors = 15, ord = o),
                              nburn = 1000, nsave = 5000)
  #saveRDS(fit_micobin3, file = "fit_micobin3_n949.rds")
}else{
  fit_micobin1 = readRDS("fit_micobin1_n947.rds")
  fit_micobin2 = readRDS("fit_micobin2_n947.rds")
  fit_micobin3 = readRDS("fit_micobin3_n947.rds")
}


fit_micobin_betasave = as.mcmc.list(list(fit_micobin1$post_save[,1:ncol(fit_micobin1$X)],
                                         fit_micobin2$post_save[,1:ncol(fit_micobin2$X)],
                                         fit_micobin3$post_save[,1:ncol(fit_micobin3$X)]))
coda::gelman.diag(fit_micobin_betasave)
bayesplot::mcmc_trace(fit_micobin_betasave)# based on standardized covariates


mean(c(as.numeric(mcmcse::multiESS(fit_micobin_betasave[[1]])),
       as.numeric(mcmcse::multiESS(fit_micobin_betasave[[2]])),
       as.numeric(mcmcse::multiESS(fit_micobin_betasave[[3]]))))
mean(c(as.numeric(fit_micobin1$t_mcmc),
       as.numeric(fit_micobin2$t_mcmc),
       as.numeric(fit_micobin3$t_mcmc)))/60


fit_micobin_betasave_all = as.matrix(fit_micobin_betasave)
fit_micobin_betasave_all_orig = coda::mcmc(t(apply(fit_micobin_betasave_all, 1, scaleback, centers = x_center,
                                                   scales = x_sd)))

mycolnames = c("Intercept", "agkffact", "bfi", "cbnf", "conif", "crophay", "fert", "manure", "pestic1997", "urbmdhi")

colnames(fit_micobin_betasave_all_orig) = mycolnames
summary(fit_micobin_betasave_all_orig)




# step 3.3. beta

# Stan + NNGP, following https://mc-stan.org/learn-stan/case-studies/nngp.html
# also see https://github.com/LuZhangstat/NNGP_STAN
source("https://raw.githubusercontent.com/LuZhangstat/NNGP_STAN/refs/heads/master/NNMatrix.R")


M = 15                 # Number of Nearest Neighbors
NN.matrix <- NNMatrix(coords = cbind(df$easting, df$northing), n.neighbors = M, n.omp.threads = 1)
str(NN.matrix)
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

RUN = TRUE
# from ?rstan::stan, even if multiple chains are used, only one seed is specified for parallel run
# for reproducability purpose I will run one chain at a time 
if(RUN){
  fit_beta1 = rstan::stan(file = "betareg_spatial_NNGP_normal.stan",
                          data = stan_data,
                          chains = 1, seed = 3,
                          iter = nburn + nsave, warmup = nburn, thin = nthin)
  #saveRDS(fit_beta1, file = "fit_beta1_n949.rds")
  fit_beta2 = rstan::stan(file = "betareg_spatial_NNGP_normal.stan",
                          data = stan_data,
                          chains = 1, seed = 33,
                          iter = nburn + nsave, warmup = nburn, thin = nthin)
  #saveRDS(fit_beta2, file = "fit_beta2_n949.rds")
  fit_beta3 = rstan::stan(file = "betareg_spatial_NNGP_normal.stan",
                          data = stan_data,
                          chains = 1, seed = 333,
                          iter = nburn + nsave, warmup = nburn, thin = nthin)
  #saveRDS(fit_beta3, file = "fit_beta3_n949.rds")
}else{
  fit_beta1 = readRDS("fit_beta1_n947.rds")
  fit_beta2 = readRDS("fit_beta2_n947.rds")
  fit_beta3 = readRDS("fit_beta3_n947.rds")
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
coda::gelman.diag(fit_beta_betasave)
bayesplot::mcmc_trace(fit_beta_betasave)# based on standardized covariates



X = as.matrix(stan_data$X[order(o),]) # stan NNGP implementation uses permuted data; permute back
all.equal(X, fit_cobin1$X, check.attributes = FALSE)
all.equal(X, fit_micobin1$X, check.attributes = FALSE)


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


mean(c(sum(get_elapsed_time(fit_beta1)),
       sum(get_elapsed_time(fit_beta2)),
       sum(get_elapsed_time(fit_beta3))))/60

mean(c(as.numeric(mcmcse::multiESS(fit_beta1_betasave)),
       as.numeric(mcmcse::multiESS(fit_beta2_betasave)),
       as.numeric(mcmcse::multiESS(fit_beta3_betasave))))

loo::waic(rbind(fit_beta1_logliksave,
                fit_beta2_logliksave,
                fit_beta3_logliksave))

fit_beta_betasave_all = as.matrix(fit_beta_betasave)
fit_beta_betasave_all_orig = coda::mcmc(t(apply(fit_beta_betasave_all, 1, scaleback, centers = x_center,
                                                scales = x_sd)))

mycolnames = c("Intercept", "agkffact", "bfi", "cbnf", "conif", "crophay", "fert", "manure", "pestic1997", "urbmdhi")

colnames(fit_beta_betasave_all_orig) = mycolnames
summary(fit_beta_betasave_all_orig)



