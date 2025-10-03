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
# 
# # step 1: remove one row with Y = 0
# df = df[df$MMI_BENT != 0,]
# dim(df) # nrow = 949

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

dbetarec_cobit <- function(x, linpred, phi, alpha, log = FALSE){
  gam = cobin::bftprime(linpred)
  theta = alpha*(1-abs(2*gam-1))
  mu = (gam - 0.5*theta)/(1-theta)
  loglik = log( ((1-theta) * dbeta(x, mu * phi, (1 - mu) * phi)) + (theta * dunif(x, 0, 1)) )
  if(log) return(loglik) else return(exp(loglik))
}

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
                         empty = TRUE) )

stan_data$phi_spatial = 0.005
stan_data$M = M
stan_data$NN_ind = NN.matrix$NN_ind
stan_data$NN_dist = NN.matrix$NN_dist
stan_data$NN_distM = NN.matrix$NN_distM

RUN = F
if(RUN){
  fit_betarec1 = rstan::stan(file = "betarecreg_spatial_NNGP.stan",
                             data = stan_data,
                             chains = 1, seed = 4,
                             iter = nburn + nsave, warmup = nburn, thin = nthin)
  saveRDS(fit_betarec1, file = "res/fit_betarec1_n950.rds")
  fit_betarec2 = rstan::stan(file = "betarecreg_spatial_NNGP.stan",
                             data = stan_data,
                             chains = 1, seed = 44,
                             iter = nburn + nsave, warmup = nburn, thin = nthin)
  saveRDS(fit_betarec2, file = "res/fit_betarec2_n950.rds")
  fit_betarec3 = rstan::stan(file = "betarecreg_spatial_NNGP.stan",
                             data = stan_data,
                             chains = 1, seed = 444,
                             iter = nburn + nsave, warmup = nburn, thin = nthin)
  saveRDS(fit_betarec3, file = "res/fit_betarec3_n950.rds")
}else{
  fit_betarec1 = readRDS("res/fit_betarec1_n950.rds")
  fit_betarec2 = readRDS("res/fit_betarec2_n950.rds")
  fit_betarec3 = readRDS("res/fit_betarec3_n950.rds")
}

out1 = rstan::extract(fit_betarec1, permute = FALSE)
fit_betarec1_betasave = cbind(out1[,1,"b_Intercept"] ,out1[,1,grep("^b\\[\\d+\\]$", dimnames(out1)$parameters)])
fit_betarec1_usave = out1[,1,grep("^u\\[\\d+\\]$", dimnames(out1)$parameters)]
fit_betarec1_phisave = out1[,1,"phi"]
fit_betarec1_alphasave = out1[,1,"alpha"]


out2 = rstan::extract(fit_betarec2, permute = FALSE)
fit_betarec2_betasave = cbind(out2[,1,"b_Intercept"] ,out2[,1,grep("^b\\[\\d+\\]$", dimnames(out2)$parameters)])
fit_betarec2_usave = out2[,1,grep("^u\\[\\d+\\]$", dimnames(out2)$parameters)]
fit_betarec2_phisave = out2[,1,"phi"]
fit_betarec2_alphasave = out2[,1,"alpha"]

out3 = rstan::extract(fit_betarec3, permute = FALSE)
fit_betarec3_betasave = cbind(out3[,1,"b_Intercept"] ,out3[,1,grep("^b\\[\\d+\\]$", dimnames(out3)$parameters)])
fit_betarec3_usave = out3[,1,grep("^u\\[\\d+\\]$", dimnames(out3)$parameters)]
fit_betarec3_phisave = out3[,1,"phi"]
fit_betarec3_alphasave = out3[,1,"alpha"]

mycolnames = c("Intercept", "agkffact", "bfi", "cbnf", "conif", "crophay", "fert", "manure", "pestic1997", "urbmdhi")
colnames(fit_betarec1_betasave) = mycolnames
colnames(fit_betarec2_betasave) = mycolnames
colnames(fit_betarec3_betasave) = mycolnames

fit_betarec_betasave = as.mcmc.list(list(as.mcmc(fit_betarec1_betasave),
                                         as.mcmc(fit_betarec2_betasave),
                                         as.mcmc(fit_betarec3_betasave)))


X = as.matrix(stan_data$X[order(o),]) # stan NNGP implementation uses permuted data; permute back


fit_betarec1_linpredsave = matrix(0, nsave, nrow(df))
fit_betarec1_logliksave = matrix(0, nsave, nrow(df))
fit_betarec2_linpredsave = matrix(0, nsave, nrow(df))
fit_betarec2_logliksave = matrix(0, nsave, nrow(df))
fit_betarec3_linpredsave = matrix(0, nsave, nrow(df))
fit_betarec3_logliksave = matrix(0, nsave, nrow(df))
for(isave in 1:nsave){
  fit_betarec1_linpredsave[isave,] = as.numeric(X %*% fit_betarec1_betasave[isave,] + fit_betarec1_usave[isave,order(o)]) # permute back random effect
  fit_betarec1_logliksave[isave,] = dbetarec_cobit(df$MMI_BENT, fit_betarec1_linpredsave[isave,], fit_betarec1_phisave[isave], fit_betarec1_alphasave[isave], log = TRUE)
  
  
  fit_betarec2_linpredsave[isave,] = as.numeric(X %*% fit_betarec2_betasave[isave,] + fit_betarec2_usave[isave,order(o)]) # permute back random effect
  fit_betarec2_logliksave[isave,] = dbetarec_cobit(df$MMI_BENT, fit_betarec2_linpredsave[isave,], fit_betarec2_phisave[isave], fit_betarec2_alphasave[isave], log = TRUE)
  
  
  fit_betarec3_linpredsave[isave,] = as.numeric(X %*% fit_betarec3_betasave[isave,] + fit_betarec3_usave[isave,order(o)]) # permute back random effect
  fit_betarec3_logliksave[isave,] = dbetarec_cobit(df$MMI_BENT, fit_betarec3_linpredsave[isave,], fit_betarec3_phisave[isave], fit_betarec3_alphasave[isave], log = TRUE)
}


mean(c(as.numeric(mcmcse::multiESS(fit_betarec1_betasave)),
       as.numeric(mcmcse::multiESS(fit_betarec2_betasave)),
       as.numeric(mcmcse::multiESS(fit_betarec3_betasave))))

mean(c(sum(get_elapsed_time(fit_betarec1)),
       sum(get_elapsed_time(fit_betarec2)),
       sum(get_elapsed_time(fit_betarec3))))/60

loo::loo(rbind(fit_betarec1_logliksave,
               fit_betarec2_logliksave,
               fit_betarec3_logliksave))

loo::waic(rbind(fit_betarec1_logliksave,
                fit_betarec2_logliksave,
                fit_betarec3_logliksave))

fit_betarec_betasave_all = as.matrix(fit_betarec_betasave)
fit_betarec_betasave_all_orig = coda::mcmc(t(apply(fit_betarec_betasave_all, 1, scaleback, centers = x_center,
                                                   scales = x_sd)))

mycolnames = c("Intercept", "agkffact", "bfi", "cbnf", "conif", "crophay", "fert", "manure", "pestic1997", "urbmdhi")

colnames(fit_betarec_betasave_all_orig) = mycolnames
summary(fit_betarec_betasave_all_orig)


fit_betarec_parsave <- mcmc.list(
  as.mcmc(cbind(
    as.matrix(fit_betarec_betasave[[1]]),
    phi  = out1[,1,"phi"],
    sigmasq = out1[,1,"sigmasq"],
    loglik = rowSums(fit_betarec1_logliksave)
  )),
  as.mcmc(cbind(
    as.matrix(fit_betarec_betasave[[2]]),
    phi  = out2[,1,"phi"],
    sigmasq = out2[,1,"sigmasq"],
    loglik = rowSums(fit_betarec2_logliksave)
  )),
  as.mcmc(cbind(
    as.matrix(fit_betarec_betasave[[3]]),
    phi  = out3[,1,"phi"],
    sigmasq = out3[,1,"sigmasq"],
    loglik = rowSums(fit_betarec3_logliksave)
  ))
)

coda::gelman.diag(fit_betarec_parsave)




# average marginal effect, logit link
nsave = 15000
fit_betarec_linpredsave_all = rbind(fit_betarec1_linpredsave, fit_betarec2_linpredsave, fit_betarec3_linpredsave)
avgslope_betarec_save = matrix(0, nsave, 9)
colnames(avgslope_betarec_save) = mycolnames[-1]
for(isave in 1:nsave){
  for(j in 1:9){
    avgslope_betarec_save[isave,j] = mean( cobin::bftprimeprime(fit_betarec_linpredsave_all[isave, ]) * fit_betarec_betasave_all_orig[isave,j+1])
  }
}
round(colMeans(avgslope_betarec_save), 3)
round(apply(avgslope_betarec_save, 2, function(x) quantile(x, 0.025)), 3)
round(apply(avgslope_betarec_save, 2, function(x) quantile(x, 0.975)), 3)


#fields::quilt.plot(df$lon, df$lat, colMeans(fit_cobin_linpredsave_all))
tt = cbind(round(colMeans(avgslope_betarec_save), 3), 
           round(apply(avgslope_betarec_save, 2, function(x) quantile(x, 0.025)), 3),
           round(apply(avgslope_betarec_save, 2, function(x) quantile(x, 0.975)), 3))

formatted <- apply(tt, 1, function(x) {
  sprintf("%.3f (%.3f, %.3f)", x[1], x[2], x[3])
})

formatted

fit_betarec_linpred_mean = colMeans(fit_betarec_linpredsave_all)
fit_betarec_linpred_sd = apply(fit_betarec_linpredsave_all, 2, sd)

fit_betarec_mu_mean = colMeans(cobin::bftprime(fit_betarec_linpredsave_all))
fit_betarec_mu_sd = apply(cobin::bftprime(fit_betarec_linpredsave_all), 2, sd)


saveRDS(list(
  coef_std_mean = colMeans(fit_betarec_betasave_all),
  coef_std_sd = apply(fit_betarec_betasave_all, 2, sd),
  coef_orig_mean = colMeans(fit_betarec_betasave_all_orig),
  coef_orig_sd = apply(fit_betarec_betasave_all_orig, 2, sd),
  coef_orig_95lo = apply(fit_betarec_betasave_all_orig, 2, function(x) quantile(x, 0.025)),
  coef_orig_95up = apply(fit_betarec_betasave_all_orig, 2, function(x) quantile(x, 0.975)),
  phi_mean = mean(c(fit_betarec1_phisave, fit_betarec2_phisave, fit_betarec3_phisave)),
  alpha_mean = mean(c(fit_betarec1_alphasave, fit_betarec2_alphasave, fit_betarec3_alphasave)),
  linpred_mean = fit_betarec_linpred_mean,
  linpred_sd = fit_betarec_linpred_sd,
  mu_mean = fit_betarec_mu_mean,
  mu_sd = fit_betarec_mu_sd),
  file = "res/summary_betarec_n950.rds")


