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
RUN = F
o = order(df$easting) # default order for NNGP
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
  saveRDS(fit_cobin1, file = "res/fit_cobin1_n947.rds")
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
  saveRDS(fit_cobin2, file = "res/fit_cobin2_n947.rds")
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
  saveRDS(fit_cobin3, file = "res/fit_cobin3_n947.rds")
}else{
  fit_cobin1 = readRDS("res/fit_cobin1_n947.rds")
  fit_cobin2 = readRDS("res/fit_cobin2_n947.rds")
  fit_cobin3 = readRDS("res/fit_cobin3_n947.rds")
}


fit_cobin_betasave = as.mcmc.list(list(fit_cobin1$post_save[,1:ncol(fit_cobin1$X)],
                                       fit_cobin2$post_save[,1:ncol(fit_cobin2$X)],
                                       fit_cobin3$post_save[,1:ncol(fit_cobin3$X)]))
library(coda)

fit_cobin_parsave <- mcmc.list(
  as.mcmc(cbind(
    as.matrix(fit_cobin_betasave[[1]]),
    lambda  = fit_cobin1$post_save[, "lambda",    drop = FALSE],
    sigmasq = fit_cobin1$post_save[, "sigma.sq",  drop = FALSE],
    loglik = rowSums(fit_cobin1$loglik_save)
  )),
  as.mcmc(cbind(
    as.matrix(fit_cobin_betasave[[2]]),
    lambda  = fit_cobin2$post_save[, "lambda",    drop = FALSE],
    sigmasq = fit_cobin2$post_save[, "sigma.sq",  drop = FALSE],
    loglik = rowSums(fit_cobin2$loglik_save)
  )),
  as.mcmc(cbind(
    as.matrix(fit_cobin_betasave[[3]]),
    lambda  = fit_cobin3$post_save[, "lambda",    drop = FALSE],
    sigmasq = fit_cobin3$post_save[, "sigma.sq",  drop = FALSE],
    loglik = rowSums(fit_cobin3$loglik_save)
  ))
)
coda::gelman.diag(fit_cobin_parsave) 

# mESS
mean(c(as.numeric(mcmcse::multiESS(fit_cobin_betasave[[1]])),
       as.numeric(mcmcse::multiESS(fit_cobin_betasave[[2]])),
       as.numeric(mcmcse::multiESS(fit_cobin_betasave[[3]]))))
# time
mean(c(as.numeric(fit_cobin1$t_mcmc),
       as.numeric(fit_cobin2$t_mcmc),
       as.numeric(fit_cobin3$t_mcmc)))/60 # in min

# PSIS-LOO
loo::loo(rbind(fit_cobin1$loglik_save,
               fit_cobin2$loglik_save,
               fit_cobin3$loglik_save))
# WAIC
loo::waic(rbind(fit_cobin1$loglik_save,
                fit_cobin2$loglik_save,
                fit_cobin3$loglik_save))

fit_cobin_betasave_all = as.matrix(fit_cobin_betasave)
fit_cobin_betasave_all_orig = coda::mcmc(t(apply(fit_cobin_betasave_all, 1, scaleback, centers = x_center,
                                                 scales = x_sd)))
mycolnames = c("Intercept", "agkffact", "bfi", "cbnf", "conif", "crophay", "fert", "manure", "pestic1997", "urbmdhi")

colnames(fit_cobin_betasave_all_orig) = mycolnames

# Table 5 entries
summary(fit_cobin_betasave_all_orig)



# save linear predictor 
fit_cobin_usave_all = rbind(fit_cobin1$post_u_save, 
                            fit_cobin2$post_u_save, 
                            fit_cobin3$post_u_save)
str(fit_cobin_usave_all)

# average marginal effect, logit link
nsave = 15000
fit_cobin_linpredsave_all = matrix(0, dim(fit_cobin_usave_all)[1], dim(fit_cobin_usave_all)[2])
avgslope_cobin_save = matrix(0, nsave, 9)
colnames(avgslope_cobin_save) = mycolnames[-1]
for(isave in 1:nsave){
  fit_cobin_linpredsave_all[isave, ] = as.numeric(fit_cobin1$X %*% fit_cobin_betasave_all[isave,] + fit_cobin_usave_all[isave,]) # permute back random effect
  for(j in 1:9){
    avgslope_cobin_save[isave,j] = mean( cobin::bftprimeprime(fit_cobin_linpredsave_all[isave, ]) * fit_cobin_betasave_all_orig[isave,j+1])
  }
}

# Table S.9
round(colMeans(avgslope_cobin_save), 3)
round(apply(avgslope_cobin_save, 2, function(x) quantile(x, 0.025)), 3)
round(apply(avgslope_cobin_save, 2, function(x) quantile(x, 0.975)), 3)

# #fields::quilt.plot(df$lon, df$lat, colMeans(fit_cobin_linpredsave_all))
# tt = cbind(round(colMeans(avgslope_cobin_save), 3), 
#            round(apply(avgslope_cobin_save, 2, function(x) quantile(x, 0.025)), 3),
#            round(apply(avgslope_cobin_save, 2, function(x) quantile(x, 0.975)), 3))
# formatted <- apply(tt, 1, function(x) {
#   sprintf("%.3f (%.3f, %.3f)", x[1], x[2], x[3])
# })
# formatted


fit_cobin_linpred_mean = colMeans(fit_cobin_linpredsave_all)
fit_cobin_linpred_sd = apply(fit_cobin_linpredsave_all, 2, sd)

fit_cobin_mu_mean = colMeans(cobin::bftprime(fit_cobin_linpredsave_all))
fit_cobin_mu_sd = apply(cobin::bftprime(fit_cobin_linpredsave_all), 2, sd)


# save


saveRDS(list(
  coef_std_mean = colMeans(fit_cobin_betasave_all),
  coef_std_sd = apply(fit_cobin_betasave_all, 2, sd),
  coef_orig_mean = colMeans(fit_cobin_betasave_all_orig),
  coef_orig_sd = apply(fit_cobin_betasave_all_orig, 2, sd),
  coef_orig_95lo = apply(fit_cobin_betasave_all_orig, 2, function(x) quantile(x, 0.025)),
  coef_orig_95up = apply(fit_cobin_betasave_all_orig, 2, function(x) quantile(x, 0.975)),
  lambda_median = median(c(fit_cobin1$post_save[,"lambda"],
                         fit_cobin2$post_save[,"lambda"],
                         fit_cobin3$post_save[,"lambda"])),
  linpred_mean = fit_cobin_linpred_mean,
  linpred_sd = fit_cobin_linpred_sd,
  mu_mean = fit_cobin_mu_mean,
  mu_sd = fit_cobin_mu_sd),
  file = "res/summary_cobin_n947.rds")






