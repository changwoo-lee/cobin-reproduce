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


# step 3 : fit the model
RUN = F
o = order(df$easting) # default order for NNGP
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
saveRDS(fit_micobin1, file = "res/fit_micobin1_n949.rds")
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
saveRDS(fit_micobin2, file = "res/fit_micobin2_n949.rds")
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
saveRDS(fit_micobin3, file = "res/fit_micobin3_n949.rds")
}else{
  fit_micobin1 = readRDS("res/fit_micobin1_n949.rds")
  fit_micobin2 = readRDS("res/fit_micobin2_n949.rds")
  fit_micobin3 = readRDS("res/fit_micobin3_n949.rds")
}


fit_micobin_betasave = as.mcmc.list(list(fit_micobin1$post_save[,1:ncol(fit_micobin1$X)],
                                         fit_micobin2$post_save[,1:ncol(fit_micobin2$X)],
                                         fit_micobin3$post_save[,1:ncol(fit_micobin3$X)]))
library(coda)

fit_micobin_parsave <- mcmc.list(
  as.mcmc(cbind(
    as.matrix(fit_micobin_betasave[[1]]),
    psi  = fit_micobin1$post_save[, "psi",    drop = FALSE],
    sigmasq = fit_micobin1$post_save[, "sigma.sq",  drop = FALSE],
    loglik = rowSums(fit_micobin1$loglik_save)
  )),
  as.mcmc(cbind(
    as.matrix(fit_micobin_betasave[[2]]),
    psi  = fit_micobin2$post_save[, "psi",    drop = FALSE],
    sigmasq = fit_micobin2$post_save[, "sigma.sq",  drop = FALSE],
    loglik = rowSums(fit_micobin2$loglik_save)
  )),
  as.mcmc(cbind(
    as.matrix(fit_micobin_betasave[[3]]),
    psi  = fit_micobin3$post_save[, "psi",    drop = FALSE],
    sigmasq = fit_micobin3$post_save[, "sigma.sq",  drop = FALSE],
    loglik = rowSums(fit_micobin3$loglik_save)
  ))
)

coda::gelman.diag(fit_micobin_parsave) # ignore loglik

mean(c(as.numeric(mcmcse::multiESS(fit_micobin_betasave[[1]])),
       as.numeric(mcmcse::multiESS(fit_micobin_betasave[[2]])),
       as.numeric(mcmcse::multiESS(fit_micobin_betasave[[3]]))))

mean(c(as.numeric(fit_micobin1$t_mcmc),
       as.numeric(fit_micobin2$t_mcmc),
       as.numeric(fit_micobin3$t_mcmc)))/60


loo::loo(rbind(fit_micobin1$loglik_save,
               fit_micobin2$loglik_save,
               fit_micobin3$loglik_save))

loo::waic(rbind(fit_micobin1$loglik_save,
                fit_micobin2$loglik_save,
                fit_micobin3$loglik_save))



fit_micobin_betasave_all = as.matrix(fit_micobin_betasave)
fit_micobin_betasave_all_orig = coda::mcmc(t(apply(fit_micobin_betasave_all, 1, scaleback, centers = x_center,
                                                   scales = x_sd)))

mycolnames = c("Intercept", "agkffact", "bfi", "cbnf", "conif", "crophay", "fert", "manure", "pestic1997", "urbmdhi")

colnames(fit_micobin_betasave_all_orig) = mycolnames
summary(fit_micobin_betasave_all_orig)


# save linear predictor 

fit_micobin_usave_all = rbind(fit_micobin1$post_u_save, 
                              fit_micobin2$post_u_save, 
                              fit_micobin3$post_u_save)
str(fit_micobin_usave_all)

# average marginal effect, logit link
nsave = 15000
fit_micobin_linpredsave_all = matrix(0, dim(fit_micobin_usave_all)[1], dim(fit_micobin_usave_all)[2])
avgslope_micobin_save = matrix(0, nsave, 9)
colnames(avgslope_micobin_save) = mycolnames[-1]
for(isave in 1:nsave){
  fit_micobin_linpredsave_all[isave, ] = as.numeric(fit_micobin1$X %*% fit_micobin_betasave_all[isave,] + fit_micobin_usave_all[isave,]) # permute back random effect
  for(j in 1:9){
    avgslope_micobin_save[isave,j] = mean( cobin::bftprimeprime(fit_micobin_linpredsave_all[isave, ]) * fit_micobin_betasave_all_orig[isave,j+1])
  }
}

tt = cbind(round(colMeans(avgslope_micobin_save), 3), 
           round(apply(avgslope_micobin_save, 2, function(x) quantile(x, 0.025)), 3),
           round(apply(avgslope_micobin_save, 2, function(x) quantile(x, 0.975)), 3))

formatted <- apply(tt, 1, function(x) {
  sprintf("%.3f (%.3f, %.3f)", x[1], x[2], x[3])
})

formatted


fit_micobin_linpred_mean = colMeans(fit_micobin_linpredsave_all)
fit_micobin_linpred_sd = apply(fit_micobin_linpredsave_all, 2, sd)

fit_micobin_mu_mean = colMeans(cobin::bftprime(fit_micobin_linpredsave_all))
fit_micobin_mu_sd = apply(cobin::bftprime(fit_micobin_linpredsave_all), 2, sd)

# save



saveRDS(list(
  coef_std_mean = colMeans(fit_micobin_betasave_all),
  coef_std_sd = apply(fit_micobin_betasave_all, 2, sd),
  coef_orig_mean = colMeans(fit_micobin_betasave_all_orig),
  coef_orig_sd = apply(fit_micobin_betasave_all_orig, 2, sd),
  coef_orig_95lo = apply(fit_micobin_betasave_all_orig, 2, function(x) quantile(x, 0.025)),
  coef_orig_95up = apply(fit_micobin_betasave_all_orig, 2, function(x) quantile(x, 0.975)),
  psi_mean = mean(c(fit_micobin1$post_save[,"psi"],
                         fit_micobin2$post_save[,"psi"],
                         fit_micobin3$post_save[,"psi"])),
  linpred_mean = fit_micobin_linpred_mean,
  linpred_sd = fit_micobin_linpred_sd,
  mu_mean = fit_micobin_mu_mean,
  mu_sd = fit_micobin_mu_sd),
  file = "res/summary_micobin_n949.rds")


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



#fit_micobin1 = readRDS("res/fit_micobin1_n949.rds")

# plug-in posterior samples
fit_micobin1$spNNGPfit$p.beta.samples = fit_micobin1$post_save[,1:ncol(fit_micobin1$X)]
# order important
fit_micobin1$spNNGPfit$p.theta.samples = fit_micobin1$post_save[,c("sigma.sq", "phi")]
fit_micobin1$spNNGPfit$p.w.samples = t(as.matrix(fit_micobin1$post_u_save))

# predict using spNNGP; takes some time
pred_micobin1 <- predict(fit_micobin1$spNNGPfit, 
                       X.0 = Xnew, 
                       coords.0 = cbind(df4$easting,df4$northing),
                       sub.sample = list(start = 1, end = fit_micobin1$nsave, thin = 1))


###################################################

transformed_coords4 <- sf::st_coordinates(usmap::usmap_transform(df4[,c("lon","lat")]))

pred_micobin1$p.w.0[is.nan(pred_micobin1$p.w.0)] = NA
pred_micobin1$p.y.0[is.nan(pred_micobin1$p.y.0)] = NA

wnew_hat_micobin = rowMeans(pred_micobin1$p.w.0, na.rm = TRUE)
wnew_sd_micobin = apply(pred_micobin1$p.w.0, 1, sd, na.rm = TRUE)

munew_hat_micobin = rowMeans(cobin::bftprime(pred_micobin1$p.y.0), na.rm = TRUE)
munew_sd_micobin = apply(cobin::bftprime(pred_micobin1$p.y.0), 1, sd, na.rm = TRUE)



drawdf4_micobin = data.frame(
  x = transformed_coords4[,1],
  y = transformed_coords4[,2],
  urbmdhi = log2(1+df4$urbmdhi),
  wnew_hat = wnew_hat_micobin,
  wnew_sd = wnew_sd_micobin,
  munew_hat = munew_hat_micobin,
  munew_sd = munew_sd_micobin
)

saveRDS(drawdf4_cobin, "res/predict_micobin.rds")





















