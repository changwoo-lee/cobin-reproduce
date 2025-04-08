#rm(list = ls())
# set path as current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df = read.csv("../mmi_lakecat.csv")
library(dplyr)
library(cobin)
packageVersion("cobin") # > 0.0.1.5
library(spam)
library(spNNGP)
# step 1: remove one row with Y = 0
#df = df[df$MMI_BENT != 0,]
# commented out since it is using original data
dim(df) # nrow = 950

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
#saveRDS(fit_micobin1, file = "fit_micobin1_n950.rds")
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
#saveRDS(fit_micobin2, file = "fit_micobin2_n950.rds")
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
#saveRDS(fit_micobin3, file = "fit_micobin3_n950.rds")
}else{
  fit_micobin1 = readRDS("fit_micobin1_n950.rds")
  fit_micobin2 = readRDS("fit_micobin2_n950.rds")
  fit_micobin3 = readRDS("fit_micobin3_n950.rds")
}

library(coda)
library(mcmcse)
library(loo)
library(bayesplot)

fit_micobin_betasave = as.mcmc.list(list(fit_micobin1$post_save[,1:ncol(fit_micobin1$X)],
                                      fit_micobin2$post_save[,1:ncol(fit_micobin2$X)],
                                      fit_micobin3$post_save[,1:ncol(fit_micobin3$X)]))
coda::gelman.diag(fit_micobin_betasave)
bayesplot::mcmc_trace(fit_micobin_betasave) # based on standardized covariates


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

#colnames(fit_beta_betasave_all_orig) = colnames(fit_cobin1$X)
#colnames(fit_cobin_betasave_all_orig) = colnames(fit_cobin1$X)
colnames(fit_micobin_betasave_all_orig) = mycolnames

summary(fit_micobin_betasave_all_orig)




# 
# #fit_micobin1 = readRDS("real2/fit_micobin1.rds")
# #fit_micobin2 = readRDS("../real2/fit_micobin2.rds")
# #fit_micobin3 = readRDS("../real2/fit_micobin3.rds")
# fit_micobin_betasave = as.mcmc.list(list(fit_micobin1$post_save[,1:ncol(fit_micobin1$X)],
#                                        fit_micobin2$post_save[,1:ncol(fit_micobin2$X)],
#                                        fit_micobin3$post_save[,1:ncol(fit_micobin3$X)]))
# 
# 
# mean(c(as.numeric(fit_micobin1$t_mcmc),
#        as.numeric(fit_micobin2$t_mcmc),
#        as.numeric(fit_micobin3$t_mcmc)))/60
# 

# 
# 
# coda::gelman.diag(fit_micobin_betasave)
# 
# bayesplot::mcmc_trace(fit_cobin_betasave)
# 
# loo::waic(rbind(fit_micobin1$loglik_save,
#                 fit_micobin2$loglik_save,
#                 fit_micobin3$loglik_save))
# 
# fit_micobin_betahat = colMeans(as.matrix(fit_micobin_betasave))
# fit_micobin_linpredhat = fit_micobin1$X%*%fit_micobin_betahat + colMeans(rbind(fit_micobin1$post_u_save, fit_micobin2$post_u_save, fit_micobin3$post_u_save))
# fit_micobin_muhat = cobin::bftprime(fit_micobin_linpredhat)
# fit_micobin_psihat = mean(c(fit_micobin1$post_save[,"psi"], fit_micobin2$post_save[,"psi"],fit_micobin3$post_save[,"psi"] ))
# 
# # quantile residual
# micobin_qresid = numeric(nrow(df))
# for(i in 1:nrow(df)){
#   micobin_qresid[i] = cobin::pmicobin(df$MMI_BENT[i], fit_micobin_linpredhat[i], fit_micobin_psihat)
# }
# 
# 
# qqnorm(qnorm(micobin_qresid), main = "Quantile residual, micobin", pch = 16)
# qqline(qnorm(micobin_qresid), col = 2)
# 
# 
# dfdf <- data.frame(
#   model = rep(c("micobin regression"), each = length(micobin_qresid)),
#   qresid = qnorm(micobin_qresid)
# )
# 
# # Generate the QQ plots using facet_wrap
# figqq = ggplot(dfdf, aes(sample = qresid)) +
#   stat_qq(shape = 16) +
#   stat_qq_line(color = "red") +
#   facet_wrap(~ model) +
#   labs(x = "Theoretical Quantiles",
#        y = "Sample Quantiles") +
#   theme_bw()
# library(dplyr)
# library(ggplot2)
# 
# n <- nrow(df)
# dfdf <- data.frame(
#   model = rep("micobin regression", each = n),
#   qresid = qnorm(micobin_qresid),
#   obs_orig = rep(1:n)
# )
# 
# # Compute theoretical quantiles by grouping and sorting by qresid,
# # but keep the original observation number in 'obs_orig'
# dfdf <- dfdf %>%
#   group_by(model) %>%
#   arrange(qresid, .by_group = TRUE) %>%
#   mutate(prob = (row_number() - 0.5) / n(),
#          theo = qnorm(prob)) %>%
#   ungroup()
# 
# # Create the QQ plot and label extreme points with their original observation number
# figqq <- ggplot(dfdf, aes(sample = qresid)) +
#   stat_qq(shape = 16) +
#   stat_qq_line(color = "red") +
#   geom_text(aes(x = theo, y = qresid,
#                 label = ifelse(theo < -2.8, obs_orig, "")),
#             hjust = +0.8, vjust = -1, size = 2.5) +
#   facet_wrap(~ model) +
#   labs(x = "Theoretical Quantiles",
#        y = "Sample Quantiles") +
#   theme_bw()
# 
# print(figqq)
# 
# df[147,]
# fit_beta_muhat[147]
# fit_cobin_muhat[147]
# fit_micobin_muhat[147]
# 
# df[34,]
# fit_beta_muhat[34]
# fit_cobin_muhat[34]
# fit_micobin_muhat[34]
# 
# df[945,]
# fit_beta_muhat[945]
# fit_cobin_muhat[945]
# fit_micobin_muhat[945]
# 
# 
# 
# 
# 
# ggsave("quantile_qqplot2.png", figqq, width = 9, height = 3)
# ggsave("quantile_qqplot2.pdf", figqq, width = 9, height = 3)
# 
# 
# 
# scaleback <- function(x, centers, scales){
#   as.numeric(c(x[1] - sum((x[-1] * centers) / scales),
#     x[-1] / scales))
# }
# 
# x_center = as.numeric(x_center)
# x_sd = as.numeric(x_sd)
# 
# fit_micobin_betasave_all = as.matrix(fit_micobin_betasave)
# 
# 
# fit_micobin_betasave_all_orig = coda::mcmc(t(apply(fit_micobin_betasave_all, 1, scaleback, centers = x_center,
#                                      scales = x_sd)))
# 
# #colnames(fit_beta_betasave_all_orig) = colnames(fit_cobin1$X)
# #colnames(fit_cobin_betasave_all_orig) = colnames(fit_cobin1$X)
# colnames(fit_micobin_betasave_all_orig) = colnames(fit_micobin1$X)
# 
# summary(fit_micobin_betasave_all_orig)
# 
# 
# #####################################################################################
# 
# 
# 
# str(fit_micobin1$spNNGPfit)
# fit_micobin1$spNNGPfit$X = fit_micobin1$X
# fit_micobin1$spNNGPfit$p.beta.samples = fit_micobin1$post_save[,1:ncol(fit_micobin1$X)]
# fit_micobin1$spNNGPfit$p.theta.samples = fit_micobin1$post_save[,c("sigma.sq", "phi")]
# fit_micobin1$spNNGPfit$p.w.samples = t(as.matrix(fit_micobin1$post_u_save))
# 
# df4 = read.csv("../real2/predictor_40000m2.csv")
# 
# # remove the duplicate
# 
# 
# df4 = df4 %>% mutate(agkffact_logstd = scale(log2(1+agkffact), center = x_center[1], scale = x_sd[1]),
#                      bfi_logstd = scale(log2(1+bfi), center = x_center[2], scale = x_sd[2]),
#                      cbnf_logstd = scale(log2(1+cbnf), center = x_center[3], scale = x_sd[3]),
#                      conif_logstd = scale(log2(1+conif), center = x_center[4], scale = x_sd[4]),
#                      crophay_logstd = scale(log2(1+crophay), center = x_center[5], scale = x_sd[5]),
#                      fert_logstd = scale(log2(1+fert), center = x_center[6], scale = x_sd[6]),
#                      manure_logstd = scale(log2(1+manure), center = x_center[7], scale = x_sd[7]),
#                      pestic1997_logstd = scale(log2(1+pestic1997), center = x_center[8], scale = x_sd[8]),
#                      urbmdhi_logstd = scale(log2(1+urbmdhi), center = x_center[9], scale = x_sd[9])
# )
# # select
# Xnew = data.frame(
#   df4$agkffact_logstd,
#   df4$bfi_logstd,
#   df4$cbnf_logstd,
#   df4$conif_logstd,
#   df4$crophay_logstd,
#   df4$fert_logstd,
#   df4$manure_logstd,
#   df4$pestic1997_logstd,
#   df4$urbmdhi_logstd
# )
# Xnew = as.matrix(cbind(1, Xnew))
# str(Xnew)
# 
# pred_micobin1 <- predict(fit_micobin1$spNNGPfit, X.0 = Xnew, coords.0 = cbind(df4$easting,df4$northing),
#                   sub.sample = list(start = 1, end = 5000, thin = 1))
# # replace NaN to NA: likely due to numerical issue with too close coordinates
# pred_micobin1$p.w.0[is.nan(pred_micobin1$p.w.0)] = NA
# pred_micobin1$p.y.0[is.nan(pred_micobin1$p.y.0)] = NA
# 
# wnew_hat = rowMeans(pred_micobin1$p.w.0, na.rm = TRUE)
# wnew_sd = apply(pred_micobin1$p.w.0, 1, sd, na.rm = TRUE)
# 
# munew_hat = rowMeans(cobin::bftprime(pred_micobin1$p.y.0), na.rm = TRUE)
# munew_sd = apply(cobin::bftprime(pred_micobin1$p.y.0), 1, sd, na.rm = TRUE)
# 
# str(wnew_hat)
# fields::quilt.plot(cbind(df4$easting,df4$northing), wnew_hat)
# fields::quilt.plot(cbind(df4$easting,df4$northing), wnew_sd)
# 
# fields::quilt.plot(cbind(df4$easting,df4$northing), munew_hat)
# fields::quilt.plot(cbind(df4$easting,df4$northing), munew_sd)
# 
# 
# saveRDS(list(wnew_hat = wnew_hat,
#              wnew_sd = wnew_sd,
#              munew_hat = munew_hat,
#              munew_sd = munew_sd), file = "real2/pred_micobin1.rds")
# 
# ####################################################
# 
# spNNGP_beta_proxy = fit_micobin1$spNNGPfit
# spNNGP_beta_proxy$X = fit_micobin1$X
# spNNGP_beta_proxy$p.beta.samples = fit_beta1_betasave  #fit_micobin1$post_save[,1:ncol(fit_micobin1$X)]
# spNNGP_beta_proxy$p.theta.samples = cbind(out[,1,"sigmasq"], 0.005)
# spNNGP_beta_proxy$p.w.samples = t(as.matrix(fit_beta1_usave[,order(o)]))
# 
# 
# #df4 = read.csv("real2/predictor_40000m2.csv")
# # remove the duplicate
# #
# #
# # df4 = df4 %>% mutate(agkffact_logstd = scale(log2(1+agkffact), center = x_center[1], scale = x_sd[1]),
# #                      bfi_logstd = scale(log2(1+bfi), center = x_center[2], scale = x_sd[2]),
# #                      cbnf_logstd = scale(log2(1+cbnf), center = x_center[3], scale = x_sd[3]),
# #                      conif_logstd = scale(log2(1+conif), center = x_center[4], scale = x_sd[4]),
# #                      crophay_logstd = scale(log2(1+crophay), center = x_center[5], scale = x_sd[5]),
# #                      fert_logstd = scale(log2(1+fert), center = x_center[6], scale = x_sd[6]),
# #                      manure_logstd = scale(log2(1+manure), center = x_center[7], scale = x_sd[7]),
# #                      pestic1997_logstd = scale(log2(1+pestic1997), center = x_center[8], scale = x_sd[8]),
# #                      urbmdhi_logstd = scale(log2(1+urbmdhi), center = x_center[9], scale = x_sd[9])
# # )
# # # select
# # Xnew = data.frame(
# #   df4$agkffact_logstd,
# #   df4$bfi_logstd,
# #   df4$cbnf_logstd,
# #   df4$conif_logstd,
# #   df4$crophay_logstd,
# #   df4$fert_logstd,
# #   df4$manure_logstd,
# #   df4$pestic1997_logstd,
# #   df4$urbmdhi_logstd
# # )
# # Xnew = as.matrix(cbind(1, Xnew))
# # str(Xnew)
# 
# pred_beta1 <- predict(spNNGP_beta_proxy, X.0 = Xnew, coords.0 = cbind(df4$easting,df4$northing),
#                          sub.sample = list(start = 1, end = 5000, thin = 1))
# # replace NaN to NA: likely due to numerical issue with too close training coordinates
# pred_beta1$p.w.0[is.nan(pred_beta1$p.w.0)] = NA
# pred_beta1$p.y.0[is.nan(pred_beta1$p.y.0)] = NA
# 
# wnew_hat = rowMeans(pred_beta1$p.w.0, na.rm = TRUE)
# wnew_sd = apply(pred_beta1$p.w.0, 1, sd, na.rm = TRUE)
# 
# munew_hat = rowMeans(cobin::bftprime(pred_beta1$p.y.0), na.rm = TRUE)
# munew_sd = apply(cobin::bftprime(pred_beta1$p.y.0), 1, sd, na.rm = TRUE)
# 
# str(wnew_hat)
# fields::quilt.plot(cbind(df4$easting,df4$northing), wnew_hat)
# fields::quilt.plot(cbind(df4$easting,df4$northing), wnew_sd)
# 
# fields::quilt.plot(cbind(df4$easting,df4$northing), munew_hat)
# fields::quilt.plot(cbind(df4$easting,df4$northing), munew_sd)
# 
# 
# saveRDS(list(wnew_hat = wnew_hat,
#              wnew_sd = wnew_sd,
#              munew_hat = munew_hat,
#              munew_sd = munew_sd), file = "real2/pred_beta1.rds")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# save.image("real2/real2_run.RData")
# out = extract(fit_beta, permute = FALSE)
# fitbeta_u_save = out[,1,grep("^u\\[\\d+\\]$", dimnames(out)$parameters)]
# 
# ((1:954)[o])[order(o)]
# 
# fields::quilt.plot(cbind(df$easting, df$northing), colMeans(fitbeta_u_save[,order(o)]))
# 
# fitbeta_save = cbind(out[,1,"b_Intercept"] ,out[,1,grep("^b\\[\\d+\\]$", dimnames(out)$parameters)])
# coda::effectiveSize(fitbeta_save)
# plot(fitbeta_save[,1],type="l")
# plot(fit_micobin$post_save[,1], type="l")
# plot(fitbeta_save[,1],type="l")
# 
# bayesplot::mcmc_intervals(betafit1, pars = c("b[1]", "b[2]", "b[3]", "b[4]", "b[5]", "b[6]", "b[7]"))
# plot(out[,1,"phi"] , type="l")
# 
# 
# 
# 
# out = extract(fit_beta, permute = FALSE)
# outt = extract(fit_beta)
# 
# fit_beta_betasave = cbind(out[,1,"b_Intercept"] ,out[,1,grep("^b\\[\\d+\\]$", dimnames(out)$parameters)])
# fit_beta_betasave = cbind(outt$b_Intercept,outt$b)
# fit_beta_usave = out[,1,grep("^u\\[\\d+\\]$", dimnames(out)$parameters)]
# fit_beta_usave = outt$u
# 
# fit_beta_phisave = out[,1,"phi"]
# fit_beta_phisave = outt$phi
# 
# 
# fields::quilt.plot(cbind(df$easting, df$northing), colMeans(fit_beta_usave[,order(o)]))
# fields::quilt.plot(cbind(df$easting, df$northing), apply(fit_beta_usave[,order(o)], 2, sd))
# 
# summary(fit_beta_phisave)
# 
# X = fit_cobin$X
# fit_cobin_betasave = fit_cobin$post_save[,1:8]
# fit_cobin_quantilesave = matrix(0, nsave, nrow(df))
# 
# fit_micobin_betasave = fit_micobin$post_save[,1:8]
# fit_micobin_quantilesave = matrix(0, nsave, nrow(df))
# 
# fit_beta_linpredsave = matrix(0, nsave, nrow(df))
# fit_beta_logliksave = matrix(0, nsave, nrow(df))
# fit_beta_quantilesave = matrix(0, nsave, nrow(df))
# for(isave in 1:nsave){
#   fit_beta_linpredsave[isave,] = as.numeric(X %*% fit_beta_betasave[isave,]) + as.numeric(fit_beta_usave[isave,order(o)])
#   fit_beta_logliksave[isave,] = betareg::dbetar(df$MMI_BENT,
#                                                 cobin::bftprime(fit_beta_linpredsave[isave,]),
#                                                 fit_beta_phisave[isave], log = TRUE)
# }
# loo::waic(fit_beta_logliksave)
# 
#  plot(df$MMI_BENT, colMeans(cobin::bftprime(fit_beta_linpredsave)))
# 
# fit_beta_betahat = colMeans(fit_beta_betasave)
# fit_beta_uhat = colMeans(fit_beta_usave[,order(o)])
# fit_beta_phihat = mean(fit_beta_phisave)
# 
# X[o]
# 
# 
# ?loo::waic
# loo::waic(fit_beta_logliksave)
# loo::waic(fit_cobin$loglik_save)
# loo::waic(fit_micobin$loglik_save)
