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

colnames(fit_micobin_betasave_all_orig) = mycolnames

summary(fit_micobin_betasave_all_orig)


