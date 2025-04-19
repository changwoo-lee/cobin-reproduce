
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

# Load libraries
library(Matrix)
library(fields)
library(spam)
library(mvnfast)
library(coda)
library(cobin)
library(matrixStats)
# load data
######################
# n500, rho 0.2
######################
source("predict_linpred_cobit.R")
dataall = readRDS("data/betarec_rho02_n500.rds")
data = dataall[[slurm_id]]
df = data.frame(y = data$y_train,
                X = data$X_train[,-1])
nburn = 1000
nsave = 5000
nthin = 1
fit <- cobin::spmicobinreg(y ~ X, df, coords = data$coords_train, NNGP = F,
                         priors = list(phi_lb = data$phi_spatial, phi_ub = data$phi_spatial,
                                       beta_intercept_scale = 100,
                                       beta_scale = 100, beta_df = Inf),
                         #nngp.control = list(n.neighbors = 30, ord = order(coords[, 1])),
                         nburn = nburn, nsave = nsave, nthin = nthin)
runtime = fit$t_mcmc
# 1. posterior summary of beta
beta_save = fit$post_save[,1:2]
u_save = fit$post_u_save

beta_mean = apply(beta_save, 2, mean)
beta_mess = mcmcse::multiESS(beta_save)
beta_summary = list(mean = beta_mean,
                    mess = beta_mess)

u_mean = apply(u_save, 2, mean)
u_mult_ess <- try(mcmcse::multiESS(u_save), silent = TRUE)
u_summary = list(mean = u_mean,
                 mult_ess = u_mult_ess)

#############################################
####### prediction ##########################
#############################################
pred_out = predict_linpred_cobit(coords_train = data$coords_train,
                                 u_save = u_save,
                                 beta_save = beta_save,
                                 sigmasq_save = fit$post_save[,"sigma.sq"],
                                 coords_test = data$coords_test,
                                 X_test = data$X_test,
                                 phi_spatial = data$phi_spatial)
# rmse
rmse_mu = sqrt(mean((colMeans(pred_out$mu_test_save) - data$mu_test)^2))

# log posterior predictive density
logpostpred = numeric(data$ntest)
for(itest in 1:data$ntest){
  logpostpred_mcmc = numeric(nsave)
  for(isave in 1:nsave){
    logpostpred_mcmc[isave] = cobin::dmicobin(data$y_test[itest], pred_out$linpred_test_save[isave,itest], fit$post_save[isave,"psi"], log = TRUE)
  }
  logpostpred[itest] = matrixStats::logSumExp(logpostpred_mcmc) - log(nsave)
}

# now combine all lists
summary_list = list(beta_summary = beta_summary,
                    u_summary = u_summary,
                    runtime = runtime,
                    rmse_mu = rmse_mu,
                    logpostpred = logpostpred,
                    avglogpostpred = mean(logpostpred))

saveRDS(summary_list, file = paste0("rho02_n500/micobinregfit_", slurm_id, ".rds"))



devtools::session_info()



