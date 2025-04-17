
rm(list = ls())
slurm <- TRUE
if(slurm){
  slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
} else{
  slurm_id <- 0

}

# Load libraries
library(brms)
library(rstan)
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
stan_data = standata(brm(y ~ X,
                         data = df,
                         family = Beta(),
                         chains = 1,
                         iter = nburn + nsave,
                         warmup = nburn, thin = nthin,
                         prior = prior("normal(0, 100)", class = "b"), empty = TRUE))
stan_data$coords = data$coords_train
stan_data$phi_spatial = data$phi_spatial

out = rstan::stan(file= "betareg_spatial.stan",
                  data = stan_data,
                  chains = 1,
                  iter = nburn + nsave, warmup = nburn, thin = nthin)

#pars <- rstan::extract(out)
pars = rstan::extract(out, permute = FALSE)
# 1. posterior summary of beta
#beta_save = cbind(pars$b_Intercept, pars$b)
beta_save = pars[,1,c("b_Intercept","b[1]")]
#u_save = pars$u
u_save = pars[,1,grep("^u\\[\\d+\\]$", dimnames(pars)$parameters)]
sigmasq_save = pars[,1,"sigmasq"]
phi_save = pars[,1,"phi"]

beta_mean = apply(beta_save, 2, mean)
beta_mess = mcmcse::multiESS(beta_save)
beta_summary = list(mean = beta_mean,
                    mess = beta_mess)

u_mean = apply(u_save, 2, mean)
u_mult_ess <- try(mcmcse::multiESS(u_save), silent = TRUE)
u_summary = list(mean = u_mean,
                 mult_ess = u_mult_ess)

runtime = as.numeric(sum(rstan::get_elapsed_time(out))) # in seconds

#############################################
####### prediction ##########################
#############################################
pred_out = predict_linpred_cobit(coords_train = data$coords_train,
                                 u_save = u_save,
                                 beta_save = beta_save,
                                 sigmasq_save = sigmasq_save,
                                 coords_test = data$coords_test,
                                 X_test = data$X_test,
                                 phi_spatial = data$phi_spatial)
# rmse
rmse_mu = sqrt(mean((colMeans(pred_out$mu_test_save) - data$mu_test)^2))

# log posterior predictive density
logpostpred = numeric(data$ntest)
for(itest in 1:data$ntest){
  logpostpred[itest] = matrixStats::logSumExp(betareg::dbetar(data$y_test[itest], as.numeric(pred_out$mu_test_save[,itest]), as.numeric(phi_save), log = TRUE)) - log(nsave)
}


# now combine all lists
summary_list = list(beta_summary = beta_summary,
                    u_summary = u_summary,
                    runtime = runtime,
                    rmse_mu = rmse_mu,
                    logpostpred = logpostpred,
                    avglogpostpred = mean(logpostpred))

saveRDS(summary_list, file = paste0("rho02_n500/betaregfit_", slurm_id, ".rds"))


devtools::session_info()



