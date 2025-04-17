rm(list = ls())
# set path as current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(cobin)
library(betareg)


# beta rectangular distribution
# see Bayes, C. L., Bazán, J. L., & Garcıa, C. (2012). A New Robust Regression Model for Proportions. Bayesian Analysis, 7(4), 841-866.
rbetarec = function(n, mu, alpha, phi){
  if(length(mu)==1) mu = rep(mu,n)
  w = (1-alpha*(1-abs(2*mu-1)))
  y = numeric(n)
  for(i in 1:n){
    if(runif(1) < w[i]){
      y[i] = betareg::rbetar(1, (mu[i] - 0.5*alpha*(1-abs(2*mu[i]-1)))/(1-alpha*(1-abs(2*mu[i]-1))), phi )
    }else{
      y[i] = runif(1)
    }
  }
  return(y)
}

# test
# mu = 0.3
# hist(rbetarec(10000, mu = mu, alpha = 0.2, phi = 10), breaks = 500)
#
#


simulate_betarec <- function(phi_spatial = 10, # spatial dependence, Matern 1.5, 0.05 or 0.1
                             sigmasq_spatial = 1, # variance of the random effects
                             ntrain, # number of training data, 200
                             ntest, # number of testing point, 50
                             seed){
  set.seed(seed)
  # Sample the coordinates at random
  coords_train <- cbind(runif(ntrain, 0, 1), runif(ntrain, 0, 1)) # spatial locations
  coords_test <- cbind(runif(ntest, 0, 1), runif(ntest, 0, 1))
  # Pre-calculate the distances for the Matern covariances
  D <- fields::rdist(rbind(coords_train, coords_test)) # (n + ntest) x (n + ntest)
  # Daa <- D[1:ntrain, 1:ntrain] # n x n
  # Dab <- D[1:ntrain, (ntrain+1):(ntrain + ntest)] # n x n
  # Dbb <- D[(ntrain+1):(ntrain + ntest), (ntrain+1):(ntrain + ntest)] # ntest x ntest

  # fixed effects
  Xtrain = cbind(rep(1, ntrain), rnorm(ntrain, 0, 3))
  Xtest = cbind(rep(1, ntest), rnorm(ntest, 0, 3))

  # random effects
  R = fields::Matern(D, smoothness = 0.5, range = 1/phi_spatial)
  u = as.numeric(mvnfast::rmvn(1, mu = rep(0,ncol(R)), sigmasq_spatial*R) )
  #u = rbridge(ntrain + ntest, phi)
  utrain = u[1:ntrain]
  utest = u[(ntrain+1):(ntrain+ntest)]

  # true beta
  beta = c(0, 1)

  linpred_train = as.numeric(Xtrain%*%beta + utrain)
  linpred_test = as.numeric(Xtest%*%beta + utest)

  mu_train = cobin::bftprime(linpred_train)
  mu_test = cobin::bftprime(linpred_test)

  y_train = rbetarec(ntrain, mu_train, alpha = 0.2, phi = 10)
  y_test = rbetarec(ntest, mu_test, alpha = 0.2, phi = 10)

  # Save the data
  data <- list(y_train = y_train,
               y_test = y_test,
               X_train = Xtrain,
               X_test = Xtest,
               mu_train = mu_train,
               mu_test = mu_test,
               u_train = utrain,
               u_test = utest,
               beta = beta,
               coords_train = coords_train,
               coords_test = coords_test,
               phi_spatial = phi_spatial,
               sigmasq_spatial = sigmasq_spatial,
               beta= beta,
               ntrain = ntrain,
               ntest = ntest,
               seed = seed,
               model = "betarec")
  return(data)
}


nsim = 100

sigmasq_spatial = 1

out1 = list()
for(isim in 1:nsim){
  out1[[isim]] = simulate_betarec(phi_spatial = 5, # rho = 0.2
                                  sigmasq_spatial = sigmasq_spatial, # variance of the random effects
                                  ntrain = 200, # number of training data
                                  ntest = 50, # number of testing point
                                  seed = isim + 5000)
}
saveRDS(out1, "data/betarec_rho02_n250.rds")

out11 = list()
for(isim in 1:nsim){
  out11[[isim]] = simulate_betarec(phi_spatial = 5, # # rho = 0.2
                                  sigmasq_spatial = sigmasq_spatial, # variance of the random effects
                                  ntrain = 400, # number of training data
                                  ntest = 100, # number of testing point
                                  seed = isim + 50000)
}
saveRDS(out11, "data/betarec_rho02_n500.rds")


out2 = list()
for(isim in 1:nsim){
  out2[[isim]] = simulate_betarec(phi_spatial = 10, # rho = 0.1
                                  sigmasq_spatial = sigmasq_spatial, # variance of the random effects
                                  ntrain = 200, # number of training data
                                  ntest = 50, # number of testing point
                                  seed = isim + 10000)
}
saveRDS(out2, "data/betarec_rho01_n250.rds")



out22 = list()
for(isim in 1:nsim){
  out22[[isim]] = simulate_betarec(phi_spatial = 10, # rho = 0.1
                                  sigmasq_spatial = sigmasq_spatial, # variance of the random effects
                                  ntrain = 400, # number of training data, 200
                                  ntest = 100, # number of testing point, 50
                                  seed = isim + 100000)
}
saveRDS(out22, "data/betarec_rho01_n500.rds")

