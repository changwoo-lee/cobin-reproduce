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

# 3-point beta mixture
rbeta3 = function(n, mu, phi){
  if(length(mu)==1) mu = rep(mu,n)
  y = numeric(n)
  mu_left = mu - pmin(mu, 1-mu)/2
  mu_right = mu + pmin(mu, 1-mu)/2
  for(i in 1:n){
    u = runif(1)
    if(u < 0.25){
      y[i] = betareg::rbetar(1, mu_left[i], phi )
    }else if(u < 0.5){
      y[i] = betareg::rbetar(1, mu_right[i], phi )
    }else{
      y[i] = betareg::rbetar(1, mu[i], phi)
    }
  }
  return(y)
}

#SAVE = FALSE
SAVE = TRUE


#########################
##### cobit link ########
#########################
xsd = 3

n = 100
dataall = list()
for(isim in 1:nsim){
  set.seed(isim+n*1000)
  
  X = matrix(rnorm(n, 0, sd = xsd), nrow = n, 1)
  eta = as.numeric(X) # beta0 = 0, beta1 = 1
  mu = as.numeric(cobin::bftprime(eta)) # cobit link
  
  # beta
  y = betareg::rbetar(n, mu, phi = rep(8,n))
  df1 = data.frame(y = y, X= X)
  
  # cobin
  y = cobin::rcobin(n, eta, rep(3,n))
  df2 = data.frame(y = y, X= X)
  
  # beta rectangular
  y = rbetarec(n, mu, alpha = 0.2, phi = 10)
  df3 = data.frame(y = y, X= X)
  
  # beta mixture
  y = rbeta3(n, mu, phi = 40)
  df4 = data.frame(y = y, X= X)
  
  dataall[[isim]] = list("beta" = df1, "cobin" = df2, "betarec" = df3, "beta3mix" = df4)
}
if(SAVE) saveRDS(dataall, file = paste0("data_cobit_n",n,".rds"))



n = 400
dataall = list()
for(isim in 1:nsim){
  set.seed(isim+n*1000)
  
  X = matrix(rnorm(n, 0, sd = xsd), nrow = n, 1)
  eta = as.numeric(X) # beta0 = 0, beta1 = 1
  mu = as.numeric(cobin::bftprime(eta)) # cobit link
  
  # beta
  y = betareg::rbetar(n, mu, phi = rep(8,n))
  df1 = data.frame(y = y, X= X)
  
  # cobin
  y = cobin::rcobin(n, eta, rep(3,n))
  df2 = data.frame(y = y, X= X)
  
  # beta rectangular
  y = rbetarec(n, mu, alpha = 0.2, phi = 10)
  df3 = data.frame(y = y, X= X)
  
  # beta mixture
  y = rbeta3(n, mu, phi = 40)
  df4 = data.frame(y = y, X= X)
  
  dataall[[isim]] = list("beta" = df1, "cobin" = df2, "betarec" = df3, "beta3mix" = df4)
}
if(SAVE) saveRDS(dataall, file = paste0("data_cobit_n",n,".rds"))



n = 1600
dataall = list()
for(isim in 1:nsim){
  set.seed(isim+n*1000)
  
  X = matrix(rnorm(n, 0, sd = xsd), nrow = n, 1)
  eta = as.numeric(X) # beta0 = 0, beta1 = 1
  mu = as.numeric(cobin::bftprime(eta)) # cobit link
  
  # beta
  y = betareg::rbetar(n, mu, phi = rep(8,n))
  df1 = data.frame(y = y, X= X)
  
  # cobin
  y = cobin::rcobin(n, eta, rep(3,n))
  df2 = data.frame(y = y, X= X)
  
  # beta rectangular
  y = rbetarec(n, mu, alpha = 0.2, phi = 10)
  df3 = data.frame(y = y, X= X)
  
  # beta mixture
  y = rbeta3(n, mu, phi = 40)
  df4 = data.frame(y = y, X= X)
  
  dataall[[isim]] = list("beta" = df1, "cobin" = df2, "betarec" = df3, "beta3mix" = df4)
}
if(SAVE) saveRDS(dataall, file = paste0("data_cobit_n",n,".rds"))



#########################
##### logit link ########
#########################

nsim = 1100
#nsim = 1000 
beta1 = 1
# logit link
xsd = 1

n = 100
dataall = list()
for(isim in 1:nsim){
  set.seed(isim+n*10000)

  X = matrix(rnorm(n, 0, sd = xsd), nrow = n, 1)
  eta = as.numeric(X) # beta0 = 0, beta1 = 1
  mu = 1/(1+exp(-eta)) # logit
  # beta
  y = betareg::rbetar(n, mu, phi = rep(8,n))
  df1 = data.frame(y = y, X= X)

  # cobin
  y = cobin::rcobin(n, bftprimeinv(mu), rep(3,n))
  df2 = data.frame(y = y, X= X)

  # beta rectangular
  y = rbetarec(n, mu, alpha = 0.2, phi = 10)
  df3 = data.frame(y = y, X= X)

  # beta mixture
  y = rbeta3(n, mu, phi = 40)
  df4 = data.frame(y = y, X= X)

  dataall[[isim]] = list("beta" = df1, "cobin" = df2, "betarec" = df3, "beta3mix" = df4)
}
if(SAVE) saveRDS(dataall, file = paste0("data_logit_n",n,".rds"))



n = 400
dataall = list()
for(isim in 1:nsim){
  set.seed(isim+n*10000)

  X = matrix(rnorm(n, 0, sd = xsd), nrow = n, 1)
  eta = as.numeric(X) # beta0 = 0, beta1 = 1
  mu = 1/(1+exp(-eta)) # logit
  # beta
  y = betareg::rbetar(n, mu, phi = rep(8,n))
  df1 = data.frame(y = y, X= X)

  # cobin
  y = cobin::rcobin(n, bftprimeinv(mu), rep(3,n))
  df2 = data.frame(y = y, X= X)

  # beta rectangular
  y = rbetarec(n, mu, alpha = 0.2, phi = 10)
  df3 = data.frame(y = y, X= X)

  # beta mixture
  y = rbeta3(n, mu, phi = 40)
  df4 = data.frame(y = y, X= X)

  dataall[[isim]] = list("beta" = df1, "cobin" = df2, "betarec" = df3, "beta3mix" = df4)
}
if(SAVE) saveRDS(dataall, file = paste0("data_logit_n",n,".rds"))



n = 1600
dataall = list()
for(isim in 1:nsim){
  set.seed(isim+n*10000)

  X = matrix(rnorm(n, 0, sd = xsd), nrow = n, 1)
  eta = as.numeric(X) # beta0 = 0, beta1 = 1
  mu = 1/(1+exp(-eta)) # logit
  # beta
  y = betareg::rbetar(n, mu, phi = rep(8,n))
  df1 = data.frame(y = y, X= X)

  # cobin
  y = cobin::rcobin(n, bftprimeinv(mu), rep(3,n))
  df2 = data.frame(y = y, X= X)

  # beta rectangular
  y = rbetarec(n, mu, alpha = 0.2, phi = 10)
  df3 = data.frame(y = y, X= X)

  # beta mixture
  y = rbeta3(n, mu, phi = 40)
  df4 = data.frame(y = y, X= X)

  dataall[[isim]] = list("beta" = df1, "cobin" = df2, "betarec" = df3, "beta3mix" = df4)
}
if(SAVE) saveRDS(dataall, file = paste0("data_logit_n",n,".rds"))





