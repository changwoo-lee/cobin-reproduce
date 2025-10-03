rm(list = ls())

# --- Set working directory to current file (if running in RStudio) ---
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  try({
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  }, silent = TRUE)
}

# --- Packages ---
pkgs <- c("cobin", "betareg", "gamlss", "gamlss.dist", "dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Please install it.", p))
  }
}
library(cobin)
library(betareg)
library(gamlss)
library(gamlss.dist)
library(dplyr)

# Make sure this file exists alongside this script
if (file.exists("betareg_cobit.R")) {
  source("betareg_cobit.R") # betareg with cobit link (requires cobin loaded)
}

safe_row <- function(expr, fallback) {
  tryCatch(expr, error = function(e) { message("Error: ", e$message); fallback })
}
# ============================================================
# t3 density under transformed response models
# ============================================================

# PDF of Y = bftprime(X) with X ~ t3(mu, sigma)
dt3_cobit <- function(y, mu = 0, sigma = 1, log = FALSE) {
  stopifnot(is.numeric(y), is.numeric(mu), is.numeric(sigma), sigma > 0)
  xinvs <- cobin::bftprimeinv(y)
  f_t3 <- function(x) {
    2 / (pi * sqrt(3) * sigma) * (1 + ((x - mu)^2) / (3 * sigma^2))^-2
  }
  dens <- f_t3(xinvs) / cobin::bftprimeprime(xinvs)
  bad <- !(y > 0 & y < 1) | !is.finite(xinvs)
  dens[bad] <- 0
  if (log) log(dens) else dens
}

# helpers
logistic_fun <- function(x) 1 / (1 + exp(-x))
logit_fun <- function(yy) log(yy) - log1p(-yy)

# PDF of Y = logistic(X) with X ~ t3(mu, sigma)
dt3_logit <- function(y, mu = 0, sigma = 1, log = FALSE) {
  stopifnot(is.numeric(y), is.numeric(mu), is.numeric(sigma), sigma > 0)
  xinvs <- logit_fun(y)
  f_t3 <- function(x) {
    2 / (pi * sqrt(3) * sigma) * (1 + ((x - mu)^2) / (3 * sigma^2))^-2
  }
  # f_Y(y) = f_X(x) / g'(x) at x = g^{-1}(y); g'(x) = logistic(x)*(1 - logistic(x))
  gprime <- logistic_fun(xinvs) * (1 - logistic_fun(xinvs))
  dens <- f_t3(xinvs) / gprime
  bad <- !(y > 0 & y < 1) | !is.finite(xinvs) | !is.finite(dens)
  dens[bad] <- 0
  if (log) log(dens) else dens
}



data_type <- "cobit"
nsim = 1000
ns = c(100, 400, 1600)
df_results <- data.frame()
# WARNING: due to boundary-proximate data, beta regression fits may fail.
# We retain only 500 successful simulations out of 1100 
for (n in ns) {
  if (data_type == "logit") {
    dataall     <- readRDS(paste0("data/data_logit_n", n, ".rds"))
    dataalltest <- readRDS(paste0("data/datatest_logit_n", n, ".rds"))
    # true mean function on test set for "logit-data"
    true_mu <- function(x) 1/(1+exp(-x))
    cross_link <- "cobit"  
  } else {
    dataall     <- readRDS(paste0("data/data_cobit_n", n, ".rds"))
    dataalltest <- readRDS(paste0("data/datatest_cobit_n", n, ".rds"))
    # true mean function on test set for "cobit-data"
    true_mu <- function(x) cobin::bftprime(x)
    cross_link <- "logit"  
  }
  
  for (isim in 1:nsim) {
    df1234     <- dataall[[isim]]
    df1234test <- dataalltest[[isim]]
    
    for (datagen in c("beta","cobin","betarec","beta3mix")) {
      df     <- df1234[[datagen]]
      dftest <- df1234test[[datagen]]
      
      if (data_type == "logit") {
        tmp1 <- safe_row({
          fit1 <- glm.cobin(y ~ X, data = df, link = "cobit", verbose = FALSE)
          mu_true <- true_mu(dftest$X)
          mu_est  <- cobin::bftprime(fit1$coefficients[1] + dftest$X*fit1$coefficients[2])
          data.frame(simulation = isim, datagen = datagen, datagenlink = "logit",
                     link = "cobit", model = "cobin", nsample = n,
                     mspe = mean((mu_est - mu_true)^2),
                     negtestLL = -mean(dcobin(dftest$y,
                                              fit1$coefficients[1] + dftest$X*fit1$coefficients[2],
                                              fit1$lambda, log = TRUE)))
        }, data.frame(simulation = isim, datagen = datagen, datagenlink = "logit",
                      link = "cobit", model = "cobin", nsample = n, mspe = NA, negtestLL = NA))
        
        tmp2 <- safe_row({
          fit2 <- betareg(y ~ X, data = df, link = "cobit")
          mu_true <- true_mu(dftest$X)
          mu_est  <- cobin::bftprime(fit2$coefficients$mean[1] + dftest$X*fit2$coefficients$mean[2])
          data.frame(simulation = isim, datagen = datagen, datagenlink = "logit",
                     link = "cobit", model = "beta", nsample = n,
                     mspe = mean((mu_est - mu_true)^2),
                     negtestLL = -mean(betareg::dbetar(dftest$y, mu_est,
                                                       fit2$coefficients$precision, log = TRUE)))
        }, data.frame(simulation = isim, datagen = datagen, datagenlink = "logit",
                      link = "cobit", model = "beta", nsample = n, mspe = NA, negtestLL = NA))
        
      } else { 
        tmp1 <- safe_row({
          fit1 <- glm.cobin(y ~ X, data = df, link = "logit", verbose = FALSE)
          mu_true <- true_mu(dftest$X)
          mu_est  <- logistic_fun(fit1$coefficients[1] + dftest$X*fit1$coefficients[2])
          data.frame(simulation = isim, datagen = datagen, datagenlink = "cobit",
                     link = "logit", model = "cobin", nsample = n,
                     mspe = mean((mu_est - mu_true)^2),
                     negtestLL = -mean(dcobin(dftest$y,
                                              cobin::bftprimeinv(mu_est),
                                              fit1$lambda, log = TRUE)))
        }, data.frame(simulation = isim, datagen = datagen, datagenlink = "cobit",
                      link = "logit", model = "cobin", nsample = n, mspe = NA, negtestLL = NA))
        
        tmp2 <- safe_row({
          fit2 <- betareg(y ~ X, data = df, link = "logit")
          mu_true <- true_mu(dftest$X)
          mu_est  <- logistic_fun(fit2$coefficients$mean[1] + dftest$X*fit2$coefficients$mean[2])
          data.frame(simulation = isim, datagen = datagen, datagenlink = "cobit",
                     link = "logit", model = "beta", nsample = n,
                     mspe = mean((mu_est - mu_true)^2),
                     negtestLL = -mean(betareg::dbetar(dftest$y, mu_est,
                                                       fit2$coefficients$precision, log = TRUE)))
        }, data.frame(simulation = isim, datagen = datagen, datagenlink = "cobit",
                      link = "logit", model = "beta", nsample = n, mspe = NA, negtestLL = NA))
      }
      
      if(data_type == "cobit"){
      tmp3 <- safe_row({
        df$ylogit <- log(df$y) - log(1 - df$y)
        fit_logit <- gamlss(ylogit ~ X, data = df,
                            sigma.formula = ~1,
                            nu.formula = ~1, nu.start = 3, nu.fix = TRUE,
                            family = TF(), trace= F)
        mu_true <- true_mu(dftest$X)
        mu_est  <- 1/(1+exp(-fit_logit$mu.coefficients[1] - dftest$X*fit_logit$mu.coefficients[2]))
        data.frame(simulation = isim, datagen = datagen,
                   datagenlink = data_type, link = cross_link,
                   model = "t3logit", nsample = n,
                   mspe = mean((mu_est - mu_true)^2),
                   negtestLL = -mean(dt3_logit(dftest$y,
                                               fit_logit$mu.coefficients[1] + dftest$X*fit_logit$mu.coefficients[2],
                                               exp(fit_logit$sigma.coefficients),
                                               log = TRUE)))
      }, data.frame(simulation = isim, datagen = datagen,
                    datagenlink = data_type, link = cross_link,
                    model = "t3logit", nsample = n, mspe = NA, negtestLL = NA))
      }else{
      tmp3 <- safe_row({
        df$ycobit <- cobin::bftprimeinv(df$y)
        fit_cobit <- gamlss(ycobit ~ X, data = df,
                            sigma.formula = ~1,
                            nu.formula = ~1, nu.start = 3, nu.fix = TRUE,
                            family = TF(), trace= F)
        mu_true <- true_mu(dftest$X)
        mu_est  <- cobin::bftprime(fit_cobit$mu.coefficients[1] + dftest$X*fit_cobit$mu.coefficients[2])
        data.frame(simulation = isim, datagen = datagen,
                   datagenlink = data_type, link = cross_link,
                   model = "t3cobit", nsample = n,
                   mspe = mean((mu_est - mu_true)^2),
                   negtestLL = -mean(dt3_cobit(dftest$y,
                                               fit_cobit$mu.coefficients[1] + dftest$X*fit_cobit$mu.coefficients[2],
                                               exp(fit_cobit$sigma.coefficients),
                                               log = TRUE)))
      }, data.frame(simulation = isim, datagen = datagen,
                    datagenlink = data_type, link = cross_link,
                    model = "t3cobit", nsample = n, mspe = NA, negtestLL = NA))
      }
      # accumulate
      df_results <- rbind(df_results, tmp1, tmp2, tmp3)
    }
    if (isim %% 100 == 0) cat(sprintf("%sdata_all: n=%d, sim=%d\n", data_type, n, isim))
  }
}

out_file <- sprintf("df_results_%sdata_wrongmean_all.rds", data_type)
saveRDS(df_results, out_file)



data_type <- "logit"
nsim = 1100
ns = c(100, 400, 1600)
df_results <- data.frame()

# WARNING: due to boundary-proximate data, many beta regression fits may fail.
# We retain only 500 successful simulations out of 1100 
for (n in ns) {
  if (data_type == "logit") {
    dataall     <- readRDS(paste0("data/data_logit_n", n, ".rds"))
    dataalltest <- readRDS(paste0("data/datatest_logit_n", n, ".rds"))
    # true mean function on test set for "logit-data"
    true_mu <- function(x) 1/(1+exp(-x))
    cross_link <- "cobit"  
  } else {
    dataall     <- readRDS(paste0("data/data_cobit_n", n, ".rds"))
    dataalltest <- readRDS(paste0("data/datatest_cobit_n", n, ".rds"))
    # true mean function on test set for "cobit-data"
    true_mu <- function(x) cobin::bftprime(x)
    cross_link <- "logit"  
  }
  
  for (isim in 1:nsim) {
    df1234     <- dataall[[isim]]
    df1234test <- dataalltest[[isim]]
    
    for (datagen in c("beta","cobin","betarec","beta3mix")) {
      df     <- df1234[[datagen]]
      dftest <- df1234test[[datagen]]
      
      if (data_type == "logit") {
        tmp1 <- safe_row({
          fit1 <- glm.cobin(y ~ X, data = df, link = "cobit", verbose = FALSE)
          mu_true <- true_mu(dftest$X)
          mu_est  <- cobin::bftprime(fit1$coefficients[1] + dftest$X*fit1$coefficients[2])
          data.frame(simulation = isim, datagen = datagen, datagenlink = "logit",
                     link = "cobit", model = "cobin", nsample = n,
                     mspe = mean((mu_est - mu_true)^2),
                     negtestLL = -mean(dcobin(dftest$y,
                                              fit1$coefficients[1] + dftest$X*fit1$coefficients[2],
                                              fit1$lambda, log = TRUE)))
        }, data.frame(simulation = isim, datagen = datagen, datagenlink = "logit",
                      link = "cobit", model = "cobin", nsample = n, mspe = NA, negtestLL = NA))
        
        tmp2 <- safe_row({
          fit2 <- betareg(y ~ X, data = df, link = "cobit")
          mu_true <- true_mu(dftest$X)
          mu_est  <- cobin::bftprime(fit2$coefficients$mean[1] + dftest$X*fit2$coefficients$mean[2])
          data.frame(simulation = isim, datagen = datagen, datagenlink = "logit",
                     link = "cobit", model = "beta", nsample = n,
                     mspe = mean((mu_est - mu_true)^2),
                     negtestLL = -mean(betareg::dbetar(dftest$y, mu_est,
                                                       fit2$coefficients$precision, log = TRUE)))
        }, data.frame(simulation = isim, datagen = datagen, datagenlink = "logit",
                      link = "cobit", model = "beta", nsample = n, mspe = NA, negtestLL = NA))
        
      } else { 
         tmp1 <- safe_row({
          fit1 <- glm.cobin(y ~ X, data = df, link = "logit", verbose = FALSE)
          mu_true <- true_mu(dftest$X)
          mu_est  <- logistic_fun(fit1$coefficients[1] + dftest$X*fit1$coefficients[2])
          data.frame(simulation = isim, datagen = datagen, datagenlink = "cobit",
                     link = "logit", model = "cobin", nsample = n,
                     mspe = mean((mu_est - mu_true)^2),
                     negtestLL = -mean(dcobin(dftest$y,
                                              cobin::bftprimeinv(mu_est),
                                              fit1$lambda, log = TRUE)))
        }, data.frame(simulation = isim, datagen = datagen, datagenlink = "cobit",
                      link = "logit", model = "cobin", nsample = n, mspe = NA, negtestLL = NA))
        
        tmp2 <- safe_row({
          fit2 <- betareg(y ~ X, data = df, link = "logit")
          mu_true <- true_mu(dftest$X)
          mu_est  <- logistic_fun(fit2$coefficients$mean[1] + dftest$X*fit2$coefficients$mean[2])
          data.frame(simulation = isim, datagen = datagen, datagenlink = "cobit",
                     link = "logit", model = "beta", nsample = n,
                     mspe = mean((mu_est - mu_true)^2),
                     negtestLL = -mean(betareg::dbetar(dftest$y, mu_est,
                                                       fit2$coefficients$precision, log = TRUE)))
        }, data.frame(simulation = isim, datagen = datagen, datagenlink = "cobit",
                      link = "logit", model = "beta", nsample = n, mspe = NA, negtestLL = NA))
      }
      
      if(data_type == "cobit"){
        tmp3 <- safe_row({
          df$ylogit <- log(df$y) - log(1 - df$y)
          fit_logit <- gamlss(ylogit ~ X, data = df,
                              sigma.formula = ~1,
                              nu.formula = ~1, nu.start = 3, nu.fix = TRUE,
                              family = TF(), trace= F)
          mu_true <- true_mu(dftest$X)
          mu_est  <- 1/(1+exp(-fit_logit$mu.coefficients[1] - dftest$X*fit_logit$mu.coefficients[2]))
          data.frame(simulation = isim, datagen = datagen,
                     datagenlink = data_type, link = cross_link,
                     model = "t3logit", nsample = n,
                     mspe = mean((mu_est - mu_true)^2),
                     negtestLL = -mean(dt3_logit(dftest$y,
                                                 fit_logit$mu.coefficients[1] + dftest$X*fit_logit$mu.coefficients[2],
                                                 exp(fit_logit$sigma.coefficients),
                                                 log = TRUE)))
        }, data.frame(simulation = isim, datagen = datagen,
                      datagenlink = data_type, link = cross_link,
                      model = "t3logit", nsample = n, mspe = NA, negtestLL = NA))
      }else{
        tmp3 <- safe_row({
          df$ycobit <- cobin::bftprimeinv(df$y)
          fit_cobit <- gamlss(ycobit ~ X, data = df,
                              sigma.formula = ~1,
                              nu.formula = ~1, nu.start = 3, nu.fix = TRUE,
                              family = TF(), trace= F)
          mu_true <- true_mu(dftest$X)
          mu_est  <- cobin::bftprime(fit_cobit$mu.coefficients[1] + dftest$X*fit_cobit$mu.coefficients[2])
          data.frame(simulation = isim, datagen = datagen,
                     datagenlink = data_type, link = cross_link,
                     model = "t3cobit", nsample = n,
                     mspe = mean((mu_est - mu_true)^2),
                     negtestLL = -mean(dt3_cobit(dftest$y,
                                                 fit_cobit$mu.coefficients[1] + dftest$X*fit_cobit$mu.coefficients[2],
                                                 exp(fit_cobit$sigma.coefficients),
                                                 log = TRUE)))
        }, data.frame(simulation = isim, datagen = datagen,
                      datagenlink = data_type, link = cross_link,
                      model = "t3cobit", nsample = n, mspe = NA, negtestLL = NA))
      }
      # accumulate
      df_results <- rbind(df_results, tmp1, tmp2, tmp3)
    }
    if (isim %% 100 == 0) cat(sprintf("%sdata_all: n=%d, sim=%d\n", data_type, n, isim))
  }
}

out_file <- sprintf("df_results_%sdata_wrongmean_all.rds", data_type)
saveRDS(df_results, out_file)



summarize_results <- function(file, keep_first_500 = TRUE) {
  df_results <- readRDS(file)
  # Clean NA/Inf and restrict to first 500 simulations if requested
  nasim <- df_results$simulation[which(is.na(df_results$mspe))]
  df_results$negtestLL[is.infinite(df_results$negtestLL)] <- NA
  df_results <- df_results %>% filter(!simulation %in% nasim)
  if (keep_first_500) {
    df_results <- df_results %>% filter(simulation <= sort(unique(simulation))[min(500, length(unique(simulation)))])
  }
  
  cat("\n== Means ==\n")
  print(df_results %>%
          group_by(datagen, model, nsample) %>%
          summarise(
            mspe = round(mean(mspe)*100, 3),
            negtestLL = round(mean(negtestLL, na.rm = TRUE), 3),
            .groups = "drop"
          ) %>% as.data.frame())
  
  cat("\n== SEs (over sims) ==\n")
  print(df_results %>%
          group_by(datagen, model, nsample) %>%
          summarise(
            mspe = round(sd(mspe)*100/sqrt(n()), 3),
            negtestLL = round(sd(negtestLL, na.rm = TRUE)/sqrt(n()), 3),
            .groups = "drop"
          ) %>% as.data.frame())
}

summarize_results("df_results_cobitdata_wrongmean_all.rds")
summarize_results("df_results_logitdata_wrongmean_all.rds")
