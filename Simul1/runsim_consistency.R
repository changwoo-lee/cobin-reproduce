rm(list = ls())

# --- Set working directory to current file (if running in RStudio) ---
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  try({
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  }, silent = TRUE)
}

# --- Packages ---
pkgs <- c("cobin", "betareg", "dplyr")
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



data_type <- "cobit"
nsim = 1000
ns = c(100, 400, 1600)
df_results <- data.frame()


for (n in ns) {
  if (data_type == "cobit") {
    dataall     <- readRDS(paste0("data/data_cobit_n", n, ".rds"))
    dataalltest <- readRDS(paste0("data/datatest_cobit_n", n, ".rds"))
    mu_from_lp  <- function(lp) cobin::bftprime(lp)
  } else { # "logit"
    dataall     <- readRDS(paste0("data/data_logit_n", n, ".rds"))
    dataalltest <- readRDS(paste0("data/datatest_logit_n", n, ".rds"))
    mu_from_lp  <- function(lp) 1/(1 + exp(-lp))
  }
  
  for (isim in 1:nsim) {
    df1234     <- dataall[[isim]]
    df1234test <- dataalltest[[isim]]
    
    for (datagen in c("beta","cobin","betarec","beta3mix")) {
      df     <- df1234[[datagen]]
      dftest <- df1234test[[datagen]]
      
      tmp_cobin <- safe_row({
        fit1 <- glm.cobin(y ~ X, data = df, link = "cobit", verbose = FALSE)
        fit1s <- summary(fit1)
        
        # negtestLL:
        if (data_type == "cobit") {
          eta_test <- fit1$coefficients[1] + dftest$X * fit1$coefficients[2]
        } else {
          mu_test_hat <- 1/(1 + exp(-(fit1$coefficients[1] + dftest$X * fit1$coefficients[2])))
          eta_test    <- cobin::bftprimeinv(mu_test_hat)
        }
        
        data.frame(
          simulation    = isim,
          datagen       = datagen,
          link          = "cobit",
          model         = "cobin",
          nsample       = n,
          beta1         = fit1$coefficients[2],
          beta1CIlength = fit1s$coefficients[2,2] * 2 * 1.96,
          beta1CIcover  = as.numeric(
            fit1s$coefficients[2,1] - 1.96 * fit1s$coefficients[2,2] < 1 &
              fit1s$coefficients[2,1] + 1.96 * fit1s$coefficients[2,2] > 1
          ),
          negtestLL     = -mean(dcobin(dftest$y, eta_test, fit1$lambda, log = TRUE))
        )
      }, data.frame(
        simulation = isim, datagen = datagen, link = "cobit", model = "cobin",
        nsample = n, beta1 = NA, beta1CIlength = NA, beta1CIcover = NA, negtestLL = NA
      ))
      
      tmp_beta <- safe_row({
        fit2  <- betareg(y ~ X, data = df, link = "cobit")
        fit2s <- summary(fit2)
        mu_test_hat <- mu_from_lp(fit2$coefficients$mean[1] + dftest$X * fit2$coefficients$mean[2])
        
        data.frame(
          simulation    = isim,
          datagen       = datagen,
          link          = "cobit",
          model         = "beta",
          nsample       = n,
          beta1         = fit2$coefficients$mean[2],
          beta1CIlength = fit2s$coefficients$mean[2,2] * 2 * 1.96,
          beta1CIcover  = as.numeric(
            fit2s$coefficients$mean[2,1] - 1.96 * fit2s$coefficients$mean[2,2] < 1 &
              fit2s$coefficients$mean[2,1] + 1.96 * fit2s$coefficients$mean[2,2] > 1
          ),
          negtestLL     = -mean(betareg::dbetar(dftest$y, mu_test_hat,
                                                fit2$coefficients$precision, log = TRUE))
        )
      }, data.frame(
        simulation = isim, datagen = datagen, link = "cobit", model = "beta",
        nsample = n, beta1 = NA, beta1CIlength = NA, beta1CIcover = NA, negtestLL = NA
      ))
      
      df_results <- rbind(df_results, tmp_cobin, tmp_beta)
    }
    
    if (isim %% 100 == 0) {
      cat(sprintf("coef_%s: n=%d, sim=%d\n", data_type, n, isim))
    }
  }
}

out_file <- sprintf("results_%sdata_all.rds", data_type)
saveRDS(df_results, out_file)




data_type <- "logit"
nsim = 1100
ns = c(100, 400, 1600)
df_results <- data.frame()
df_results <- data.frame()


for (n in ns) {
  if (data_type == "cobit") {
    dataall     <- readRDS(paste0("data/data_cobit_n", n, ".rds"))
    dataalltest <- readRDS(paste0("data/datatest_cobit_n", n, ".rds"))
    mu_from_lp  <- function(lp) cobin::bftprime(lp)
  } else { # "logit"
    dataall     <- readRDS(paste0("data/data_logit_n", n, ".rds"))
    dataalltest <- readRDS(paste0("data/datatest_logit_n", n, ".rds"))
    mu_from_lp  <- function(lp) 1/(1 + exp(-lp))
  }
  
  for (isim in 1:nsim) {
    df1234     <- dataall[[isim]]
    df1234test <- dataalltest[[isim]]
    
    for (datagen in c("beta","cobin","betarec","beta3mix")) {
      df     <- df1234[[datagen]]
      dftest <- df1234test[[datagen]]
      
      tmp_cobin <- safe_row({
        fit1 <- glm.cobin(y ~ X, data = df, link = "logit", verbose = FALSE)
        fit1s <- summary(fit1)
        
        # negtestLL:
        if (data_type == "cobit") {
          eta_test <- fit1$coefficients[1] + dftest$X * fit1$coefficients[2]
        } else {
          mu_test_hat <- 1/(1 + exp(-(fit1$coefficients[1] + dftest$X * fit1$coefficients[2])))
          eta_test    <- cobin::bftprimeinv(mu_test_hat)
        }
        
        data.frame(
          simulation    = isim,
          datagen       = datagen,
          link          = "logit",
          model         = "cobin",
          nsample       = n,
          beta1         = fit1$coefficients[2],
          beta1CIlength = fit1s$coefficients[2,2] * 2 * 1.96,
          beta1CIcover  = as.numeric(
            fit1s$coefficients[2,1] - 1.96 * fit1s$coefficients[2,2] < 1 &
              fit1s$coefficients[2,1] + 1.96 * fit1s$coefficients[2,2] > 1
          ),
          negtestLL     = -mean(dcobin(dftest$y, eta_test, fit1$lambda, log = TRUE))
        )
      }, data.frame(
        simulation = isim, datagen = datagen, link = "logit", model = "cobin",
        nsample = n, beta1 = NA, beta1CIlength = NA, beta1CIcover = NA, negtestLL = NA
      ))
      
      tmp_beta <- safe_row({
        fit2  <- betareg(y ~ X, data = df, link = "logit")
        fit2s <- summary(fit2)
        mu_test_hat <- mu_from_lp(fit2$coefficients$mean[1] + dftest$X * fit2$coefficients$mean[2])
        
        data.frame(
          simulation    = isim,
          datagen       = datagen,
          link          = "logit",
          model         = "beta",
          nsample       = n,
          beta1         = fit2$coefficients$mean[2],
          beta1CIlength = fit2s$coefficients$mean[2,2] * 2 * 1.96,
          beta1CIcover  = as.numeric(
            fit2s$coefficients$mean[2,1] - 1.96 * fit2s$coefficients$mean[2,2] < 1 &
              fit2s$coefficients$mean[2,1] + 1.96 * fit2s$coefficients$mean[2,2] > 1
          ),
          negtestLL     = -mean(betareg::dbetar(dftest$y, mu_test_hat,
                                                fit2$coefficients$precision, log = TRUE))
        )
      }, data.frame(
        simulation = isim, datagen = datagen, link = "logit", model = "beta",
        nsample = n, beta1 = NA, beta1CIlength = NA, beta1CIcover = NA, negtestLL = NA
      ))
      
      df_results <- rbind(df_results, tmp_cobin, tmp_beta)
    }
    
    if (isim %% 100 == 0) {
      cat(sprintf("coef_%s: n=%d, sim=%d\n", data_type, n, isim))
    }
  }
}

out_file <- sprintf("results_%sdata_all.rds", data_type)
saveRDS(df_results, out_file)



summarize_coef_bias_rmse <- function(file, keep_first_1000 = TRUE) {
  df <- readRDS(file)
  
  # Drop failed sims and clean negtestLL
  nasim <- df$simulation[is.na(df$beta1)]
  df$negtestLL[is.infinite(df$negtestLL)] <- NA
  df <- dplyr::filter(df, !simulation %in% nasim)
  
  # Optionally keep only first <= 1000 unique simulation IDs
  if (keep_first_1000) {
    K <- min(1000, length(unique(df$simulation)))
    df <- dplyr::filter(df, simulation <= sort(unique(simulation))[K])
  }
  
  out <- df %>%
    dplyr::group_by(datagen, model, nsample) %>%
    dplyr::summarise(
      bias = round(mean(beta1) - 1, 3),
      rmse = round(sqrt(mean((beta1 - 1)^2)), 3),
      .groups = "drop"
    ) %>% as.data.frame()
  
  print(out)
  invisible(out)
}


summarize_coef_ci_cover_ll <- function(file, keep_first_1000 = TRUE) {
  df <- readRDS(file)
  
  # Drop failed sims and clean negtestLL
  nasim <- df$simulation[is.na(df$beta1)]
  df$negtestLL[is.infinite(df$negtestLL)] <- NA
  df <- dplyr::filter(df, !simulation %in% nasim)
  
  # Optionally keep only first <= 1000 unique simulation IDs
  if (keep_first_1000) {
    K <- min(1000, length(unique(df$simulation)))
    df <- dplyr::filter(df, simulation <= sort(unique(simulation))[K])
  }
  out <- df %>%
    dplyr::group_by(datagen, model, nsample) %>%
    dplyr::summarise(
      nsim = n(),
      cilength  = round(mean(beta1CIlength), 3),
      cilength_se  =round(sd(beta1CIlength)/sqrt(n()), 3),
      cover     = round(mean(beta1CIcover), 3),
      NTLL = round(mean(negtestLL, na.rm = TRUE),3),
      NTLL_se = round(sd(negtestLL, na.rm = T)/sqrt(n()),3),
      .groups   = "drop"
    )  %>% as.data.frame()
  
  print(out)
  invisible(out)
}

summarize_coef_bias_rmse("results_cobitdata_all.rds",  keep_first_1000 = 1000)
summarize_coef_bias_rmse("results_logitdata_all.rds",  keep_first_1000 = 1000)

# Monte Carlo standard errors of bias and rmse can be separately calculated based on
# Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. Statistics in medicine, 38(11), 2074-2102.

summarize_coef_ci_cover_ll("results_cobitdata_all.rds",  keep_first_1000 = 1000)
summarize_coef_ci_cover_ll("results_logitdata_all.rds",  keep_first_1000 = 1000)

# Monte Carlo standard errors of cover can be separately calculated based on
# Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. Statistics in medicine, 38(11), 2074-2102.
