
rm(list = ls())

# set path as current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(cobin)
library(betareg)
source("betareg_cobit.R") # betareg with cobit link. only works when cobin library is loaded


df_results <- data.frame()
# A. cobit link
#nsim = 1100
nsim = 1000
for(n in c(100, 400, 1600)){
  dataall = readRDS(paste0("data_cobit_n",n,".rds"))
  for(isim in 1:nsim){
    df1234 = dataall[[isim]]

    for(datagen in c("beta","cobin","betarec","beta3mix")){
      df = df1234[[datagen]]
      tryCatch({
        fit1 = glm.cobin(y ~ X, data = df, link = "cobit", verbose = F)
      }, error = function(e) {
        message("Error in glm.cobin in iteration ", isim, " and n =", n,": ", e$message)
        fit1 <<- list(coefficients = c(NA, NA))
        return(NULL)
      })
      df_temp1 <- data.frame("simulation" = isim,
                             "datagen" = datagen,
                             "link" = "cobit",
                             "model" = "cobin",
                             "nsample" = n,
                             "beta1" = fit1$coefficients[2])
      tryCatch({
        fit2 = betareg(y ~ X, data = df, link = "cobit")
      }, error = function(e) {
        message("Error in betareg in iteration ", isim, " and n =", n,": ", e$message)
        fit2 <<- list(coefficients = list(mean = c(NA, NA)))
        return(NULL)
      })
      df_temp2 <- data.frame("simulation" = isim,
                             "datagen" = datagen,
                             "link" = "cobit",
                             "model" = "beta",
                             "nsample" = n,
                             "beta1" = fit2$coefficients$mean[2])

      df_results <- rbind(df_results, df_temp1, df_temp2)
    }
    if(isim %% 100 == 0){
      # print n and isim
      print(paste0("n = ",n, "isim = ",isim))
    }
  }
}

#saveRDS(df_results,"df_results_cobit.rds")
#df_results = readRDS("df_results_cobit.rds")

nasim = df_results$simulation[which(is.na(df_results$beta1))]
nasim
df_results = df_results %>% filter(!simulation %in% nasim) %>%
  filter(simulation <= 1000 + length(nasim))
df_results$simulation[which(is.na(df_results$beta1))]

df_results %>%
  group_by(datagen, model, nsample) %>%
  summarise(bias = round(mean(beta1) -1,3),
            rmse = round(sqrt(mean((beta1 -1)^2)),3)) %>% as.data.frame()

df_results <- data.frame()



# A. logit link
nsim = 1100
for(n in c(100, 400, 1600)){
  dataall = readRDS(paste0("data_logit_n",n,".rds"))
  for(isim in 1:nsim){
    df1234 = dataall[[isim]]

    for(datagen in c("beta","cobin","betarec","beta3mix")){
      df = df1234[[datagen]]
      tryCatch({
        fit1 = glm.cobin(y ~ X, data = df, link = "logit", verbose = F)
      }, error = function(e) {
        message("Error in glm.cobin in iteration ", isim, " and n =", n,": ", e$message)
        fit1 <<- list(coefficients = c(NA, NA))
        return(NULL)
      })
      df_temp1 <- data.frame("simulation" = isim,
                             "datagen" = datagen,
                             "link" = "logit",
                             "model" = "cobin",
                             "nsample" = n,
                             "beta1" = fit1$coefficients[2])
      tryCatch({
        fit2 = betareg(y ~ X, data = df, link = "logit")
      }, error = function(e) {
        message("Error in betareg in iteration ", isim, " and n =", n,": ", e$message)
        fit2 <<- list(coefficients = list(mean = c(NA, NA)))
        return(NULL)
      })
      df_temp2 <- data.frame("simulation" = isim,
                             "datagen" = datagen,
                             "link" = "logit",
                             "model" = "beta",
                             "nsample" = n,
                             "beta1" = fit2$coefficients$mean[2])

      df_results <- rbind(df_results, df_temp1, df_temp2)
    }
    if(isim %% 100 == 0){
      # print n and isim
      print(paste0("n = ",n, "isim = ",isim))
    }
  }
}
#saveRDS(df_results,"df_results_logit.rds")
#df_results = readRDS("df_results_logit.rds")

nasim = df_results$simulation[which(is.na(df_results$beta1))]
nasim
df_results = df_results %>% filter(!simulation %in% nasim) %>%
  filter(simulation <= 1000 + length(nasim))
df_results$simulation[which(is.na(df_results$beta1))]

df_results %>%
  group_by(datagen, model, nsample) %>%
  summarise(bias = round(mean(beta1) -1,3),
            rmse = round(sqrt(mean((beta1 -1)^2)),3)) %>% as.data.frame()

