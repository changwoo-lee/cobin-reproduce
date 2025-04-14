#rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# data: n949
df = read.csv("mmi_lakecat.csv")
df = df[df$MMI_BENT != 0,]
# also remove two rows with comid = 9201925, 22845861
#excludeidx = which(df$comid %in% c(9201925, 22845861))
#df[excludeidx,]
#df = df[-excludeidx,]

dim(df) # nrow = 949
o =  order(df$easting)
nsave = 5000
# quantile residual

fit_cobin1 = readRDS("results_main_n949/fit_cobin1_n949.rds")
fit_cobin2 = readRDS("results_main_n949/fit_cobin2_n949.rds")
fit_cobin3 = readRDS("results_main_n949/fit_cobin3_n949.rds")


fit_cobin_betasave = as.mcmc.list(list(fit_cobin1$post_save[,1:ncol(fit_cobin1$X)],
                                       fit_cobin2$post_save[,1:ncol(fit_cobin2$X)],
                                       fit_cobin3$post_save[,1:ncol(fit_cobin3$X)]))




fit_micobin1 = readRDS("results_main_n949/fit_micobin1_n949.rds")
fit_micobin2 = readRDS("results_main_n949/fit_micobin2_n949.rds")
fit_micobin3 = readRDS("results_main_n949/fit_micobin3_n949.rds")


fit_micobin_betasave = as.mcmc.list(list(fit_cobin1$post_save[,1:ncol(fit_micobin1$X)],
                                       fit_cobin2$post_save[,1:ncol(fit_micobin2$X)],
                                       fit_cobin3$post_save[,1:ncol(fit_micobin3$X)]))

fit_beta1 = readRDS("results_main_n949/fit_beta1_n949.rds")
fit_beta2 = readRDS("results_main_n949/fit_beta2_n949.rds")
fit_beta3 = readRDS("results_main_n949/fit_beta3_n949.rds")

out1 = rstan::extract(fit_beta1, permute = FALSE)
fit_beta1_betasave = cbind(out1[,1,"b_Intercept"] ,out1[,1,grep("^b\\[\\d+\\]$", dimnames(out1)$parameters)])
fit_beta1_usave = out1[,1,grep("^u\\[\\d+\\]$", dimnames(out1)$parameters)]
fit_beta1_phisave = out1[,1,"phi"]

out2 = rstan::extract(fit_beta2, permute = FALSE)
fit_beta2_betasave = cbind(out2[,1,"b_Intercept"] ,out2[,1,grep("^b\\[\\d+\\]$", dimnames(out2)$parameters)])
fit_beta2_usave = out2[,1,grep("^u\\[\\d+\\]$", dimnames(out2)$parameters)]
fit_beta2_phisave = out2[,1,"phi"]

out3 = rstan::extract(fit_beta3, permute = FALSE)
fit_beta3_betasave = cbind(out3[,1,"b_Intercept"] ,out3[,1,grep("^b\\[\\d+\\]$", dimnames(out3)$parameters)])
fit_beta3_usave = out3[,1,grep("^u\\[\\d+\\]$", dimnames(out3)$parameters)]
fit_beta3_phisave = out3[,1,"phi"]

X = fit_cobin1$X

fit_beta1_linpredsave = matrix(0, nsave, nrow(df))
fit_beta1_logliksave = matrix(0, nsave, nrow(df))
fit_beta2_linpredsave = matrix(0, nsave, nrow(df))
fit_beta2_logliksave = matrix(0, nsave, nrow(df))
fit_beta3_linpredsave = matrix(0, nsave, nrow(df))
fit_beta3_logliksave = matrix(0, nsave, nrow(df))
for(isave in 1:nsave){
  fit_beta1_linpredsave[isave,] = as.numeric(X %*% fit_beta1_betasave[isave,] + fit_beta1_usave[isave,order(o)]) # permute back random effect
  fit_beta1_logliksave[isave,] = betareg::dbetar(df$MMI_BENT,
                                                 cobin::bftprime(fit_beta1_linpredsave[isave,]),
                                                 fit_beta1_phisave[isave], log = TRUE)
  fit_beta2_linpredsave[isave,] = as.numeric(X %*% fit_beta2_betasave[isave,] + fit_beta2_usave[isave,order(o)]) # permute back random effect
  fit_beta2_logliksave[isave,] = betareg::dbetar(df$MMI_BENT,
                                                 cobin::bftprime(fit_beta2_linpredsave[isave,]),
                                                 fit_beta2_phisave[isave], log = TRUE)
  fit_beta3_linpredsave[isave,] = as.numeric(X %*% fit_beta3_betasave[isave,] + fit_beta3_usave[isave,order(o)]) # permute back random effect
  fit_beta3_logliksave[isave,] = betareg::dbetar(df$MMI_BENT,
                                                 cobin::bftprime(fit_beta3_linpredsave[isave,]),
                                                 fit_beta3_phisave[isave], log = TRUE)
}

#################################################################


fit_cobin_betahat = colMeans(as.matrix(fit_cobin_betasave))
fit_cobin_linpredhat = fit_cobin1$X%*%fit_cobin_betahat + colMeans(rbind(fit_cobin1$post_u_save, fit_cobin2$post_u_save, fit_cobin3$post_u_save))
fit_cobin_muhat = cobin::bftprime(fit_cobin_linpredhat)
fit_cobin_lambdahat = median(c(fit_cobin1$post_save[,"lambda"], fit_cobin2$post_save[,"lambda"],fit_cobin3$post_save[,"lambda"] ))

fit_micobin_betahat = colMeans(as.matrix(fit_micobin_betasave))
fit_micobin_linpredhat = fit_cobin1$X%*%fit_micobin_betahat + colMeans(rbind(fit_micobin1$post_u_save, fit_micobin2$post_u_save, fit_micobin3$post_u_save))
fit_micobin_muhat = cobin::bftprime(fit_micobin_linpredhat)
fit_micobin_psihat = mean(c(fit_micobin1$post_save[,"psi"], fit_micobin2$post_save[,"psi"],fit_micobin3$post_save[,"psi"] ))

fit_beta_muhat = cobin::bftprime(colMeans(rbind(fit_beta1_linpredsave,
                                                fit_beta2_linpredsave,
                                                fit_beta3_linpredsave)))
fit_beta_phihat = mean(c(fit_beta1_phisave,fit_beta2_phisave, fit_beta3_phisave))




# quantile residual
cobin_qresid = numeric(nrow(df))
micobin_qresid = numeric(nrow(df))
for(i in 1:nrow(df)){
  cobin_qresid[i] = cobin::pcobin(df$MMI_BENT[i], fit_cobin_linpredhat[i], fit_cobin_lambdahat)
  micobin_qresid[i] = cobin::pmicobin(df$MMI_BENT[i], fit_micobin_linpredhat[i], fit_micobin_psihat)
}
beta_qresid = betareg::pbetar(df$MMI_BENT, fit_beta_muhat, fit_beta_phihat)




par(mfrow = c(1,3))
qqnorm(qnorm(beta_qresid), main = "Quantile residual, betareg", pch = 16)
abline(a=0,b=1)
idx12_beta = order(qnorm(beta_qresid))[1:2]
idx12_beta

qqnorm(qnorm(cobin_qresid), main = "Quantile residual, cobin", pch = 16)
abline(a=0,b=1)
idx12_cobin = order(qnorm(cobin_qresid))[1:2]
idx12_cobin

qqnorm(qnorm(micobin_qresid), main = "Quantile residual, micobin", pch = 16)
abline(a=0,b=1)
idx12_micobin = order(qnorm(micobin_qresid))[1:2]
idx12_micobin

# idx12 = idx12_beta
# 
# df[idx12,"comid"] # unique identifer
# df[idx12,]
# 
# 
# 
library(ggplot2)

n <- nrow(df)
n
myshape1 = rep(20,n)

myshape1[idx12[1]] = 15
myshape1[idx12[2]] = 17

myshape = rep(myshape1, 3)

myshape

drawdf <- data.frame(
  model = rep(c("beta regression", "cobin regression", "micobin regression"), each = n),
  qresid = c(qnorm(beta_qresid), qnorm(cobin_qresid), qnorm(micobin_qresid)),
  myshape = myshape
)

drawdf <- drawdf %>%
  group_by(model) %>%
  arrange(qresid, .by_group = TRUE) %>%
  mutate(prob = (row_number() - 0.5) / n(),
         theo = qnorm(prob)) %>%
  ungroup()

figqq <- ggplot(drawdf, aes(x = theo, y = qresid)) +
  geom_point(aes(shape = myshape)) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~ model) +
  labs(x = "Standard normal quantile",
       y = "Quantile residuals") +
  theme_bw() +
 # scale_shape_identity() + 
  scale_shape_identity(
    guide = "legend",
    breaks = c(15, 17), 
    labels = c("\nA lake with MMI = 0.020\n(Jones pond, Anson, NC)\n",
               "\nA lake with MMI = 0.021\n(Ferguson lake, Saline, AR)\n")
  ) + theme(legend.title = element_blank())

print(figqq)

ggsave("qresidplot.pdf", figqq, width = 9, height = 2.5)




















