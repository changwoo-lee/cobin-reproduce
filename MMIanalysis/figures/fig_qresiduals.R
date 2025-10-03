rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# data: n949
df = read.csv("../mmi_lakecat.csv")
df = df[df$MMI_BENT != 0,]

dim(df) # nrow = 949
o =  order(df$easting)
nsave = 5000


source("qresid_functions.R")

res_beta = readRDS(file = "../n949/res/summary_beta_n949.rds")
res_cobin = readRDS(file = "../n949/res/summary_cobin_n949.rds")
res_betarec = readRDS(file = "../n949/res/summary_betarec_n949.rds")
res_micobin = readRDS(file = "../n949/res/summary_micobin_n949.rds")

# quantile residual

mu_beta = cobin::bftprime(res_beta$linpred_mean)
phi_beta = res_beta$phi_mean

eta_cobin = res_cobin$linpred_mean
lambda_cobin = res_cobin$lambda_median

mu_betarec = cobin::bftprime(res_betarec$linpred_mean)
phi_betarec = res_betarec$phi_mean
alpha_betarec = res_betarec$alpha_mean

eta_micobin = res_micobin$linpred_mean
psi_micobin = res_micobin$psi_mean

beta_qresid = qres_beta(y = df$MMI_BENT, mu = mu_beta, phi = phi_beta)
cobin_qresid = qres_cobin(y = df$MMI_BENT, eta = eta_cobin, lambda = lambda_cobin)
betarec_qresid = qres_betarec(y = df$MMI_BENT, mu = mu_betarec, alpha = alpha_betarec, phi = phi_betarec)
micobin_qresid = qres_micobin(y = df$MMI_BENT, eta = eta_micobin, psi = psi_micobin)


idx12 = c(147,34)
 # df[idx12,"comid"] # unique identifer
# df[idx12,]
library(ggplot2)

n <- nrow(df)
n
myshape1 = rep(20,n)

myshape1[idx12[1]] = 15
myshape1[idx12[2]] = 17

myshape = rep(myshape1, 4)

myshape

drawdf <- data.frame(
  model = rep(c(" Beta regression", " Cobin regression", "Betarec regression", "Micobin regression"), each = n),
  qresid = c(beta_qresid, cobin_qresid, betarec_qresid, micobin_qresid),
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
  facet_wrap(~ model, nrow = 1) +
  labs(x = "Standard normal quantile",
       y = "Quantile residuals") +
  theme_bw() +
 # scale_shape_identity() + 
  scale_shape_identity(
    guide = "legend",
    breaks = c(15, 17), 
    labels = c("\nA lake with \nMMI = 0.020\n(Jones pond, \nAnson, NC)\n",
               "\nA lake with \nMMI = 0.021\n(Ferguson lake, \nSaline, AR)\n")
  ) + theme(legend.title = element_blank())

print(figqq)

#ggsave("qresidplot_4.pdf", figqq, width = 9, height = 2.5)




















