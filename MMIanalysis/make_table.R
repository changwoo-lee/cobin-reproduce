
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# n = 949
res = readRDS(file = "n949/res/summary_beta_n949.rds")
paste(round(res$coef_orig_mean, 3),"(", round(res$coef_orig_95lo, 3),",",round(res$coef_orig_95up, 3),")" )

res = readRDS(file = "n949/res/summary_cobin_n949.rds")
paste(round(res$coef_orig_mean, 3),"(", round(res$coef_orig_95lo, 3),",",round(res$coef_orig_95up, 3),")" )

res = readRDS(file = "n949/res/summary_betarec_n949.rds")
paste(round(res$coef_orig_mean, 3),"(", round(res$coef_orig_95lo, 3),",",round(res$coef_orig_95up, 3),")" )

res = readRDS(file = "n949/res/summary_micobin_n949.rds")
paste(round(res$coef_orig_mean, 3),"(", round(res$coef_orig_95lo, 3),",",round(res$coef_orig_95up, 3),")" )


# n = 947
res = readRDS(file = "n947/res/summary_beta_n947.rds")
paste(round(res$coef_orig_mean, 3),"(", round(res$coef_orig_95lo, 3),",",round(res$coef_orig_95up, 3),")" )

res = readRDS(file = "n947/res/summary_cobin_n947.rds")
paste(round(res$coef_orig_mean, 3),"(", round(res$coef_orig_95lo, 3),",",round(res$coef_orig_95up, 3),")" )

res = readRDS(file = "n947/res/summary_betarec_n947.rds")
paste(round(res$coef_orig_mean, 3),"(", round(res$coef_orig_95lo, 3),",",round(res$coef_orig_95up, 3),")" )

res = readRDS(file = "n947/res/summary_micobin_n947.rds")
paste(round(res$coef_orig_mean, 3),"(", round(res$coef_orig_95lo, 3),",",round(res$coef_orig_95up, 3),")" )

# n = 950
res = readRDS(file = "n950/res/summary_betarec_n950.rds")
paste(round(res$coef_orig_mean, 3),"(", round(res$coef_orig_95lo, 3),",",round(res$coef_orig_95up, 3),")" )

res = readRDS(file = "n950/res/summary_micobin_n950.rds")
paste(round(res$coef_orig_mean, 3),"(", round(res$coef_orig_95lo, 3),",",round(res$coef_orig_95up, 3),")" )


# difference between 949 and 947

# beta
res1 = readRDS(file = "n949/res/summary_beta_n949.rds")
res2 = readRDS(file = "n947/res/summary_beta_n947.rds")
round(sqrt(sum((res1$coef_orig_mean - res2$coef_orig_mean)^2)),3)

# cobin
res1 = readRDS(file = "n949/res/summary_cobin_n949.rds")
res2 = readRDS(file = "n947/res/summary_cobin_n947.rds")
round(sqrt(sum((res1$coef_orig_mean - res2$coef_orig_mean)^2)),3)

# betarec
res1 = readRDS(file = "n949/res/summary_betarec_n949.rds")
res2 = readRDS(file = "n947/res/summary_betarec_n947.rds")
round(sqrt(sum((res1$coef_orig_mean - res2$coef_orig_mean)^2)),3)

# micobin
res1 = readRDS(file = "n949/res/summary_micobin_n949.rds")
res2 = readRDS(file = "n947/res/summary_micobin_n947.rds")
round(sqrt(sum((res1$coef_orig_mean - res2$coef_orig_mean)^2)),3)




# betarec
res1 = readRDS(file = "n950/res/summary_betarec_n950.rds")
res2 = readRDS(file = "n947/res/summary_betarec_n947.rds")
round(sqrt(sum((res1$coef_orig_mean - res2$coef_orig_mean)^2)),3)

# micobin
res1 = readRDS(file = "n950/res/summary_micobin_n950.rds")
res2 = readRDS(file = "n947/res/summary_micobin_n947.rds")
round(sqrt(sum((res1$coef_orig_mean - res2$coef_orig_mean)^2)),3)








