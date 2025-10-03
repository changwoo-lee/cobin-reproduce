
data = readRDS("../data/data_outliersim.rds")
source("qresid_functions.R")
df_outlier = data[[44]]$beta_y10m3
df = df_outlier[-1,]
library(cobin)
library(rstan)
library(brms)

nsave = 5000
nburn = 1000
nthin = 1
######## beta regression ##########
# whole dataset
stan_data_o = brms::standata(brm(y ~ ., data = df_outlier, 
                                 family = Beta(), empty = TRUE) )
fit_beta_o = rstan::stan("../betareg_fixedeffect.stan",
                         data = stan_data_o,
                         chains = 1, seed = 4,
                         iter = nburn + nsave, warmup = nburn, thin = nthin)

betahat_o = summary(fit_beta_o)$summary[c("b_Intercept","b[1]","b[2]"),"mean"]
phi_o = summary(fit_beta_o)$summary["phi","mean"]

beta_qresid = qres_beta(df_outlier$y, 
                        mu = cobin::bftprime(cbind(1, df_outlier$x1, df_outlier$x2)%*%betahat_o),
                        phi = phi_o)

qq_beta  <- draw_qqplot(beta_qresid, model_name = "Beta regression", highlight_idx = 1,
                        highlight_labels = c("outlier"))
qq_beta


######## beta rect regression ##########
# whole dataset
fit_betarec_o = rstan::stan("../betarecreg_fixedeffect.stan",
                            data = stan_data_o,
                            chains = 1, seed = 44,
                            iter = nburn + nsave, warmup = nburn, thin = nthin)
betahat_o = summary(fit_betarec_o)$summary[c("b_Intercept","b[1]","b[2]"),"mean"]
phi_o = summary(fit_betarec_o)$summary["phi","mean"]
alpha_o = summary(fit_betarec_o)$summary["alpha","mean"]
betarec_qresid = qres_betarec(df_outlier$y, 
                              mu = cobin::bftprime(cbind(1, df_outlier$x1, df_outlier$x2)%*%betahat_o),
                              phi = phi_o,
                              alpha = 0.1)
qq_betarec  <- draw_qqplot(betarec_qresid, model_name = "Betarec regression", highlight_idx = 1,
                           highlight_labels = c("outlier"))
alpha_o
qq_betarec

######## cobin regression ##########
set.seed(6)
fit_cobin_o = cobin::cobinreg(y ~ ., data = df_outlier, link = "cobit",
                              nburn = nburn, nsave = nsave, nthin = nthin)

betahat_o = colMeans(fit_cobin_o$post_save[,1:3])
lambdahat_o = apply(fit_cobin_o$post_save[,"lambda",drop=F], 2, median)

cobin_qresid = qres_cobin(df_outlier$y, 
                          cbind(1, df_outlier$x1, df_outlier$x2)%*%betahat_o,
                          lambda = lambdahat_o)
qq_cobin  <- draw_qqplot(cobin_qresid, model_name = "Cobin regression", highlight_idx = 1,
                         highlight_labels = c("outlier"))
qq_cobin

######## micobin regression ##########
set.seed(55)
fit_micobin_o = cobin::micobinreg(y ~ ., data = df_outlier, link = "cobit", 
                                  nburn = nburn, nsave = nsave, nthin = nthin)
betahat_o = colMeans(fit_micobin_o$post_save[,1:3])
psihat_o = apply(fit_micobin_o$post_save[,"psi",drop=F], 2, mean)
micobin_qresid = qres_micobin(df_outlier$y, 
                              cbind(1, df_outlier$x1, df_outlier$x2)%*%betahat_o,
                              psi = psihat_o)
qq_micobin  <- draw_qqplot(micobin_qresid, model_name = "Micobin regression", highlight_idx = 1,
                           highlight_labels = c("outlier"))
qq_micobin


# 
# draw_qqplot_faceted <- function(qresid_list,
#                                 model_names = NULL,
#                                 highlight_idx_list = NULL,   
#                                 highlight_labels = c("Highlight 1", "Highlight 2"),
#                                 qres_is_prob = FALSE,
#                                 ylim_vals = c(-6, 6)) {
#   stopifnot(is.list(qresid_list), length(qresid_list) > 0)
#   
#   if (is.null(model_names)) {
#     model_names <- paste("Model", seq_along(qresid_list))
#   }
#   stopifnot(length(model_names) == length(qresid_list))
#   
#   # Build stacked data
#   make_block <- function(qres, model, hi_idx) {
#     # z-scale residuals
#     zres <- if (qres_is_prob) {
#       p <- pmin(pmax(qres, 1e-12), 1 - 1e-12)
#       qnorm(p)
#     } else {
#       as.numeric(qres)
#     }
#     n <- length(zres)
#     
#     shape_vec <- rep(16L, n)
#     if (!is.null(hi_idx)) {
#       hi_idx <- hi_idx[hi_idx >= 1 & hi_idx <= n]
#       if (length(hi_idx) >= 1) shape_vec[hi_idx[1]] <- 15L
#       if (length(hi_idx) >= 2) shape_vec[hi_idx[2]] <- 17L
#     }
#     
#     data.frame(model = model, qresid = zres, shape = shape_vec, stringsAsFactors = FALSE)
#   }
#   
#   df <- do.call(rbind, Map(make_block, qresid_list, model_names, 
#                            if (is.null(highlight_idx_list)) vector("list", length(qresid_list)) else highlight_idx_list))
#   
#   # Compute theo quantiles per model (your recipe)
#   suppressPackageStartupMessages({
#     library(dplyr); library(ggplot2)
#   })
#   drawdf <- df %>%
#     group_by(model) %>%
#     arrange(qresid, .by_group = TRUE) %>%
#     mutate(prob = (row_number() - 0.5) / n(),
#            theo = qnorm(prob)) %>%
#     ungroup()
#   
#   # Build plot
#   ggplot(drawdf, aes(x = theo, y = qresid)) +
#     geom_point(aes(shape = factor(shape))) +
#     geom_abline(intercept = 0, slope = 1, color = "red") +
#     facet_wrap(~ model, scales = "fixed", nrow = 1) +
#     labs(x = "Standard normal quantile",
#          y = "Quantile residuals") +
#     theme_bw() +
#     scale_shape_manual(
#       values = c("16" = 16, "15" = 17, "17" = 17),
#       breaks = c("15", "17"), 
#       labels = c(
#         "15" = if (length(highlight_labels) >= 1) highlight_labels[1] else "Highlight 1",
#         "17" = if (length(highlight_labels) >= 2) highlight_labels[2] else "Highlight 2"
#       ),
#       guide = "legend"
#     ) +
#     coord_cartesian(ylim = ylim_vals)  +
#     theme(legend.title = element_blank(), legend.position = "right")
# }
# 
# 
# p_faceted <- draw_qqplot_faceted(
#   qresid_list = list(beta_qresid, cobin_qresid, betarec_qresid, micobin_qresid),
#   model_names = c(" Beta regression", " Cobin regression", "Betarec regression", "Micobin regression"),
#   highlight_idx_list = list(c(1), c(1), c(1), c(1)),  # highlight the first obs in each panel
#   highlight_labels = c("outlier"),
#   qres_is_prob = FALSE,        # set TRUE if you pass CDF values instead of z-residuals
#   ylim_vals = c(-6.1, 3)       # match your previous limit if desired
# )
# 
# print(p_faceted)
# ggsave("fig_sim2data_qqplot_faceted.pdf", p_faceted, width = 9, height = 2.5)
# 
# 
