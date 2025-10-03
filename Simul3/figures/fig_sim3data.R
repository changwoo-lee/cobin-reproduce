rm(list = ls())

library(ggplot2)
library(cobin)
library(betareg)
library(ggpubr)

ngrid   = 5000
xgrid   = seq(0.00001, 0.99999, length = ngrid)

mu_vals_base <- c(
  cobin::bftprime(-6 + qnorm(0.01)*3),  # 1%
  cobin::bftprime(-6),                  # median
  cobin::bftprime(-6 + qnorm(0.99)*3)   # 99% 
)

mu_all    <- c(mu_vals_base, 0.5)

# legend labels (aligned to mu_all order)
mu_labels <- c("1% quantile", "50% quantile", "99% quantile", "outlier")

# lock legend ordering
mu_labels <- factor(mu_labels, levels = mu_labels)

# color specification (named vector)
cols <- c(
  "1% quantile"        = "#F8766D",
  "50% quantile"       = "#00BA38",
  "99% quantile"       = "#619CFF",
  "outlier" = "grey40"  # forest green
)

phi_beta <- 17
df_beta <- do.call(rbind, lapply(seq_along(mu_all), function(i){
  mu <- mu_all[i]
  data.frame(
    x = xgrid,
    y = dbeta(xgrid, mu * phi_beta, (1 - mu) * phi_beta),
    mu_label = mu_labels[i]
  )
}))

p_beta <- ggplot(df_beta, aes(x = x, y = y, color = mu_label)) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(x = 0.001, y = -0.3), shape = 17, size = 2, stroke = 1, inherit.aes = FALSE) +  # "x" mark at (0.01,0)
  scale_color_manual(values = cols, drop = FALSE) +
  labs(y = "density", x = "y", color = NULL) +
  facet_wrap(~"Beta distribution") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key.width = unit(1.6, "lines"))

p_beta

lambda_cobin <- 6
df_cobin <- do.call(rbind, lapply(seq_along(mu_all), function(i){
  mu  <- mu_all[i]
  eta <- bftprimeinv(mu)
  data.frame(
    x = xgrid,
    y = dcobin(xgrid, eta, lambda_cobin),
    mu_label = mu_labels[i]
  )
}))

p_cobin <- ggplot(df_cobin, aes(x = x, y = y, color = mu_label)) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(x = 0.001, y = -0.5), shape = 17, size = 2, stroke = 1, inherit.aes = FALSE) +  # "x" mark at (0.01,0)
  scale_color_manual(values = cols, drop = FALSE) +
  labs(y = "density", x = "y", color = NULL) +
  facet_wrap(~"Cobin distribution") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key.width = unit(1.6, "lines"))

final_plot <- ggarrange(p_beta, p_cobin, ncol = 2, common.legend = TRUE, legend = "right")
final_plot
#ggsave("fig_sim3data.pdf", final_plot, width = 10, height = 3)


