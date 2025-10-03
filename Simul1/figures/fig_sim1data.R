rm(list = ls())

library(ggplot2)
library(cobin)
library(betareg)
library(ggpubr)

dbetarec <- function(x, mu, phi, alpha, log = FALSE){
  gam   = mu
  theta = alpha * (1 - abs(2 * gam - 1))
  mu_adj = (gam - 0.5 * theta) / (1 - theta)
  loglik = log( (1 - theta) * dbeta(x, mu_adj * phi, (1 - mu_adj) * phi) +
                  theta * dunif(x, 0, 1) )
  if (log) loglik else exp(loglik)
}

dbeta3 <- function(x, mu, phi, log = FALSE){
  mu_left  = mu - pmin(mu, 1 - mu) / 2
  mu_right = mu + pmin(mu, 1 - mu) / 2
  comp_left  = betareg::dbetar(x, mu_left,  phi)
  comp_right = betareg::dbetar(x, mu_right, phi)
  comp_mid   = betareg::dbetar(x, mu,       phi)
  dens = 0.25 * comp_left + 0.25 * comp_right + 0.50 * comp_mid
  if (log) log(dens) else dens
}

ngrid   = 5000
xgrid   = seq(0.00001, 0.99999, length = ngrid)

mu_vals = c(cobin::bftprime(qnorm(0.01)*3), 0.5, cobin::bftprime( qnorm(0.99)*3))
mu_labels = c("1% quantile", "50% quantile", "99% quantile")

phi_beta = 8
df_beta <- do.call(rbind, lapply(seq_along(mu_vals), function(i){
  mu = mu_vals[i]
  data.frame(
    x = xgrid,
    y = dbeta(xgrid, mu * phi_beta, (1 - mu) * phi_beta),
    mu_label = mu_labels[i]
  )
}))

p_beta <- ggplot(df_beta, aes(x = x, y = y, color = mu_label)) +
  geom_line() +  ylim(0, 6) +
  labs(y = "density", x = "y") +
  facet_wrap(~"Beta (cobit)") +
  theme_bw() +
  theme(legend.title = element_blank())

lambda_cobin = 3
df_cobin <- do.call(rbind, lapply(seq_along(mu_vals), function(i){
  mu  = mu_vals[i]
  eta = bftprimeinv(mu)
  data.frame(
    x = xgrid,
    y = dcobin(xgrid, eta, lambda_cobin),
    mu_label = mu_labels[i]
  )
}))

p_cobin <- ggplot(df_cobin, aes(x = x, y = y, color = mu_label)) +
  geom_line() +  ylim(0, 6) +
  labs(y = "density", x = "y") +
  facet_wrap(~"Cobin (cobit)") +
  theme_bw() +
  theme(legend.title = element_blank())

phi_betarec = 10
alpha = 0.2
df_betarec <- do.call(rbind, lapply(seq_along(mu_vals), function(i){
  mu = mu_vals[i]
  data.frame(
    x = xgrid,
    y = dbetarec(xgrid, mu, phi_betarec, alpha),
    mu_label = mu_labels[i]
  )
}))

p_betarec <- ggplot(df_betarec, aes(x = x, y = y, color = mu_label)) +
  geom_line() +  ylim(0, 6) +
  labs(y = "density", x = "y") +
  facet_wrap(~"Betarec (cobit)") +
  theme_bw() +
  theme(legend.title = element_blank())

phi_beta3 = 40
df_beta3 <- do.call(rbind, lapply(seq_along(mu_vals), function(i){
  mu = mu_vals[i]
  data.frame(
    x = xgrid,
    y = dbeta3(xgrid, mu, phi_beta3),
    mu_label = mu_labels[i]
  )
}))

p_beta3 <- ggplot(df_beta3, aes(x = x, y = y, color = mu_label)) +
  geom_line() +  ylim(0, 6) +
  labs(y = "density", x = "y") +
  facet_wrap(~"Beta3mix (cobit)") +
  theme_bw() +
  theme(legend.title = element_blank())

final_plot1 = ggarrange(p_beta, p_cobin, p_betarec, p_beta3,
                       ncol = 4, common.legend = TRUE, legend = "right")
print(final_plot1)

#################################################################################

ngrid   = 5000
xgrid   = seq(0.0001, 0.9999, length = ngrid)
invlogit = function(x) exp(x) / (1 + exp(x))
mu_vals = c(invlogit(qnorm(0.01)), 0.5, invlogit(qnorm(0.99)))
mu_labels = c("1% quantile", "50% quantile", "99% quantile")

phi_beta = 8
df_beta <- do.call(rbind, lapply(seq_along(mu_vals), function(i){
  mu = mu_vals[i]
  data.frame(
    x = xgrid,
    y = dbeta(xgrid, mu * phi_beta, (1 - mu) * phi_beta),
    mu_label = mu_labels[i]
  )
}))

p_beta <- ggplot(df_beta, aes(x = x, y = y, color = mu_label)) +
  geom_line() +  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10))+
  labs(y = "density", x = "y") +
  facet_wrap(~"Beta (logit)") +
  theme_bw() +
  theme(legend.title = element_blank())

lambda_cobin = 3
df_cobin <- do.call(rbind, lapply(seq_along(mu_vals), function(i){
  mu  = mu_vals[i]
  eta = bftprimeinv(mu)
  data.frame(
    x = xgrid,
    y = dcobin(xgrid, eta, lambda_cobin),
    mu_label = mu_labels[i]
  )
}))

p_cobin <- ggplot(df_cobin, aes(x = x, y = y, color = mu_label)) +
  geom_line() + scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10))+
  labs(y = "density", x = "y") +
  facet_wrap(~"Cobin (logit)") +
  theme_bw() +
  theme(legend.title = element_blank())

phi_betarec = 10
alpha = 0.2
df_betarec <- do.call(rbind, lapply(seq_along(mu_vals), function(i){
  mu = mu_vals[i]
  data.frame(
    x = xgrid,
    y = dbetarec(xgrid, mu, phi_betarec, alpha),
    mu_label = mu_labels[i]
  )
}))

p_betarec <- ggplot(df_betarec, aes(x = x, y = y, color = mu_label)) +
  geom_line() + scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10))+
  labs(y = "density", x = "y") +
  facet_wrap(~"Betarec (logit)") +
  theme_bw() +
  theme(legend.title = element_blank())

phi_beta3 = 40
df_beta3 <- do.call(rbind, lapply(seq_along(mu_vals), function(i){
  mu = mu_vals[i]
  data.frame(
    x = xgrid,
    y = dbeta3(xgrid, mu, phi_beta3),
    mu_label = mu_labels[i]
  )
}))

p_beta3 <- ggplot(df_beta3, aes(x = x, y = y, color = mu_label)) +
  geom_line() + scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10))+
  labs(y = "density", x = "y") +
  facet_wrap(~"Beta3mix (logit)") +
  theme_bw() +
  theme(legend.title = element_blank())

final_plot2 = ggarrange(p_beta, p_cobin, p_betarec, p_beta3,
                        ncol = 4, common.legend = TRUE, legend = "right")
print(final_plot2)
final_plot = ggarrange(final_plot1, final_plot2, nrow = 2, ncol = 1, legend = "right")
final_plot
#ggsave("sim1_data_densities.pdf", final_plot, width = 12, height = 5)
