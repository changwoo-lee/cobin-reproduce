
library(ggplot2)

# log unnormalized density: lambda * (x - n * bft(x))
log_unnorm <- function(x, n, lambda) {
  lambda * (x - n * cobin::bft(x))
}

# Normalizing constant Z
compute_Z <- function(n, lambda) {
  integrate(
    f    = function(x) exp(log_unnorm(x, n, lambda)),
    lower = -Inf, upper = Inf,
    rel.tol = 1e-10
  )$value
}

# Proper pdf using the normalizer
d_bft <- function(x, n, lambda, Z) {
  # stable: exp(log f - log Z)
  exp(log_unnorm(x, n, lambda) - log(Z))
}

n_vals      <- c(20, 40, 80)
lambda_vals <- c(1, 2)

norm_tbl <- do.call(
  rbind,
  lapply(n_vals, function(n) {
    data.frame(
      n = n,
      lambda = lambda_vals,
      Z = sapply(lambda_vals, function(lam) compute_Z(n, lam))
    )
  })
)

x_grid <- seq(-120, 0, length.out = 2000)

dens_df_list <- list()
idx <- 1
for (n in n_vals) {
  for (lam in lambda_vals) {
    Zval <- norm_tbl$Z[norm_tbl$n == n & norm_tbl$lambda == lam]
    dens_df_list[[idx]] <- data.frame(
      x = x_grid,
      n = factor(n),
      lambda = factor(lam),
      density = d_bft(x_grid, n, lam, Zval)
    )
    idx <- idx + 1
  }
}
dens_df <- do.call(rbind, dens_df_list)

gg = ggplot(dens_df, aes(x = x, y = density, color = n, linetype = lambda)) +
  geom_line() +
  labs(
    x = expression(beta),
    y = "Density",
    color = "n",
    linetype = expression(lambda)
  ) +
  coord_cartesian(xlim = c(-120, 0), ylim = c(0, 0.13)) +
  theme_minimal(base_size = 12)
gg
#ggsave("target.pdf",gg,width = 10, height = 3.5)


library(cobin)

n = 20
set.seed(1)
y20_lambda1 = cobin::rcobin(n, -n, lambda = 1)
summary(y20_lambda1)
y20_lambda10 = cobin::rcobin(n, -n, lambda = 10)
y20_lambda10 = y20_lambda10 + (1/n - mean(y20_lambda10))
summary(y20_lambda10)


n = 40
y40_lambda1 = cobin::rcobin(n, -n, lambda = 1)
y40_lambda1 = y40_lambda1 + (1/n - mean(y40_lambda1))
summary(y40_lambda1)
y40_lambda10 = cobin::rcobin(n, -n, lambda = 10)
y40_lambda10 = y40_lambda10 + (1/n - mean(y40_lambda10))
summary(y40_lambda10)

set.seed(3)
n = 80
y80_lambda1 = cobin::rcobin(n, -n, lambda = 1)
y80_lambda1 = y80_lambda1 + (1/n - mean(y80_lambda1))
summary(y80_lambda1)
y80_lambda10 = cobin::rcobin(n, -n, lambda = 10)
y80_lambda10 = y80_lambda10 + (1/n - mean(y80_lambda10))
summary(y80_lambda10)

nsave = 100000
nburn = 1000
set.seed(1)
# y20, lambda = 1
fit_y20_lambda1 <- cobin::cobinreg(
  y ~ 1,
  data = data.frame(y = y20_lambda1),
  nsave = nsave, nburn = nburn,
  lambda_fixed = 1
)
# y20, lambda = 10
fit_y20_lambda10 <- cobin::cobinreg(
  y ~ 1,
  data = data.frame(y = y20_lambda10),
  nsave = nsave, nburn = nburn,
  lambda_fixed = 10
)

# y40, lambda = 1
fit_y40_lambda1 <- cobin::cobinreg(
  y ~ 1,
  data = data.frame(y = y40_lambda1),
  nsave = nsave, nburn = nburn,
  lambda_fixed = 1
)
fit_y40_lambda10 <- cobin::cobinreg(
  y ~ 1,
  data = data.frame(y = y40_lambda10),
  nsave = nsave, nburn = nburn,
  lambda_fixed = 10
)
# y80, lambda = 1
fit_y80_lambda1 <- cobin::cobinreg(
  y ~ 1,
  data = data.frame(y = y80_lambda1),
  nsave = nsave, nburn = nburn,
  lambda_fixed = 1
)
# y80, lambda = 10
fit_y80_lambda10 <- cobin::cobinreg(
  y ~ 1,
  data = data.frame(y = y80_lambda10),
  nsave = nsave, nburn = nburn,
  lambda_fixed = 10
)



library(dplyr)

get_trace_df <- function(fit, n, lambda) {
  data.frame(
    iter = 1:5000,
    value = fit$post_save[1:5000, 1],   # first column: intercept samples
    n = factor(n),
    lambda = as.character(lambda)
  )
}

# Combine traces
trace_df <- bind_rows(
  get_trace_df(fit_y20_lambda1, 20, 1),
  get_trace_df(fit_y20_lambda10, 20, 10),
  get_trace_df(fit_y40_lambda1, 40, 1),
  get_trace_df(fit_y40_lambda10, 40, 10),
  get_trace_df(fit_y80_lambda1, 80, 1),
  get_trace_df(fit_y80_lambda10, 80, 10)
)
#trace_df <- trace_df %>% rename(value = var1)
library(ggplot2)

g12 = ggplot(trace_df, aes(x = iter, y = value, color = n)) +
  geom_line(linewidth = 0.1) +
  facet_wrap(~ lambda, labeller = label_bquote(lambda == .(lambda))) +
  labs(
    x = "Iteration",
    y = expression(beta),
    color = "n"
  ) +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
g12


get_acf_df <- function(fit, n, lambda, lag_max = 50) {
  acf_obj <- acf(fit$post_save[, 1], lag.max = lag_max, plot = FALSE)
  data.frame(
    lag = acf_obj$lag,
    acf = acf_obj$acf,
    n = factor(n),
    lambda = as.character(lambda)
  )
}


acf_df <- bind_rows(
  get_acf_df(fit_y20_lambda1, 20, 1),
  get_acf_df(fit_y20_lambda10, 20, 10),
  get_acf_df(fit_y40_lambda1, 40, 1),
  get_acf_df(fit_y40_lambda10, 40, 10),
  get_acf_df(fit_y80_lambda1, 80, 1),
  get_acf_df(fit_y80_lambda10, 80, 10)
)


ggplot(acf_df, aes(x = lag, y = acf, color = n)) +
  geom_line() +
  facet_wrap(~ lambda, labeller = label_bquote(lambda == .(lambda)), scales = "free_y") +
  labs(
    x = "Lag",
    y = "ACF",
    color = "n"
  ) +
  theme_minimal(base_size = 12)


g3 = ggplot(acf_df, aes(x = lag, y = acf, color = n, linetype = lambda)) +
  geom_line() +
  labs(
    x = "Lag",
    y = "ACF",
    color = "n",
    linetype = expression(lambda)
  ) +
  scale_linetype_manual(values = c("1" = "solid", "10" = "dashed")) +
  theme_minimal(base_size = 12) 

# save
library(patchwork)

g123 <- (g12 | g3) + plot_layout(widths = c(5, 2))
g123
#ggsave("boundaryproximate_trace.pdf",g123,width = 10, height = 3.5)
