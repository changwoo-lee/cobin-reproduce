
# Quantile residuals for Beta regression
qres_beta <- function(y, mu, phi, eps = 1e-12) {
  p <- betareg::pbetar(y, mu, phi)
  p <- pmin(pmax(p, eps), 1 - eps)
  qnorm(p)
}

# Quantile residuals for cobin regression
qres_cobin <- function(y, eta, lambda, eps = 1e-12) {
  p = numeric(length(y))
  for(i in seq_along(y)) {
    p[i] <- cobin::pcobin(y[i], eta[i], lambda)
  }
  p <- pmin(pmax(p, eps), 1 - eps)
  qnorm(p)
}


## ---------- Beta-rectangular CDF ----------
pbetarec <- function(y, mu, phi, alpha) {
  # recycle to common length
  y   <- as.numeric(y)
  mu  <- rep_len(mu,  length(y))
  phi <- rep_len(phi, length(y))
  alpha <- rep_len(alpha, length(y))
  
  # parameters as in dbetarec()
  gam   <- mu
  theta <- alpha * (1 - abs(2 * gam - 1))
  mu_adj <- (gam - 0.5 * theta) / (1 - theta)
  a <- mu_adj * phi
  b <- (1 - mu_adj) * phi
  
  F_beta <- pbeta(y, a, b)
  F_unif <- punif(y,0,1)
  
  # mixture CDF on [0,1]
  F_mix <- (1 - theta) * F_beta + theta * F_unif
  
  # enforce 0/1 outside support
  F_mix[y <= 0] <- 0
  F_mix[y >= 1] <- 1
  
  F_mix
}

# Quantile residuals for beta rectangular regression
qres_betarec <- function(y, mu, phi, alpha, eps = 1e-12) {
  p <- pbetarec(y, mu, phi, alpha)
  p <- pmin(pmax(p, eps), 1 - eps)
  qnorm(p)
}

# Quantile residuals for micobin regression
qres_micobin <- function(y, eta, psi, eps = 1e-12) {
  p = numeric(length(y))
  for(i in seq_along(y)) {
    p[i] <- cobin::pmicobin(y[i], eta[i], psi)
  }
  p <- pmin(pmax(p, eps), 1 - eps)
  qnorm(p)
}


draw_qqplot <- function(qresid,
                        model = NULL,                
                        model_name = "Model",      
                        highlight_idx = NULL,       
                        highlight_labels = NULL,    
                        shape_outlier = c(15, 17)) {       
  
  # Dependencies
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install.packages('dplyr')")
  library(dplyr)
  library(ggplot2)
  
  
  # Build base data.frame
  n <- length(qresid)
  if (is.null(model)) {
    model <- rep(model_name, n)
  } else {
    if (length(model) != n) stop("length(model) must equal length(qresid)")
  }
  
  shape_vec <- rep(16L, n)
  if (!is.null(highlight_idx)) {
    if (length(highlight_idx) >= 1) shape_vec[highlight_idx[1]] <- 15L
    if (length(highlight_idx) >= 2) shape_vec[highlight_idx[2]] <- 17L
  }
  
  df <- data.frame(model = factor(model),
                   qresid = qresid,
                   shape  = shape_vec,
                   stringsAsFactors = FALSE)
  
  drawdf <- df %>%
    dplyr::group_by(model) %>%
    dplyr::arrange(qresid, .by_group = TRUE) %>%
    dplyr::mutate(prob = (dplyr::row_number() - 0.5) / dplyr::n(),
                  theo = qnorm(prob)) %>%
    dplyr::ungroup()
  
  # Build plot
  p <- ggplot(drawdf, aes(x = theo, y = qresid)) +
    geom_point(aes(shape = factor(shape))) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    facet_wrap(~ model) +
    labs(x = "Standard normal quantile",
         y = "Quantile residuals") +
    theme_bw() +
    scale_shape_manual(
      values = c("16" = 16, "15" = shape_outlier[1], "17" = shape_outlier[2]),
      breaks = c("15", "17"),
      labels = c(
        "15" = if (!is.null(highlight_labels) && length(highlight_labels) >= 1)
          highlight_labels[1] else "Highlight 1",
        "17" = if (!is.null(highlight_labels) && length(highlight_labels) >= 2)
          highlight_labels[2] else "Highlight 2"
      ),
      guide = "legend"
    ) + ylim(-6.1,3)+
    theme(legend.title = element_blank())
  
  return(p)
}

