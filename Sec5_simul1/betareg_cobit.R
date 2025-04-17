# source: betareg R package
# including "cobit" link 
library(cobin) # required
library(Formula)
betareg <- function (formula, data, subset, na.action, weights, offset,
          link = c("logit", "probit", "cloglog", "cauchit", "log",
                   "loglog", "cobit"), link.phi = NULL, type = c("ML", "BC", "BR"),
          dist = NULL, nu = NULL, control = betareg.control(...), model = TRUE,
          y = TRUE, x = FALSE, ...)
{
  cl <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if (length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~1)
    simple_formula <- TRUE
  }
  else {
    if (length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1L:2L))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
  .add_predvars_and_dataClasses <- function(terms, model.frame) {
    rval <- terms
    nval <- if (inherits(model.frame, "terms"))
      model.frame
    else terms(model.frame, dot = control$dot)
    ovar <- sapply(as.list(attr(rval, "variables")), deparse)[-1]
    nvar <- sapply(as.list(attr(nval, "variables")), deparse)[-1]
    if (!all(ovar %in% nvar))
      stop(paste("The following terms variables are not part of the model.frame:",
                 paste(ovar[!(ovar %in% nvar)], collapse = ", ")))
    ix <- match(ovar, nvar)
    if (!is.null(attr(rval, "predvars")))
      warning("terms already had 'predvars' attribute, now replaced")
    attr(rval, "predvars") <- attr(nval, "predvars")[1L +
                                                       c(0L, ix)]
    if (!is.null(attr(rval, "dataClasses")))
      warning("terms already had 'dataClasses' attribute, now replaced")
    attr(rval, "dataClasses") <- attr(nval, "dataClasses")[ix]
    return(rval)
  }
  mt <- .add_predvars_and_dataClasses(mt, mf)
  mtX <- .add_predvars_and_dataClasses(mtX, mf)
  mtZ <- .add_predvars_and_dataClasses(mtZ, mf)
  if (length(Y) < 1)
    stop("empty model")
  if (any(Y < 0 | Y > 1))
    stop("invalid dependent variable, all observations must be in [0, 1]")
  n <- length(Y)
  weights <- model.weights(mf)
  if (is.null(weights))
    weights <- 1
  if (length(weights) == 1)
    weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  expand_offset <- function(offset) {
    if (is.null(offset))
      offset <- 0
    if (length(offset) == 1)
      offset <- rep.int(offset, n)
    as.vector(offset)
  }
  offsetX <- expand_offset(model.offset(model.part(formula,
                                                   data = mf, rhs = 1L, terms = TRUE)))
  offsetZ <- expand_offset(model.offset(model.part(formula,
                                                   data = mf, rhs = 2L, terms = TRUE)))
  if (!is.null(cl$offset))
    offsetX <- offsetX + expand_offset(mf[, "(offset)"])
  offset <- list(mu = offsetX, phi = offsetZ)
  type <- match.arg(type, c("ML", "BC", "BR"))
  if (!is.null(dist)) {
    dist <- tolower(as.character(dist))
    if (dist == "xb")
      dist <- "xbeta"
    if (dist == "xbx")
      dist <- "xbetax"
    dist <- match.arg(dist, c("beta", "xbeta", "xbetax"))
    if (dist == "beta") {
      if (any((Y <= 0 | Y >= 1)[weights > 0]))
        stop("dependent variable not suitable for 'beta' distribution, all observations must be in (0, 1)")
    }
    if (dist == "xbeta" && is.null(nu)) {
      warning(sprintf("estimation of 'nu' with 'xbeta' distribution is not feasible, using '%s' instead",
                      dist <- if (any((Y <= 0 | Y >= 1)[weights > 0]))
                        "xbetax"
                      else "beta"))
    }
  }
  else {
    if (is.null(nu)) {
      dist <- if (all((Y > 0 & Y < 1)[weights > 0]))
        "beta"
      else "xbetax"
    }
    else {
      dist <- "xbeta"
    }
  }
  if (dist != "beta" && type != "ML") {
    warning(sprintf("only 'ML' estimation is available for '%s' distribution",
                    dist))
    type <- "ML"
  }
  if (is.character(link))
    link <- match.arg(link)
  if (is.null(link.phi))
    link.phi <- if (simple_formula)
      "identity"
  else "log"
  if (is.character(link.phi))
    link.phi <- match.arg(link.phi, c("identity", "log",
                                      "sqrt"))
  rval <- betareg.fit(X, Y, Z, weights, offset, link, link.phi,
                      type, control, dist, nu)
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(mu = mtX, phi = mtZ, full = mt)
  rval$levels <- list(mu = .getXlevels(mtX, mf), phi = .getXlevels(mtZ,
                                                                   mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(mu = attr(X, "contrasts"), phi = attr(Z,
                                                               "contrasts"))
  if (model)
    rval$model <- mf
  if (y)
    rval$y <- Y
  if (x)
    rval$x <- list(mu = X, phi = Z)
  if (is.null(rval$dist) || (rval$dist == "beta")) {
    for (n in intersect(names(rval), betareg:::fix_names_mu_phi)) names(rval[[n]])[1L:2L] <- c("mean",
                                                                                     "precision")
  }
  class(rval) <- "betareg"
  return(rval)
}













betareg.fit <- function (x, y, z = NULL, weights = NULL, offset = NULL, link = "logit",
          link.phi = "log", type = "ML", control = betareg.control(),
          dist = NULL, nu = NULL)
{
  if (is.null(dist)) {
    if (is.null(nu)) {
      dist <- if (all(y > 0 & y < 1))
        "beta"
      else "xbetax"
    }
    else {
      dist <- "xbeta"
    }
  }
  if (dist != "beta") {
    if (type != "ML")
      stop(sprintf("only 'ML' estimation implemented for '%s'",
                   dist))
    control$hessian <- TRUE
    control$fsmaxit <- 0
  }
  estnu <- if (dist != "beta" && is.null(nu))
    TRUE
  else FALSE
  n <- NROW(x)
  k <- NCOL(x)
  if (is.null(weights))
    weights <- rep.int(1, n)
  nobs <- sum(weights > 0)
  if (is.null(offset))
    offset <- rep.int(0, n)
  if (!is.list(offset))
    offset <- list(mu = offset, phi = rep.int(0, n))
  if (is.null(z)) {
    m <- 1L
    z <- matrix(1, ncol = m, nrow = n)
    colnames(z) <- "(Intercept)"
    rownames(z) <- rownames(x)
    phi_const <- TRUE
  }
  else {
    m <- NCOL(z)
    if (m < 1L)
      stop("dispersion regression needs to have at least one parameter")
    phi_const <- (m == 1L) && isTRUE(all.equal(as.vector(z[,
                                                           1L]), rep.int(1, n)))
  }
  if (is.character(link)) {
    linkstr <- link
    if (linkstr != "loglog" & linkstr != "cobit") {
      linkobj <- make.link(linkstr)
      linkobj$d2mu.deta <- betareg:::make.d2mu.deta(linkstr)
    }
    if(linkstr == "loglog") {
      linkobj <- structure(list(linkfun = function(mu) -log(-log(mu)),
                                linkinv = function(eta) pmax(pmin(exp(-exp(-eta)),
                                                                  1 - .Machine$double.eps), .Machine$double.eps),
                                mu.eta = function(eta) {
                                  eta <- pmin(eta, 700)
                                  pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
                                }, d2mu.deta = function(eta) pmax(exp(-exp(-eta) -
                                                                        eta) * expm1(-eta), .Machine$double.eps), valideta = function(eta) TRUE,
                                name = "loglog"), class = "link-glm")
    }
    if(linkstr == "cobit") {
      linkobj <- structure(list(linkfun = bftprimeinv, # g
                                linkinv = bftprime, # ginv
                                mu.eta = bftprimeprime, #ginvprime
                                d2mu.deta = bftprimeprimeprime,
                                valideta = function(eta) TRUE,
                                name = "cobit"), class = "link-glm")
    }
  }
  else {
    linkobj <- link
    linkstr <- link$name
    if (type != "ML") {
      if (is.null(linkobj$dmu.deta) & is.null(linkobj$d2mu.deta)) {
        warning("link needs to provide d2mu.deta component for BC/BR")
      }
      else {
        if (is.null(linkobj$d2mu.deta))
          linkobj$d2mu.deta <- linkobj$dmu.deta
      }
    }
  }
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta
  d2mu.deta <- linkobj$d2mu.deta
  if (is.character(link.phi)) {
    phi_linkstr <- link.phi
    phi_linkobj <- make.link(phi_linkstr)
    phi_linkobj$d2mu.deta <- betareg:::make.d2mu.deta(phi_linkstr)
  }
  else {
    phi_linkobj <- link.phi
    phi_linkstr <- link.phi$name
    if (type != "ML") {
      if (is.null(phi_linkobj$dmu.deta) & is.null(phi_linkobj$d2mu.deta)) {
        warning("link.phi needs to provide d2mu.deta component for BC/BR")
      }
      else {
        if (is.null(phi_linkobj$d2mu.deta))
          phi_linkobj$d2mu.deta <- phi_linkobj$dmu.deta
      }
    }
  }
  phi_linkfun <- phi_linkobj$linkfun
  phi_linkinv <- phi_linkobj$linkinv
  phi_mu.eta <- phi_linkobj$mu.eta
  phi_d2mu.deta <- phi_linkobj$d2mu.deta
  ystar <- qlogis(y)
  ocontrol <- control
  phi_full <- control$phi
  method <- control$method
  gradient <- control$gradient
  hessian <- control$hessian
  start <- control$start
  fsmaxit <- control$fsmaxit
  fstol <- control$fstol
  quad <- control$quad
  control$phi <- control$method <- control$gradient <- control$hessian <- control$start <- control$fsmaxit <- control$fstol <- control$quad <- NULL
  if (is.null(gradient))
    gradient <- (dist == "beta")
  if (is.null(start)) {
    auxreg <- lm.wfit(x, if (dist == "beta")
      linkfun(y)
      else linkfun((y * (n - 1) + 0.5)/n), weights, offset = offset[[1L]])
    beta <- auxreg$coefficients
    yhat <- linkinv(auxreg$fitted.values)
    dlink <- 1/mu.eta(linkfun(yhat))
    res <- auxreg$residuals
    res[weights <= 0] <- 0
    sigma2 <- sum(weights * res^2)/((sum(weights) - k) *
                                      (dlink)^2)
    phi_y <- weights * yhat * (1 - yhat)/(sum(weights) *
                                            sigma2) - 1/n
    phi <- rep.int(0, ncol(z))
    phi[1L] <- suppressWarnings(phi_linkfun(sum(phi_y)))
    if (!isTRUE(phi_linkinv(phi[1L]) > 0)) {
      warning("no valid starting value for precision parameter found, using 1 instead")
      phi[1L] <- 1
    }
    start <- list(mu = beta, phi = phi)
    if (estnu)
      start$nu <- log(mean(y <= 0 | y >= 1))
  }
  if (is.list(start))
    start <- do.call("c", start)
  indices01 <- (y > 0) & (y < 1)
  indices0 <- (y <= 0)
  indices1 <- (y >= 1)
  fitfun <- function(par, deriv = 0L) {
    beta <- par[seq.int(length.out = k)]
    gamma <- par[seq.int(length.out = m) + k]
    nu <- if (estnu)
      exp(par[k + m + 1])
    else nu
    eta <- as.vector(x %*% beta + offset[[1L]])
    phi_eta <- as.vector(z %*% gamma + offset[[2L]])
    mu <- linkinv(eta)
    phi <- phi_linkinv(phi_eta)
    shape1 <- mu * phi
    shape2 <- (1 - mu) * phi
    if (deriv >= 1L) {
      d1 <- digamma(shape1)
      d2 <- digamma(shape2)
      mustar <- d1 - d2
    }
    else {
      d1 <- d2 <- mustar <- NULL
    }
    if (deriv >= 2L) {
      psi1 <- trigamma(shape1)
      psi2 <- trigamma(shape2)
    }
    else {
      psi1 <- psi2 <- NULL
    }
    if (deriv >= 3L) {
      b12 <- beta(shape1, shape2)
      d12 <- digamma(phi)
    }
    else {
      b12 <- d12 <- NULL
    }
    list(shape1 = shape1, shape2 = shape2, d1 = d1, d2 = d2,
         b12 = b12, d12 = d12, beta = beta, gamma = gamma,
         nu = nu, eta = eta, phi_eta = phi_eta, mu = mu, phi = phi,
         mustar = mustar, psi1 = psi1, psi2 = psi2)
  }
  if (dist == "beta") {
    loglikfun <- function(par, fit = NULL) {
      if (is.null(fit))
        fit <- fitfun(par)
      with(fit, {
        if (any(!is.finite(phi)) | any(shape1 > 1e+300) |
            any(shape2 > 1e+300)) {
          NaN
        }
        else {
          ll <- suppressWarnings(dbeta(y, shape1, shape2,
                                       log = TRUE))
          ll[weights <= 0] <- 0
          if (any(!is.finite(ll)))
            NaN
          else sum(weights * ll)
        }
      })
    }
    gradfun <- function(par, sum = TRUE, fit = NULL) {
      if (is.null(fit))
        fit <- fitfun(par, deriv = 1L)
      rval <- with(fit, {
        cbind(phi * (ystar - mustar) * mu.eta(eta) *
                weights * x, (mu * (ystar - mustar) + log(1 -
                                                            y) - digamma((1 - mu) * phi) + digamma(phi)) *
                phi_mu.eta(phi_eta) * weights * z)
      })
      rval[weights <= 0, ] <- 0
      if (sum)
        colSums(rval)
      else rval
    }
    hessfun <- function(par, inverse = FALSE, fit = NULL) {
      if (is.null(fit))
        fit <- fitfun(par, deriv = 2L)
      with(fit, {
        a <- psi1 + psi2
        b <- psi1 * mu^2 + psi2 * (1 - mu)^2 - trigamma(phi)
        wbb <- phi^2 * a * mu.eta(eta)^2
        wpp <- b * phi_mu.eta(phi_eta)^2
        wbp <- phi * (mu * a - psi2) * mu.eta(eta) *
          phi_mu.eta(phi_eta)
        kbb <- if (k > 0L)
          crossprod(sqrt(weights) * sqrt(wbb) * x)
        else crossprod(x)
        kpp <- if (m > 0L)
          crossprod(sqrt(weights) * sqrt(wpp) * z)
        else crossprod(z)
        kbp <- if (k > 0L & m > 0L)
          crossprod(weights * wbp * x, z)
        else crossprod(x, z)
        K <- cbind(rbind(kbb, t(kbp)), rbind(kbp, kpp))
        if (inverse)
          chol2inv(chol(K))
        else K
      })
    }
    biasfun <- function(par, fit = NULL, vcov = NULL) {
      if (is.null(fit))
        fit <- fitfun(par, deriv = 2L)
      InfoInv <- if (is.null(vcov))
        try(hessfun(par, inverse = TRUE), silent = TRUE)
      else vcov
      mu <- fit$mu
      with(fit, {
        kappa2 <- psi1 + psi2
        D1 <- mu.eta(eta)
        D2 <- phi_mu.eta(phi_eta)
        D1dash <- d2mu.deta(eta)
        D2dash <- phi_d2mu.deta(phi_eta)
        dpsi1 <- psigamma(mu * phi, 2)
        dpsi2 <- psigamma((1 - mu) * phi, 2)
        kappa3 <- dpsi1 - dpsi2
        psi3 <- psigamma(phi, 1)
        dpsi3 <- psigamma(phi, 2)
        PQsum <- function(t) {
          if (t <= k) {
            Xt <- x[, t]
            bb <- if (k > 0L)
              crossprod(x, weights * phi^2 * D1 * (phi *
                                                     D1^2 * kappa3 + D1dash * kappa2) * Xt *
                          x)
            else crossprod(x)
            bg <- if ((k > 0L) & (m > 0L))
              crossprod(x, weights * phi * D1^2 * D2 *
                          (mu * phi * kappa3 + phi * dpsi2 + kappa2) *
                          Xt * z)
            else crossprod(x, z)
            gg <- if (m > 0L)
              crossprod(z, weights * phi * D1 * D2^2 *
                          (mu^2 * kappa3 - dpsi2 + 2 * mu * dpsi2) *
                          Xt * z) + crossprod(z, weights * phi *
                                                D1 * D2dash * (mu * kappa2 - psi2) *
                                                Xt * z)
            else crossprod(z)
          }
          else {
            Zt <- z[, t - k]
            bb <- if (k > 0L)
              crossprod(x, weights * phi * D2 * (phi *
                                                   D1^2 * mu * kappa3 + phi * D1^2 * dpsi2 +
                                                   D1dash * mu * kappa2 - D1dash * psi2) *
                          Zt * x)
            else crossprod(x)
            bg <- if ((k > 0L) & (m > 0L))
              crossprod(x, weights * D1 * D2^2 * (phi *
                                                    mu^2 * kappa3 + phi * (2 * mu - 1) *
                                                    dpsi2 + mu * kappa2 - psi2) * Zt * z)
            else crossprod(x, z)
            gg <- if (m > 0L)
              crossprod(z, weights * D2^3 * (mu^3 * kappa3 +
                                               (3 * mu^2 - 3 * mu + 1) * dpsi2 - dpsi3) *
                          Zt * z) + crossprod(z, weights * D2dash *
                                                D2 * (mu^2 * kappa2 + (1 - 2 * mu) *
                                                        psi2 - psi3) * Zt * z)
            else crossprod(z)
          }
          pq <- rbind(cbind(bb, bg), cbind(t(bg), gg))
          sum(diag(InfoInv %*% pq))/2
        }
        if (inherits(InfoInv, "try-error")) {
          bias <- adjustment <- rep.int(NA_real_, k +
                                          m)
        }
        else {
          adjustment <- sapply(1:(k + m), PQsum)
          bias <- -InfoInv %*% adjustment
        }
        list(bias = bias, adjustment = adjustment)
      })
    }
  }
  else {
    dfun <- if (dist == "xbeta") {
      function(x, mu, phi, nu, ...) dxbeta(x, mu = mu,
                                           phi = phi, nu = nu, ...)
    }
    else {
      quadrule <- betareg:::quadtable(nquad = quad)
      function(x, mu, phi, nu, ...) dxbetax(x, mu = mu,
                                            phi = phi, nu = nu, quad = quadrule, ...)
    }
    loglikfun <- function(par, fit = NULL) {
      if (is.null(fit))
        fit <- fitfun(par)
      with(fit, {
        if (any(!is.finite(phi)) | any(phi < 0))
          NaN
        else {
          rval <- dfun(y, mu = mu, phi = phi, nu = nu,
                       log = TRUE)
          sum(weights * rval)
        }
      })
    }
    gradfun_xbeta <- function(par, sum = TRUE, fit = NULL) {
      if (is.null(fit))
        fit <- fitfun(par, deriv = 3L)
      with(fit, {
        ynu <- (y + nu)/(1 + 2 * nu)
        ystarnu <- qlogis(ynu)
        grad_l01 <- cbind(phi * (ystarnu - mustar) *
                            mu.eta(eta) * weights * x, (mu * (ystarnu -
                                                                mustar) + log(1 - ynu) - digamma((1 - mu) *
                                                                                                   phi) + digamma(phi)) * phi_mu.eta(phi_eta) *
                            weights * z, ((mu * phi - 1)/(y + nu) + ((1 -
                                                                        mu) * phi - 1)/(1 - y + nu) - 2 * (phi - 1)/(1 +
                                                                                                                       2 * nu)) * weights * nu)
        nu_low <- nu/(1 + 2 * nu)
        nu_upp <- (1 + nu)/(1 + 2 * nu)
        dlow <- dbeta(nu_low, shape1, shape2)
        plow <- pbeta(nu_low, shape1, shape2)
        dupp <- dbeta(nu_upp, shape1, shape2)
        pupp <- pbeta(nu_upp, shape1, shape2)
        Fs1 <- h3f2(shape1, shape2, nu_low, n, maxiter = 10000,
                    eps = 0)
        Fs2 <- h3f2(shape1, shape2, nu_upp, n, maxiter = 10000,
                    eps = 0)
        Fs3 <- h3f2(shape2, shape1, nu_low, n, maxiter = 10000,
                    eps = 0)
        Fs4 <- h3f2(shape2, shape1, nu_upp, n, maxiter = 10000,
                    eps = 0)
        delta1low <- plow * (d12 - d1 + log(nu_low)) -
          nu_low^shape1 * Fs1/(shape1^2 * b12)
        delta2low <- (1 - plow) * (d2 - d12 - log(nu_upp)) +
          nu_upp^shape2 * Fs4/(shape2^2 * b12)
        delta1upp <- pupp * (d12 - d1 + log(nu_upp)) -
          nu_upp^shape1 * Fs2/(shape1^2 * b12)
        delta2upp <- (1 - pupp) * (d2 - d12 - log(nu_low)) +
          nu_low^shape2 * Fs3/(shape2^2 * b12)
        grad_l0 <- cbind(phi * mu.eta(eta) * (delta1low -
                                                delta2low) * weights * x/plow, phi_mu.eta(phi_eta) *
                           (delta1low * mu + delta2low * (1 - mu)) * weights *
                           z/plow, dlow/(plow * (1 + 2 * nu)^2) * weights *
                           nu)
        grad_l1 <- cbind(phi * mu.eta(eta) * (delta2upp -
                                                delta1upp) * weights * x/(1 - pupp), -phi_mu.eta(phi_eta) *
                           (delta1upp * mu + delta2upp * (1 - mu)) * weights *
                           z/(1 - pupp), dupp/((1 - pupp) * (1 + 2 * nu)^2) *
                           weights * nu)
        grad_l0[!indices0, ] <- 0
        grad_l1[!indices1, ] <- 0
        grad_l01[!indices01, ] <- 0
        out <- grad_l0 + grad_l1 + grad_l01
        out <- if (estnu)
          out
        else out[, 1:(k + m)]
        if (sum) {
          colSums(out)
        }
        else {
          out
        }
      })
    }
    if (dist == "xbeta") {
      gradfun <- gradfun_xbeta
    }
    if (dist == "xbetax") {
      gradfun <- function(par, sum = TRUE, fit = NULL) {
        if (is.null(fit)) {
          fit <- fitfun(par, deriv = 3L)
        }
        with(fit, {
          dens <- apply(quadrule, 1, function(rule) {
            e <- rule[1] * nu
            rule[2] * dxbeta(y, mu, phi, nu = e, log = FALSE)
          })
          tdens <- rowSums(dens)
          obsders <- lapply(seq.int(quad), function(inds) {
            current_nu <- fit$nu <- quadrule[inds, 1] *
              nu
            par[k + m + 1] <- log(current_nu)
            out <- gradfun_xbeta(par, fit = fit, sum = FALSE)
            dens[, inds] * out/tdens
          })
          out <- Reduce("+", obsders)
          if (sum) {
            colSums(out)
          }
          else {
            out
          }
        })
      }
    }
  }
  if (method == "nlminb") {
    stopifnot(requireNamespace("numDeriv"))
    if ("maxit" %in% control) {
      if (is.null(control$iter.max))
        control$iter.max <- control$maxit
      control$maxit <- NULL
    }
    if ("reltol" %in% control) {
      if (is.null(control$rel.tol))
        control$rel.tol <- control$reltol
      control$reltol <- NULL
    }
    opt <- nlminb(start = start, objective = function(par,
                                                      ...) -loglikfun(par, ...), gradient = if (gradient)
                                                        function(par, ...) -gradfun(par, ...)
                  else NULL, control = control)
    opt$hessian <- numDeriv::hessian(loglikfun, opt$par)
  }
  else {
    opt <- optim(par = start, fn = loglikfun, gr = if (gradient)
      gradfun
      else NULL, method = method, hessian = hessian, control = control)
  }
  par <- opt$par
  if (type == "BR" & fsmaxit <= 0)
    warning("BR cannot be performed with fsmaxit <= 0")
  step <- .Machine$integer.max
  iter <- 0
  if (fsmaxit > 0 & !(hessian & type == "ML")) {
    for (iter in 1:fsmaxit) {
      stepPrev <- step
      stepFactor <- 0
      testhalf <- TRUE
      while (testhalf & stepFactor < 11) {
        fit <- fitfun(par, deriv = 2L)
        scores <- gradfun(par, fit = fit)
        InfoInv <- try(hessfun(par, fit = fit, inverse = TRUE))
        if (failedInv <- inherits(InfoInv, "try-error")) {
          warning("failed to invert the information matrix: iteration stopped prematurely")
          break
        }
        bias <- if (type == "BR")
          biasfun(par, fit = fit, vcov = InfoInv)$bias
        else 0
        par <- par + 2^(-stepFactor) * (step <- InfoInv %*%
                                          scores - bias)
        stepFactor <- stepFactor + 1
        testhalf <- drop(crossprod(stepPrev) < crossprod(step))
      }
      if (failedInv | (all(abs(step) < fstol))) {
        break
      }
    }
  }
  if ((fsmaxit == 0 & opt$convergence > 0) | (iter >= fsmaxit &
                                              fsmaxit > 0)) {
    converged <- FALSE
    warning("optimization failed to converge")
  }
  else {
    converged <- TRUE
  }
  if (type == "BC") {
    bias <- as.vector(biasfun(par)$bias)
    par <- par - bias
  }
  else {
    bias <- rep.int(NA_real_, k + m + estnu)
  }
  fit <- fitfun(par, deriv = 3L)
  beta <- fit$beta
  gamma <- fit$gamma
  eta <- fit$eta
  mu <- fit$mu
  phi <- fit$phi
  nu <- if (!estnu)
    nu
  else as.vector(exp(par[k + m + 1]))
  ll <- loglikfun(par, fit = fit)
  if (gradient) {
    ef <- gradfun(par, fit = fit, sum = FALSE)
  }
  else {
    stopifnot(requireNamespace("numDeriv"))
    ef <- numDeriv::grad(loglikfun, par)
  }
  vcov <- if (hessian && (type == "ML"))
    solve(-as.matrix(opt$hessian))
  else hessfun(fit = fit, inverse = TRUE)
  wcor <- function(x, y, weights = NULL) {
    if (is.null(weights) || identical(rep.int(1, length(x)),
                                      weights))
      return(cor(x, y))
    x <- x[weights > 0]
    y <- y[weights > 0]
    w <- weights[weights > 0]/sum(weights)
    x <- x - sum(x * w)
    x <- x/sqrt(sum(w * x^2))
    y <- y - sum(y * w)
    y <- y/sqrt(sum(w * y^2))
    sum(w * x * y)
  }
  pseudor2 <- if (dist != "beta" || var(eta[weights > 0]) *
                  var(ystar[weights > 0]) <= 0)
    NA
  else wcor(eta, linkfun(y), weights)^2
  names(beta) <- colnames(x)
  names(gamma) <- if (phi_const & phi_linkstr == "identity")
    "(phi)"
  else colnames(z)
  rownames(vcov) <- colnames(vcov) <- names(bias) <- c(colnames(x),
                                                       if (phi_const & phi_linkstr == "identity") "(phi)" else paste("(phi)",
                                                                                                                     colnames(z), sep = "_"), if (estnu) "Log(nu)" else NULL)
  marg_e <- switch(dist, beta = mu, xbeta = vapply(seq.int(n),
                                                   function(i) mean_xbeta(mu[i], phi[i], nu), 0), xbetax = vapply(seq.int(n),
                                                                                                                  function(i) mean_xbetax(mu[i], phi[i], nu, quad), 0))
  rval <- list(coefficients = list(mu = beta, phi = gamma),
               residuals = y - marg_e, fitted.values = structure(marg_e,
                                                                 .Names = names(y)), type = type, dist = dist, optim = opt,
               method = method, control = ocontrol, scoring = iter,
               start = start, weights = if (identical(as.vector(weights),
                                                      rep.int(1, n))) NULL else weights, offset = list(mu = if (identical(offset[[1L]],
                                                                                                                          rep.int(0, n))) NULL else offset[[1L]], phi = if (identical(offset[[2L]],
                                                                                                                                                                                      rep.int(0, n))) NULL else offset[[2L]]), n = n, nobs = nobs,
               df.null = nobs - 2 - estnu, df.residual = nobs - k -
                 m - estnu, phi = phi_full, nu = nu, loglik = ll,
               vcov = vcov, bias = bias, pseudo.r.squared = pseudor2,
               link = list(mu = linkobj, phi = phi_linkobj), converged = converged,
               grad = ef)
  if (estnu)
    rval$coefficients$nu <- c(`Log(nu)` = log(nu))
  if (dist == "beta") {
    for (n in intersect(names(rval), betareg:::fix_names_mu_phi)) names(rval[[n]])[1L:2L] <- c("mean",
                                                                                     "precision")
  }
  return(rval)
}
