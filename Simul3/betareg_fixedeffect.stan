// generated with brms 2.22.0  (edited: no random effects)
functions {
  /* the inverse cobit function */
  vector inv_cobit(vector x) {
    return 1/(1-exp(-x))-1/x;
  }
  real unifsq_lpdf(real x){
    return -log(2*8.74) - 0.5*log(x);
  }
}
data {
  int<lower=1> N;           // total number of observations
  vector[N] Y;              // response variable in (0,1)
  int<lower=1> K;           // number of population-level effects
  matrix[N, K] X;           // population-level design matrix
  int<lower=1> Kc;          // number of population-level effects after centering
  // removed: coords, phi_spatial (GP)
  int prior_only;           // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;         // centered version of X without an intercept
  vector[Kc] means_X;       // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  real Intercept;           // temporary intercept for centered predictors
  vector[Kc] b;             // regression coefficients
  real<lower=0, upper=76.3876> phi;        // precision parameter
  // removed: sd_1, z_1 (group effects), u (GP latent), sigma (GP scale)
}
transformed parameters {
  real lprior = 0;          // prior contributions to the log posterior
  lprior += normal_lpdf(b | 0, 100);
  lprior += normal_lpdf(Intercept | 0, 100);
  lprior += unifsq_lpdf(phi);
  // removed: cauchy_lpdf(sigma | 0, 1)
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Xc * b;

    // removed random effects

    mu = inv_cobit(mu);
    target += beta_lpdf(Y | mu * phi, (1 - mu) * phi);
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
