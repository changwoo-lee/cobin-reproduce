// generated with brms 2.22.0
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
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  array[N] vector[2] coords;
  real phi_spatial;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> phi;  // precision parameter
  // vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  // array[M_1] vector[N_1] z_1;  // standardized group-level effects
  vector[N] u;
  real<lower = 0> sigma;
  //real<lower = 0.09, upper = 0.11> phi_spatial;
}
transformed parameters {
  real sigmasq = square(sigma);
  //vector[N * (N - 1) / 2] vdist = get_vdist(coords);
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(b | 0, 100);
  lprior += normal_lpdf(Intercept | 0, 100);
  lprior += unifsq_lpdf(phi);
  lprior += cauchy_lpdf(sigma | 0, 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Xc * b;

    matrix[N,N] KernelL = cholesky_decompose(gp_exponential_cov(coords, sigma, 1/phi_spatial));
    //K = add_diag(K, 1e-9); // for numerical stability
    u ~ multi_normal_cholesky(rep_vector(0.0, N), KernelL);
    mu += u;
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
