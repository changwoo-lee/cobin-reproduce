
predict_linpred_cobit <- function(coords_train, u_save, beta_save, sigmasq_save,
                            coords_test, X_test, phi_spatial) {
  nsave <- nrow(beta_save)
  eta_test_save = matrix(0, nrow = nsave, ncol = nrow(X_test))
  D_nn <- fields::rdist(coords_train)
  D_nt <- fields::rdist(coords_train, coords_test)
  D_tt <- fields::rdist(coords_test)
  # Progress bar
  Rnn <- fields::Matern(D_nn, range = 1/phi_spatial, smoothness = 0.5)
  Rnt <- fields::Matern(D_nt, range = 1/phi_spatial, smoothness = 0.5)
  Rtt <- fields::Matern(D_tt, range = 1/phi_spatial, smoothness = 0.5)
  Rnninv_Rnt <- solve(Rnn, Rnt)

  pb <- txtProgressBar(style = 3)
  for (isave in 1:nsave) {
    # Check progress
    setTxtProgressBar(pb, isave / (nsave))
    # Calculate kernels

    # Prediction
    u_train <- as.numeric(u_save[isave,])
    temp_mean <- t(Rnninv_Rnt) %*% u_train
    temp_var <- sigmasq_save[isave] * (Rtt - t(Rnninv_Rnt) %*% Rnt)
    u_test = as.numeric(mvnfast::rmvn(1, mu = temp_mean, sigma = temp_var))
    eta_test_save[isave,] = X_test%*%beta_save[isave,] + u_test
  }
  return(list(linpred_test_save = eta_test_save,
              mu_test_save = cobin::bftprime(eta_test_save)))
}
