## update dSVA extension for this setting
## what if we assume we know q = 1
## allow an intercept
dsva_for_sim <- function(Y, Theta, n_comp = 0, intercept = TRUE) {
  n <- ncol(Y)
  m <- nrow(Y)
  q <- n_comp # rename n_comp to q
  
  ## add an intercept if chosen
  if (intercept) {
    X <- model.matrix(~1 + Theta)
  } else {
    X <- Theta
  }
  
  # if no hidden factor is estimated, return the P_hat from the simplified model
  if (q == 0) {
    message("Use PA to estimate the number of hidden factors (under development).")
  }
  
  ## step 1: obtain the canonical model residual (using linear regression here instead of PNNLS)
  svd_X <- svd(X)
  U_x <- svd_X$u
  Sigma_x <- svd_X$d
  V_x <- svd_X$v
  B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y   
  R <- Y - X %*% B_star_hat
  
  ## step 2: svd on the residual space
  # q = 2
  svd_R <- svd(R)
  if (q == 1) {
    U_q <- as.matrix(svd_R$u[, 1]) 
    Psi_hat <- rbind(svd_R$d[1] %*% t(svd_R$v)[1, ])
  } else {
    U_q <- svd_R$u[, seq(q)]
    Psi_hat <- diag(svd_R$d[seq(q)]) %*% t(svd_R$v)[seq(q), ]
  }
  
  ## step 3: estimate the surrogate variable
  # Psi_hat <- t(U_q) %*% R
  M_jn <- matrix(1/n, n, n)
  D_jn <- diag(1, n) - M_jn
  
  ## TODO: see Seunggeung's code to simplify this (you can reduce some computation)
  if (q == 1) {
    C <- 1/(Psi_hat %*% D_jn %*% t(Psi_hat))
    if (C[1, 1] < 1e-15) message("q_hat = 1 and Psi_hat %*% D_jn %*% t(Psi_hat) is smaller than 1e-15.")
    Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) * C[1, 1]
  } else {
    C <- Psi_hat %*% D_jn %*% t(Psi_hat)
    if (1/kappa(C) < 1e-15) message("q_hat > 1 and Psi_hat %*% D_jn %*% t(Psi_hat) is close to singular, i.e. its reverse condition number < 1e-15.")
    ## use svd to invert C to avoid invertibility problems
    svd_C <- svd(C)
    U_C <- svd_C$u
    Sigma_C <- svd_C$d
    V_C <- svd_C$v
    C_inv <- svd_C$u %*% diag(1/Sigma_C) %*% t(svd_C$v)
    Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) %*% C_inv
  }

  ## step 4: fitting the model again with the surrogate variable (this time using PNNLS)
  X_new <- cbind(Gamma_hat, X)
  B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = ncol(X_new) - ncol(Theta), sum = 1)$x})
  P_hat <- B_hat[-seq(ncol(X_new) - ncol(Theta)), ]

  ## return the P_hat
  P_hat
}


