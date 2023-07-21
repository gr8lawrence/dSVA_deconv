## update dSVA extension for this setting
## what if we assume we know q = 1
dsva_ext_know_q <- function(Y, X, q = 1) {
  n <- ncol(Y)
  m <- nrow(Y)
  
  # Y <- true_data$Y
  # X <- true_data$X
  # q = 1
  ## step 1: obtain the canonical model residual (using linear regression here instead of PNNLS)
  svd_X <- svd(X)
  U_x <- svd_X$u
  Sigma_x <- svd_X$d
  V_x <- svd_X$v
  M_x <- U_x %*% t(U_x)
  B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y   
  R <- Y - X %*% B_star_hat

  # if no hidden factor is estimated, return the P_hat from the simplified model
  if (q == 0) {
    message("Assume no latent variables.")
    return(B_star_hat)
  }
  
  ## step 2: svd on the residual space
  if (q == 1) {
    U_q <- as.matrix(svd(R)$u[, seq(q)]) 
  } else {
    U_q <- svd(R)$u[, seq(q)]
  }
  
  ## step 3: estimate the surrogate variable
  # Psi_hat <- apply(R, 2, function(r) {lsei::lsei(a = U_q, b = r)})
  Psi_hat <- t(U_q) %*% Y 
  M_jn <- matrix(1/n, n, n)
  D_jn <- diag(1, n) - M_jn
  
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
    Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) %*% svd_C$v
  }
  
  # if (q == 1) {
  #   if (Psi_hat %*% D_jn %*% as.matrix(Psi_hat) < 1e-15) {
  #     message("q_hat = 1 and Psi_hat %*% D_jn %*% as.matrix(Psi_hat) is singular, Gamma_hat is assigned to equal U_q.")
  #     Gamma_hat <- U_q
  #   } else {
  #     Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% as.matrix(Psi_hat) %*% (Psi_hat %*% D_jn %*% as.matrix(Psi_hat))^(-1)  
  #   }
  # } else {
  #   if (kappa(Psi_hat %*% D_jn %*% t(Psi_hat)) < 1e-15) {
  #     message("q_hat > 1 and Psi_hat %*% D_jn %*% t(Psi_hat) is singular, Gamma_hat is assigned to equal U_q.")
  #     Gamma_hat <- U_q
  #   } else {
  #     Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) %*% solve(Psi_hat %*% D_jn %*% t(Psi_hat))
  #   }
  # }

  ## step 4: fitting the model again with the surrogate variable (this time using PNNLS)
  X_new <- cbind(Gamma_hat, X)
  B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = q, sum = 1)$x})
  P_hat <- B_hat[-seq(q), ]

  ## return the P_hat
  return(list(P_hat = P_hat,
              q_hat = q))
}

## update dSVA extension for this setting
dsva_ext <- function(Y, X, ...) {
  n <- ncol(Y)
  m <- nrow(Y)
  
  # Y <- true_data$Y
  # X <- true_data$X
  
  ## step 1: obtain the canonical model residual (using linear regression here instead of PNNLS)
  svd_X <- svd(X)
  U_x <- svd_X$u
  Sigma_x <- svd_X$d
  V_x <- svd_X$v
  M_x <- U_x %*% t(U_x)
  B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y   
  R <- Y - X %*% B_star_hat
  
  # Y <- true_data$Y
  # X <- true_data$X
  # however, using this R to estimate q will lead to a lot of biases (won't estimate correctly most of the times)
  # estimate the pNNLS residual instead
  # B2 <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, sum = 1)$x})
  # # we project the residual to be orthogonal to X
  # R2 <- (diag(1, m) - M_x) %*% (Y - X %*% B2) # R2 is actually best
  # R3 <- Y - X %*% B2
  # however, using this R to estimate q will lead to a lot of biases (won't estimate correctly most of the times)
  
  # now use Horn's PA to estimate the number of factors q in the residual
  q_ls <- paran::paran(x = R, ...)
  q <- q_ls$Retained 
  # R3 <- Y - X %*% B2 
  # if no hidden factor is estimated, return the P_hat from the simplified model
  if (q == 0) {
    message("Assume no latent variables.")
    return(B_star_hat)
  }
  
  ## step 2: svd on the residual space
  if (q == 1) {
    U_q <- as.matrix(svd(R)$u[, seq(q)]) 
  } else {
    U_q <- svd(R)$u[, seq(q)]
  }
  
  ## step 3: estimate the surrogate variable
  # Psi_hat <- apply(R, 2, function(r) {lsei::lsei(a = U_q, b = r)})
  Psi_hat <- t(U_q) %*% Y 
  M_jn <- matrix(1/n, n, n)
  D_jn <- diag(1, n) - M_jn
  
  if (q == 1) {
    # C <- 1/(Psi_hat %*% D_jn %*% as.matrix(Psi_hat))
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
    Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) %*% svd_C$v
  }
  
  ## step 4: fitting the model again with the surrogate variable (this time using PNNLS)
  X_new <- cbind(Gamma_hat, X)
  B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = q, sum = 1)$x})
  P_hat <- B_hat[-seq(q), ]
  
  ## return the P_hat
  return(list(P_hat = P_hat,
              q_hat = q))
}


## estimate the eigenvalues

