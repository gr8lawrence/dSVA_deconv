## update dSVA extension for this setting
## what if we assume we know q = 1
## allow an intercept

dsva_for_sim <- function(Y, Theta, n_comp = 0, alg = c("nnls", "pnnls"), intercept = TRUE, test = FALSE) {
  n <- ncol(Y)
  m <- nrow(Y)
  
  ## add an intercept if chosen
  if (intercept) {
    X <- model.matrix(~1 + Theta)
  } else {
    X <- Theta
  }
  
  # if no hidden factor is estimated, return the P_hat from the simplified model
  # if (q == 0) {
  #   message("Use PA to estimate the number of hidden factors (under development).")
  # }
  if (n_comp == 0) {
    message("Specifying number of latent factors = 0.")
    P_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, sum = 1)$x})
    return(P_hat)
  }
  
  ## step 1: obtain the canonical model residual (using linear regression here instead of PNNLS)
  svd_X <- svd(X)
  U_x <- svd_X$u
  Sigma_x <- svd_X$d
  V_x <- svd_X$v
  B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y   
  R <- Y - X %*% B_star_hat
  
  ## step 2: svd on the residual space
  svd_R <- svd(R)
  if (n_comp == 1) {
    U_q <- as.matrix(svd_R$u[, 1]) 
    Psi_hat <- rbind(svd_R$d[1] %*% t(svd_R$v)[1, ])
  } else {
    U_q <- svd_R$u[, seq(n_comp)]
    Psi_hat <- diag(svd_R$d[seq(n_comp)]) %*% t(svd_R$v)[seq(n_comp), ]
  }
  
  ## step 3: estimate the surrogate variable
  # Psi_hat <- t(U_q) %*% R
  M_jn <- matrix(1/n, n, n)
  D_jn <- diag(1, n) - M_jn
  if (n_comp == 1) {
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

  ## step 4: fitting the model again with the surrogate variable (using NNLS or PNNLS)
  X_new <- cbind(Gamma_hat, X)
  
  if (alg == "pnnls") {
    B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = ncol(X_new) - ncol(Theta), sum = 1)$x})
    P_hat <- B_hat[-seq(ncol(X_new) - ncol(Theta)), ]
  } else if (alg == "nnls") {
    B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = ncol(X_new) - ncol(Theta), sum = NULL)$x})
    P_tilde <- B_hat[-seq(ncol(X_new) - ncol(Theta)), ]
    P_hat <- apply(P_tilde, 2, function(x) x/sum(x))
  }
  
  ## return the P_hat
  if (!test) {
    return(P_hat)
  } else {
    Y_lat_hat <- Y - Theta %*% P_hat
    return(list(P_hat = P_hat,
                Y_lat_hat = Y_lat_hat))
  }
}


## a first pass method to find the number of components
estimate_n_comp <- function(Y, Theta, intercept = TRUE) {
  
  ## add an intercept 
  if (intercept) {
    X <- model.matrix(~1 + Theta)
  } else {
    X <- Theta
  }
  
  ## perform first pass regression
  Bhat <- apply(Y, 2, function(y) coef(lm(y ~ X - 1)))
  R <- Y - X %*% Bhat
  
  ## use the gap of log singular values
  d_R <- log(prcomp(t(R))$sdev^2)
  diff <- d_R[seq(length(d_R) - 2)] - d_R[seq(2, length(d_R) - 1)]
  n_comp <- which(diff == max(diff))
  n_comp
}


## code for tests
# Y <- true_data$Y
# Theta <- true_data$X
# intercept <- TRUE

# estimate_n_comp(true_data$Y, true_data$X)
