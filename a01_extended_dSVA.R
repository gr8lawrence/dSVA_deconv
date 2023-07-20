## write the algorithm (requires the 'lsei' package)

# ## generate a model for testing
# m = 15
# n = 30
# K = 3
# p = 2
# X <- matrix(rchisq(m * K, df = 10), m, K)
# B <- matrix(rchisq(K * n, df = 5), K, n)
# B <- apply(B, 2, function(b) b/sum(b))
# Z <- matrix(rnorm(m * p, 2), m, p)
# D <- matrix(rnorm(p * n), p, n)
# E <- matrix(rnorm(m * n, 0, 0.01), m, n)
# Y <- X %*% B + Z %*% D + E
# sum(Y < 0) # make sure this is 0
# Y[Y < 0] <- 0
# 
# ## assume we know q
# q <- ncol(Z)

# B_star_hat2 <- apply(Y, 2, function(y) {lsei::lsei(a = X, b = y)})
# R2 <- Y - X %*% B_star_hat2 
# t(X) %*% R2 # roughly 0

## Y: bulk expression matrix
## X: signature matrix
## q: number of latent variables for adjustment
## solving P with sum-to-one constraint
dsva_ext <- function(Y, X, q) {
  n <- ncol(Y)
  m <- nrow(Y)
  
  ## step 1: obtain the canonical model residual 
  B_star_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, sum = 1)$x})
  # B_star_hat <- apply(Y, 2, function(y) {lsei::nnls(a = X, b = y)$x})
  # M_x <- X %*% solve(t(X) %*% X) %*% t(X)
  U_x <- svd(X)$u
  M_x <- U_x %*% t(U_x) 
  # R <- Y - X %*% B_star_hat # challenge No.1: the residual space isn't orthogonal to columns of X
  R <- (diag(1, m) - M_x) %*% (Y - X %*% B_star_hat) # we project the residual to be orthogonal to X
  
  ## step 2: svd on the residual space
  U_q <- svd(R)$u[, seq(q)] 
  
  ## step 3: estimate the surrogate variable
  Psi_hat <- apply(R, 2, function(r) {lsei::lsei(a = U_q, b = r)})
  # J <- rep(1, n)
  # M_jn <- J %*% solve(t(J) %*% J) %*% t(J) 
  M_jn <- matrix(1/n, n, n)
  D_jn <- diag(1, n) - M_jn
  Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) %*% solve(Psi_hat %*% D_jn %*% t(Psi_hat))  
  
  ## step 4: fitting the model again with the surrogate variable
  X_new <- cbind(Gamma_hat, X)
  B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = q, sum = 1)$x})
  P_hat <- B_hat[-seq(q), ]
  P_hat
}

## solving P without sum-to-one constraint
dsva_ext_nn <- function(Y, X, q) {
  n <- ncol(Y)
  m <- nrow(Y)
  
  ## step 1: obtain the canonical model residual 
  B_star_hat <- apply(Y, 2, function(y) {lsei::nnls(a = X, b = y)$x})
  U_x <- svd(X)$u
  M_x <- U_x %*% t(U_x) 
  R <- (diag(1, m) - M_x) %*% (Y - X %*% B_star_hat) # we project the residual to be orthogonal to X
  
  ## step 2: svd on the residual space
  U_q <- svd(R)$u[, seq(q)] 
  
  ## step 3: estimate the surrogate variable
  Psi_hat <- apply(R, 2, function(r) {lsei::lsei(a = U_q, b = r)})
  M_jn <- matrix(1/n, n, n)
  D_jn <- diag(1, n) - M_jn
  Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) %*% solve(Psi_hat %*% D_jn %*% t(Psi_hat))  
  
  ## step 4: fitting the model again with the surrogate variable
  X_new <- cbind(Gamma_hat, X)
  B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = q, sum = NULL)$x})
  P_tilde <- B_hat[-seq(q), ]
  P_hat <- apply(P_tilde, 2, function(x) x/sum(x))
  P_hat
}

# mean(abs(B - P_hat))
# mean(abs(B - B_star_hat))

## Use dsva to impute cell-type specific expression adjusting for the hidden covariates
## X is the proportion matrix: each row represents a samples
dsva_ext2 <- function(Y, X, q) {
  n <- nrow(Y)
  m <- ncol(Y)
  
  ## step 1: obtain the canonical model residual 
  B_star_hat <- apply(Y, 2, function(y) {lsei::nnls(a = X, b = y)$x})
  # B_star_hat <- apply(Y, 2, function(y) {lsei::nnls(a = X, b = y)$x})
  # U_x <- svd(X)$u
  # M_x <- U_x %*% t(U_x) 
  # R <- Y - X %*% B_star_hat # when B_star is not forced to have sum to 1 constraint, the residual is orthogonal
  R <- (diag(1, n) - M_x) %*% (Y - X %*% B_star_hat) # we project the residual to be orthogonal to X
  
  ## step 2: svd on the residual space
  U_q <- svd(R)$u[, seq(q)] 
  
  ## step 3: estimate the surrogate variable
  Psi_hat <- apply(R, 2, function(r) {lsei::lsei(a = U_q, b = r)})
  # J <- rep(1, n)
  # M_jn <- J %*% solve(t(J) %*% J) %*% t(J) 
  M_jn <- matrix(1/m, m, m)
  D_jn <- diag(1, m) - M_jn
  Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) %*% solve(Psi_hat %*% D_jn %*% t(Psi_hat))  
  
  ## step 4: fitting the model again with the surrogate variable
  X_new <- cbind(Gamma_hat, X)
  B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = q, sum = NULL)$x})
  Theta_hat <- B_hat[-seq(q), ]
  Theta_hat
}

## alternatively solving the regession problem Y = XB + ZD + E where X and D are known
solve_reg <- function(Y, X, D, tol = 1e-5, max_iter = 1e4) {
  
  ## test data
  true_data <- dSVA_model_sim(m, n, K, p, p_sig, lambda)
  Y <- true_data$Y
  X <- true_data$X
  D <- true_data$D
  tol <- 1e-5
  max_iter <- 1e4
  
  ## initialize a B matrix of proportions and an empty Z matrix
  B <- t(extraDistr::rdirichlet(n = ncol(Y), alpha = rep(1, ncol(X))))
  Z <- matrix(nrow = nrow(Y), ncol = nrow(D))
  Y1 <- Y - X %*% B
  for (j in 1:nrow(Z)) Z[j, ] <- c(lsei::lsei(a = t(D), b = Y1[j, ]))
  
  ## record the estimated Y
  Y_tilde <- X %*% B + Z %*% D 
  
  t <- 0
  if (t <= max_iter) {
    if (t == max_iter) {
      warning("Maximal number of iteratinos reached.") 
      break
      }
    t <- t + 1
    ## then solve B with Z
    Y2 <- Y - Z %*% D
    for (i in 1:ncol(B)) B[, i] <- lsei::pnnls(a = X, b = Y2[, i], sum = 1)$x
    
    ## solve Z with B
    Y1 <- Y - X %*% B
    for (j in 1:nrow(Z)) Z[j, ] <- c(lsei::lsei(a = t(D), b = Y1[j, ]))
    
    ## calculate the tolerance
    Y_hat <- X %*% B + Z %*% D 
    d <- norm(Y_hat - Y_tilde, "F")^2/(nrow(Y) * ncol(Y))
    if (d < tol) break
    Y_tilde <- Y_hat
  }
  B - true_data$P_star
  return(list(B = B, Z = Z))
}

## update dSVA extension for this setting
dsva_ext2 <- function(Y, X, ...) {
  n <- ncol(Y)
  m <- nrow(Y)
  
  ## step 1: obtain the canonical model residual 
  B_star_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, sum = 1)$x})
  U_x <- svd(X)$u
  M_x <- U_x %*% t(U_x) 
  # we project the residual to be orthogonal to X
  R <- (diag(1, m) - M_x) %*% (Y - X %*% B_star_hat) 
  
  # estimate q, the number of hidden factors on R, using Horn's method
  q_ls <- paran::paran(x = R, ...)
  q <- q_ls$Retained
  
  # if no hidden factor is estimated, return the P_hat from the simplified model
  if (q == 0) {
    message("Horn's method detected no latent variables.")
    return(B_star_hat)
  }
  
  ## step 2: svd on the residual space
  if (q == 1) {
    U_q <- as.matrix(svd(R)$u[, seq(q)]) 
  } else {
    U_q <- svd(R)$u[, seq(q)]
  }
  
  ## step 3: estimate the surrogate variable
  Psi_hat <- apply(R, 2, function(r) {lsei::lsei(a = U_q, b = r)})
  M_jn <- matrix(1/n, n, n)
  D_jn <- diag(1, n) - M_jn
  
  if (q == 1) {
    if (Psi_hat %*% D_jn %*% as.matrix(Psi_hat) > 1e-15) {
      message("q_hat = 1 and Psi_hat %*% D_jn %*% as.matrix(Psi_hat) is singular, Gamma_hat is assigned to equal U_q.")
      Gamma_hat <- U_q
    } else {
      Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% as.matrix(Psi_hat) %*% solve(Psi_hat %*% D_jn %*% as.matrix(Psi_hat))  
    }
  } else {
    if (kappa(Psi_hat %*% D_jn %*% t(Psi_hat)) > 1e-15) {
      message("q_hat > 1 and Psi_hat %*% D_jn %*% t(Psi_hat) is singular, Gamma_hat is assigned to equal U_q.")
      Gamma_hat <- U_q
    } else {
      Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) %*% solve(Psi_hat %*% D_jn %*% t(Psi_hat))
    }
  }
  # Gamma_hat 
  
  ## step 4: fitting the model again with the surrogate variable
  X_new <- cbind(Gamma_hat, X)
  B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = q, sum = 1)$x})
  P_hat <- B_hat[-seq(q), ]
  
  ## return the P_hat
  return(list(P_hat = P_hat,
              q_hat = q))
}

## update dSVA extension for this setting
## what if we assume we know q = 1
dsva_ext2_know_q <- function(Y, X, q = 1) {
  
  n <- ncol(Y)
  m <- nrow(Y)
  
  # true_data <- dSVA_model_sim_2(m, n, K, p, p_sig, lambda, gamma)
  # Y <- true_data$Y
  # X <- true_data$X
  # Z <- true_data$Z
  # D <- true_data$D

  ## step 1: obtain the canonical model residual 
  # B_star_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, sum = 1)$x})
  
  svd_X <- svd(X)
  U_x <- svd_X$u
  Sigma_x <- svd_X$d
  V_x <- svd_X$v
  M_x <- U_x %*% t(U_x)
  B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y   
  # B_star_hat <- apply(Y, 2, function(y) {lsei::lsei(a = X, b = y)})
  
  # U_x <- svd(X)$u
  # M_x <- U_x %*% t(U_x) 
  # we project the residual to be orthogonal to X
  # R <- (diag(1, m) - M_x) %*% (Y - X %*% B_star_hat) 
  R <- Y - X %*% B_star_hat
  # M_x <- U_x %*% t(U_x) 
  # 
  # (diag(1, m) - M_x) %*% Y - R
  
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
  Psi_hat <- apply(R, 2, function(r) {lsei::lsei(a = U_q, b = r)})
  M_jn <- matrix(1/n, n, n)
  D_jn <- diag(1, n) - M_jn
  
  if (q == 1) {
    if (Psi_hat %*% D_jn %*% as.matrix(Psi_hat) < 1e-15) {
      message("q_hat = 1 and Psi_hat %*% D_jn %*% as.matrix(Psi_hat) is singular, Gamma_hat is assigned to equal U_q.")
      Gamma_hat <- U_q
    } else {
      Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% as.matrix(Psi_hat) %*% (Psi_hat %*% D_jn %*% as.matrix(Psi_hat))^(-1)  
    }
  } else {
    if (kappa(Psi_hat %*% D_jn %*% t(Psi_hat)) < 1e-15) {
      message("q_hat > 1 and Psi_hat %*% D_jn %*% t(Psi_hat) is singular, Gamma_hat is assigned to equal U_q.")
      Gamma_hat <- U_q
    } else {
      Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) %*% solve(Psi_hat %*% D_jn %*% t(Psi_hat))
    }
  }
  # Gamma_hat 
  
  ## step 4: fitting the model again with the surrogate variable
  X_new <- cbind(Gamma_hat, X)
  B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = q, sum = 1)$x})
  # B_hat <- apply(Y, 2, function(y) {lsei::lsei(a = X_new, b = y)})
  # D_hat <- B_hat[seq(q), ]
  P_hat <- B_hat[-seq(q), ]
  
  # Gamma_hat %*% D_hat - Z %*% D
  
  ## return the P_hat
  return(list(P_hat = P_hat,
              q_hat = q))
}

## another version of the dSVA with intercept 
dsva_ext3_know_q <- function(Y, Theta, q = 1) {
  
  # Y <- true_data$Y
  # Theta <- true_data$X
  # q <- 1
  
  n <- ncol(Y)
  m <- nrow(Y)
  
  ## build the regression coefficient matrix
  X <- cbind(diag(1, m, m), Theta)
  
  ## step 1: obtain the canonical model residual 
  B_star_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, k = m, sum = 1)$x})
  U_x <- svd(X)$u
  M_x <- U_x %*% t(U_x) 
  # we project the residual to be orthogonal to X
  R <- (diag(1, m) - M_x) %*% (Y - X %*% B_star_hat) 
  
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
  Psi_hat <- apply(R, 2, function(r) {lsei::lsei(a = U_q, b = r)})
  M_jn <- matrix(1/n, n, n)
  D_jn <- diag(1, n) - M_jn
  
  if (q == 1) {
    if (Psi_hat %*% D_jn %*% as.matrix(Psi_hat) > 1e-15) {
      message("q_hat = 1 and Psi_hat %*% D_jn %*% as.matrix(Psi_hat) is singular, Gamma_hat is assigned to equal U_q.")
      Gamma_hat <- U_q
    } else {
      Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% as.matrix(Psi_hat) %*% solve(Psi_hat %*% D_jn %*% as.matrix(Psi_hat))  
    }
  } else {
    if (kappa(Psi_hat %*% D_jn %*% t(Psi_hat)) > 1e-15) {
      message("q_hat > 1 and Psi_hat %*% D_jn %*% t(Psi_hat) is singular, Gamma_hat is assigned to equal U_q.")
      Gamma_hat <- U_q
    } else {
      Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) %*% solve(Psi_hat %*% D_jn %*% t(Psi_hat))
    }
  }
  # Gamma_hat 
  
  ## step 4: fitting the model again with the surrogate variable (the first m + q variables are not NN-restricted)
  X_new <- cbind(Gamma_hat, X)
  B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = m + q, sum = 1)$x})
  P_hat <- B_hat[-seq(m + q), ]
  
  ## return the P_hat
  return(list(P_hat = P_hat,
              q_hat = q))
}