dSVA_deconv_sim <- function(m, n, K, p, p_sig, p_class, rho, lambda) {
  ## generate signature matrices
  ## (1) no signature (2) 50% signatures (3) all signatures
  X <- matrix(nrow = m, ncol = K)
  D_exp <- diag(rchisq(n = m, df = 200)) # mean expression value of each gene
  for (i in 1:m) {
    ii <- extraDistr::rbern(1, p_sig)
    X[i, ] <- sample(ii * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - ii) * extraDistr::rdirichlet(1, rep(1, K)))
  }
  X <- D_exp %*% X 
  # U_X <- svd(X)$u ## get the SVD of X
  
  ## generate the categorical and continuous latent variables
  Z <- matrix(nrow = m, ncol = p)
  Z[, 1] <- extraDistr::rbern(m, p_class)
  Z[, 2] <- rnorm(m, mean = 0, sd = 1)
  
  ## generate proportions
  coefs_star <- extraDistr::rdirichlet(n = n, alpha = c(seq(1, K), seq(rho, p)))
  P_star <- matrix(nrow = K, ncol = n)
  for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = coefs_star[i, seq(K)])
  D_star <- matrix(nrow = p, ncol = n)
  for (i in 1:n) D_star[, i] <- MASS::mvrnorm(n = 1, mu = c(0, 0), Sigma = diag(coefs_star[i, seq(K + 1, K + p)]))
  
  ## generate the measurement error
  E <- matrix(nrow = m, ncol = n)
  Y_star <- X %*% P_star 
  for (i in 1:m) E[i, ] <- rnorm(n = n, mean = 0, sd = 0.1)
  Y <- Y_star + Z %*% D_star + E  
  Y[Y < 0] <- 0
  
  ## gather the results into a list
  list(X = X,
       Z = Z,
       P_star = P_star,
       D_star = D_star,
       E = E,
       Y = Y)
}

