## This file contains all the official simulation functions
get_signatures <- function(m, n, K, p_sig = 0.5, lambda = 3, chi_df = 200) {
  X <- matrix(NA_real_, nrow = m, ncol = K)
  W <- extraDistr::rbern(m, p_sig)
  D_exp <- diag(rchisq(n = m, df = chi_df)) # mean expression value of each gene
  for (i in 1:m) X[i, ] <- sample(W[i] * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - W[i]) * extraDistr::rdirichlet(1, rep(1, K)))
  X <- D_exp %*% X
  list(X = X, W = W)
}

