dSVA_deconv_gene <- function(m, n, K, p, p_sig, p_class, rho, lambda, mu) {
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
  # Z[, 1] <- extraDistr::rbern(m, p_class)
  Z[, 1] <- c(rep(0, floor(m/2)), rep(1, ceiling(m/2)))
  Z[, 2] <- c(rnorm(floor(m/2), mean = -mu, sd = 1), rnorm(floor(m/2), mean = mu, sd = 1))

  ## generate proportions
  coefs_star <- extraDistr::rdirichlet(n = n, alpha = c(seq(1, K), seq(rho, p)))
  P_star <- matrix(nrow = K, ncol = n)
  for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = coefs_star[i, seq(K)])
  D_star <- matrix(nrow = p, ncol = n)
  for (i in 1:n) D_star[, i] <- MASS::mvrnorm(n = 1, mu = c(0, 0), Sigma = diag(coefs_star[i, seq(K + 1, K + p)]))

  ## generate the measurement error
  E <- matrix(nrow = m, ncol = n)
  sigma_sq <- statmod::rinvgauss(m, 10, 9)
  for (i in 1:m) E[i, ] <- rnorm(n = n, mean = 0, sd = sqrt(sigma_sq[i]))
  
  Y <- X %*% P_star + Z %*% D_star + E
  
  # X <- matrix(rchisq(m * K, df = 10), m, K)
  # B <- matrix(rchisq(K * n, df = 5), K, n)
  # P_star <- apply(B, 2, function(b) b/sum(b))
  # Z <- matrix(abs(rnorm(m * p, mean = 2)), m, p)
  # Z[, 1] <- extraDistr::rbern(m, p_class)
  # D_star <- matrix(abs(rnorm(p * n)), p, n)
  # E <- matrix(rnorm(m * n, 0, 0.01), m, n)
  # Y <- X %*% P_star + Z %*% D_star + E
  # Y[Y < 0] <- 0
  
  ## gather the results into a list
  list(X = X,
       Z = Z,
       P_star = P_star,
       D_star = D_star,
       E = E,
       Y = Y)
}


dSVA_deconv_sim_people <- function(m, n, K, p, p_sig, p_class, rho, lambda, mu) {
  ## generate signature matrices
  ## (1) no signature (2) 50% signatures (3) all signatures
  
  ## True signature matrix
  X <- matrix(nrow = m, ncol = K)
  D_exp <- diag(rchisq(n = m, df = 200)) # mean expression value of each gene
  for (i in 1:m) {
    ii <- extraDistr::rbern(1, p_sig)
    X[i, ] <- sample(ii * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - ii) * extraDistr::rdirichlet(1, rep(1, K)))
  }
  X <- D_exp %*% X
  # U_X <- svd(X)$u ## get the SVD of X
  
  ## generate the categorical and continuous latent variables
  Z <- matrix(nrow = n, ncol = p)
  # Z[, 1] <- extraDistr::rbern(m, p_class)
  Z[, 1] <- c(rep(0, floor(n/2)), rep(1, ceiling(n/2)))
  Z[, 2] <- c(rnorm(floor(n/2), mean = -mu, sd = 1), rnorm(floor(n/2), mean = mu, sd = 1))
  
  ## generate proportions
  # coefs_star <- extraDistr::rdirichlet(n = n, alpha = c(seq(1, K), seq(rho, p)))
  P_alpha <- extraDistr::rdirichlet(n = 1, alpha = seq(K))
  P_star <- t(extraDistr::rdirichlet(n = n, alpha = P_alpha))
  D_star <- matrix(nrow = p, ncol = m)
  D_alpha <- extraDistr::rdirichlet(n = 1, alpha = rep(1, p))
  D_star <- t(MASS::mvrnorm(n = m, mu = c(0, 0), Sigma = diag(c(D_alpha))))
  
  ## generate the measurement error
  E <- matrix(nrow = n, ncol = m)
  sigma_sq <- statmod::rinvgauss(n, 10, 9)
  for (i in 1:n) E[i, ] <- rnorm(n = m, mean = 0, sd = sqrt(sigma_sq[i]))
  
  Y <- t(P_star) %*% t(X) + Z %*% D_star + E
  
  ## gather the results into a list
  list(X_star = t(X),
       Z = Z,
       P = t(P_star),
       D_star = D_star,
       E = E,
       Y = Y)
}

## function for column-wise correlation and MSE between P1 and P2
my_cor <- function(P1, P2) diag(cor(P1, P2))
my_mse <- function(P1, P2) mean((P1 - P2)^2)
