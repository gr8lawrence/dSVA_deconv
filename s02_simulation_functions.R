source("s_sim_functions.R")

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
  Z[, 2] <- c(rnorm(floor(m/2), mean = -mu, sd = 10), rnorm(floor(m/2), mean = mu, sd = 10))

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
  Y[Y < 0] <- runif(n = sum(Y < 0), min = 0, max = 1e-16)
  
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

## row sums of X are always one; no scaling 
dSVA_deconv_gene2 <- function(m, n, K, p, p_sig, p_class, rho, lambda, mu) {
  ## generate signature matrices
  ## (1) no signature (2) 50% signatures (3) all signatures
  
  X <- matrix(nrow = m, ncol = K)
  for (i in 1:m) {
    ii <- extraDistr::rbern(1, p_sig)
    X[i, ] <- sample(ii * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - ii) * extraDistr::rdirichlet(1, rep(1, K)))
  }
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
  sigma_sq <- statmod::rinvgauss(m, 1, 1)
  for (i in 1:m) E[i, ] <- rnorm(n = n, mean = 0, sd = sqrt(sigma_sq[i]))
  
  Y <- X %*% P_star + Z %*% D_star + E
  Y[Y < 0] <- runif(n = sum(Y < 0), min = 0, max = 1e-16)

  ## gather the results into a list
  list(X = X,
       Z = Z,
       P_star = P_star,
       D_star = D_star,
       E = E,
       Y = Y)
}

## column sums of X are always one; unable to pick up the scaling by constrained NNLS
dSVA_deconv_gene3 <- function(m, n, K, p, p_sig, p_class, rho, lambda, mu) {
  ## generate signature matrices
  ## (1) no signature (2) 50% signatures (3) all signatures
  
  D_exp <- diag(rchisq(n = m, df = 200)) # mean expression value of each gene
  X <- matrix(nrow = m, ncol = K)
  for (i in 1:m) {
    ii <- extraDistr::rbern(1, p_sig)
    X[i, ] <- sample(ii * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - ii) * extraDistr::rdirichlet(1, rep(1, K)))
  }
  X <- D_exp %*% X
  X <- apply(X, 2, function(x) x/sum(x))
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
  sigma_sq <- statmod::rinvgauss(m, 1, 1)
  for (i in 1:m) E[i, ] <- rnorm(n = n, mean = 0, sd = sqrt(sigma_sq[i]))
  
  Y <- X %*% P_star + Z %*% D_star + E
  Y[Y < 0] <- runif(n = sum(Y < 0), min = 0, max = 1e-16)
  
  ## gather the results into a list
  list(X = X,
       Z = Z,
       P_star = P_star,
       D_star = D_star,
       E = E,
       Y = Y)
}

## if we transpose the model
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
  Z[, 2] <- c(rnorm(floor(n/2), mean = -mu, sd = 10), rnorm(floor(n/2), mean = mu, sd = 10))
  
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
  Y[Y < 0] <- runif(n = sum(Y < 0), min = 0, max = 1e-16)
  
  # var(c(t(P_star) %*% t(X)))
  # var(c(Z %*% D_star))
  
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
my_ccc <- function(P1, P2) {
  ccc <- vector(mode = "double", length = ncol(P1))
  for (i in 1:ncol(P1)) ccc[i] <- DescTools::CCC(x = P1[, i], y = P2[, i])$rho.c[1, 1]
  ccc
}
my_mse <- function(P1, P2) mean((P1 - P2)^2)

# P1 = P2 = matrix(1:9, 3, 3)
# DescTools::CCC(x = P1[, 2], y = P2[, 1])$rho.c[1, 1]

## finding a function where the simulated latent variables can carry significant variation
## for simplicity, assuming each gene has a total expression of 1
## m >> n in this case when latent variables are on samples
## p_sig: the proportions of marker genes
## lambda: how distinguished the abundances of cell types are

dSVA_model_sim <- function(m, n, K, p = 2, p_sig = 0.5, lambda = 1, mu = 5, sigma = 5) {
  
  # test parameters
  n = 20
  m = 1000
  K = 5
  p = 3
  p_sig = 0.5
  lambda = 1

  ## generate signature matrices
  ## (1) no signature (2) 50% signatures (3) all signatures
  X <- matrix(nrow = m, ncol = K)
  for (i in 1:m) {
    ii <- extraDistr::rbern(1, p_sig)
    X[i, ] <- sample(ii * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - ii) * extraDistr::rdirichlet(1, rep(1, K)))
  }  
  
  ## generate proportions
  P_star <- matrix(nrow = K, ncol = n)
  for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = lambda^seq(K))
  
  ## generate the categorical and continuous latent variables in the D matrix
  Z <- matrix(nrow = m, ncol = p)
  D <- matrix(nrow = p, ncol = n)
  D[1, ] <- c(rep(0, floor(n/2)), rep(1, ceiling(n/2)))
  D[2, ] <- rep(c(rep(0, floor(n/4)), rep(1, ceiling(n/4))), 2)
  D[3, ] <- 1 - D[2, ]
  Z[, 1] <- rnorm(m, 0, sigma/K^2)
  Z[, 2] <- rnorm(m, -mu/K^2, sigma/K^2)
  Z[, 3] <- rnorm(m, mu/K^2, sigma/K^2)
  
  # Z[, 1] <- c(rep(0, floor(m/2)), rep(1, ceiling(m/2)))
  # Z[, 2] <- rep(c(rep(0, floor(m/4)), rep(1, ceiling(m/4))), 2)
  # Z[, 3] <- 1 - Z[, 2]
  # D[1, ] <- rnorm(n, 0, 1/K^2)
  # D[2, ] <- rnorm(n, -1/K^2, 1/K^2)
  # D[3, ] <- rnorm(n, 1/K^2, 1/K^2)
  
  ## generate 
  Y_reg <- X %*% P_star
  Y_lat <- Z %*% D
  
  # Y_reg + Y_lat
  # var(c(Y_reg))/var(c(Y_lat))
  
  ## generate the measurement error
  # E <- matrix(nrow = m, ncol = n)
  # sigma_sq <- statmod::rinvgauss(m, 1, 1)
  # for (i in 1:m) E[i, ] <- rnorm(n = n, mean = 0, sd = sqrt(sigma_sq[i]))
  # for (i in 1:m) E[i, ] <- rnorm(n = n, mean = 0, sd = sqrt(0.01))
  Y <- Y_reg + Y_lat 
  Y[Y < 0] <- runif(n = sum(Y < 0), min = 0, max = 1e-16)

  ## gather the results into a list
  list(X = X,
       Z = Z,
       P_star = P_star,
       D = D,
       # E = E,
       Y = Y)
}

## a new version of the simulation
## use relative contributions to expression of genes of latent variables - may need an intercept
dSVA_model_sim_2 <- function(m, n, K, p = 1, p_sig = 0.5, lambda = 3, gamma = 2) {
  
  # test parameters
  # n = 20
  # m = 1000
  # K = 5
  # p = 1
  # p_sig = 0.5
  # mu = 1
  # sigma = 3
  # lambda = 3

  ## generate signature matrices
  ## (1) no signature (2) 50% signatures (3) all signatures
  X <- matrix(nrow = m, ncol = K)
  W <- extraDistr::rbern(m, p_sig)
  D_exp <- diag(rchisq(n = m, df = 200)) # mean expression value of each gene
  for (i in 1:m) X[i, ] <- sample(W[i] * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - W[i]) * extraDistr::rdirichlet(1, rep(1, K)))
  X <- D_exp %*% X
  
  ## generate proportions
  P_star <- matrix(nrow = K, ncol = n)
  for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = rep(1, K))
  # for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = lambda^seq(K))
  
  ## generate the categorical and continuous latent variables in the D matrix
  Z <- matrix(nrow = m, ncol = p)
  D <- matrix(nrow = p, ncol = n)
  D[1, ] <- c(rep(0, floor(n/2)), rep(1, ceiling(n/2)))
  # Z[, 1] <- rinvgauss(n = m, mean = mu, dispersion = sigma)
  Z[, 1] <- W * rchisq(n = m, df = gamma * 200) + (1 - W) * rchisq(n = m, df = 200)
  
  ## generate 
  Y_reg <- X %*% P_star
  Y_lat <- Z %*% D
  
  ## generate the measurement error
  # E <- matrix(nrow = m, ncol = n)
  # sigma_sq <- statmod::rinvgauss(m, 1, 1)
  # for (i in 1:m) E[i, ] <- rnorm(n = n, mean = 0, sd = sqrt(sigma_sq[i]))
  # for (i in 1:m) E[i, ] <- rnorm(n = n, mean = 0, sd = sqrt(0.01))
  Y <- Y_reg + Y_lat 
  # Y[Y < 0] <- runif(n = sum(Y < 0), min = 0, max = 1e-16)
  # 
  ## gather the results into a list
  ls <- list(X = X,
             Z = Z,
             P_star = P_star,
             D = D,
             Y = Y)
  return(ls)
}

## this version strengthens the expression of the first cell type
dSVA_model_sim_ct_strengthen<- function(m, n, K, p = 1, p_sig = 0.5, lambda = 3, gamma = 2) {
  
  # test parameters
  # n = 20
  # m = 1000
  # K = 5
  # p = 1
  # p_sig = 0.5
  # mu = 1
  # sigma = 3
  # lambda = 3
  
  ## generate signature matrices
  ## (1) no signature (2) 50% signatures (3) all signatures
  X <- matrix(nrow = m, ncol = K)
  W <- extraDistr::rbern(m, p_sig)
  D_exp <- diag(rchisq(n = m, df = 200)) # mean expression value of each gene
  for (i in 1:m) X[i, ] <- sample(W[i] * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - W[i]) * extraDistr::rdirichlet(1, rep(1, K)))
  X <- D_exp %*% X
  
  ## generate proportions
  P_star <- matrix(nrow = K, ncol = n)
  for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = rep(1, K))
  # for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = lambda^seq(K))
  
  ## generate the categorical and continuous latent variables in the D matrix
  Z <- matrix(nrow = m, ncol = p)
  D <- matrix(nrow = p, ncol = n)
  D[1, ] <- c(rep(0, floor(n/2)), rep(1, ceiling(n/2)))
  Z[, 1] <- gamma * X[, 1]
  # Z[, 1] <- W * rchisq(n = m, df = gamma * 200) + (1 - W) * rchisq(n = m, df = 200)
  
  ## generate 
  Y_reg <- X %*% P_star
  Y_lat <- Z %*% D
  
  ## generate the measurement error
  # E <- matrix(nrow = m, ncol = n)
  # sigma_sq <- statmod::rinvgauss(m, 1, 1)
  # for (i in 1:m) E[i, ] <- rnorm(n = n, mean = 0, sd = sqrt(sigma_sq[i]))
  # for (i in 1:m) E[i, ] <- rnorm(n = n, mean = 0, sd = sqrt(0.01))
  Y <- Y_reg + Y_lat 
  # Y[Y < 0] <- runif(n = sum(Y < 0), min = 0, max = 1e-16)
  # 
  ## gather the results into a list
  ls <- list(X = X,
             Z = Z,
             P_star = P_star,
             D = D,
             Y = Y)
  return(ls)
}

## build a new model with small increments on genes
dSVA_model_sim <- function(m, n, K, p = 1, p_sig = 0.5, lambda = 3, gamma = 2, chi_df = 200) {
  
  ## generate signature matrices
  X <- matrix(nrow = m, ncol = K)
  W <- extraDistr::rbern(m, p_sig)
  D_exp <- diag(rchisq(n = m, df = chi_df)) # mean expression value of each gene
  for (i in 1:m) X[i, ] <- sample(W[i] * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - W[i]) * extraDistr::rdirichlet(1, rep(1, K)))
  X <- D_exp %*% X
  
  ## generate proportions
  P_star <- matrix(nrow = K, ncol = n)
  for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = rep(1, K))
  
  ## generate the categorical and continuous latent variables in the D matrix
  Z <- matrix(nrow = m, ncol = p)
  D <- matrix(nrow = p, ncol = n)
  D[1, ] <- c(rep(0, floor(n/2)), rep(1, ceiling(n/2)))
  Z[, 1] <- rnorm(n, mu, sqrt(chi_df))
  
  ## generating the regression part and the latent variable part
  Y_reg <- X %*% P_star
  Y_lat <- Z %*% D
  
  ## generate the measurement error)
  Y <- Y_reg + Y_lat 
  
  ## gather the results into a list
  ls <- list(X = X,
             Z = Z,
             P_star = P_star,
             D = D,
             Y = Y)
  return(ls)
}

## build a new model with small deviations
## Question to ask: would small additions on genes lead to larger deviations of results?
dSVA_model_small_dev <- function(m, n, K, p = 1, p_sig = 0.5, lambda = 3, gamma = 2, chi_df = 200) {
  
  ## generate signature matrices
  X <- get_signatures(m, n, K, p_sig, lambda, chi_df)$X
  
  ## generate proportions
  P_star <- matrix(nrow = K, ncol = n)
  for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = rep(1, K))
  
  ## generate the categorical and continuous latent variables in the D matrix
  Z <- matrix(nrow = m, ncol = p)
  D <- matrix(nrow = p, ncol = n)
  D[1, ] <- c(rep(0, floor(n/2)), rep(1, ceiling(n/2)))
  Z[, 1] <- rnorm(n, 0, sqrt(chi_df/3))
  
  ## generating the regression part and the latent variable part
  Y_reg <- X %*% P_star
  Y_lat <- Z %*% D
  
  ## generate the measurement error)
  Y <- Y_reg + Y_lat 
  
  ## force neg values to nonneg (zero)
  Y[Y < 0] <- 0
  
  ## gather the results into a list
  ls <- list(X = X,
             Z = Z,
             P_star = P_star,
             D = D,
             Y = Y)
  return(ls)
}

## build a new model with an intercept
dSVA_model_sim_intercept <- function(m, n, K, p = 1, p_sig = 0.5, lambda = 3, gamma = 2, chi_df = 200) {
  
  ## generate signature matrices
  sig_ls <- get_signatures(m, n, K, p_sig, lambda, chi_df)
  X <- sig_ls$X
  W <- sig_ls$W
  
  ## generate proportions
  P_star <- matrix(nrow = K, ncol = n)
  for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = rep(1, K))
  
  ## generate the categorical and continuous latent variables in the D matrix
  Z <- matrix(nrow = m, ncol = p)
  D <- matrix(nrow = p, ncol = n)
  D[1, ] <- c(rep(0, floor(n/2)), rep(1, ceiling(n/2)))
  Z[, 1] <- W * rchisq(n = m, df = gamma * chi_df) + (1 - W) * rchisq(n = m, df = chi_df)
  
  ## generating the regression part and the latent variable part
  Y_reg <- X %*% P_star
  Y_lat <- Z %*% D
  
  ## generate the measurement error)
  Y <- Y_reg + Y_lat 
  
  ## gather the results into a list
  ls <- list(X = X,
             Z = Z,
             P_star = P_star,
             D = D,
             Y = Y)
  return(ls)
}
