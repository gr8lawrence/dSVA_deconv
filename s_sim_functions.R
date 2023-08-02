## This file contains all the official simulation functions
my_cor <- function(P1, P2) mean(diag(cor(P1, P2)))
my_ccc <- function(P1, P2) {
  ccc <- vector(mode = "double", length = ncol(P1))
  for (i in 1:ncol(P1)) ccc[i] <- DescTools::CCC(x = P1[, i], y = P2[, i])$rho.c[1, 1]
  mean(ccc)
}
my_mse <- function(P1, P2) mean((P1 - P2)^2)

get_signatures <- function(m, n, K, p_sig = 0.5, lambda = 3, chi_df = 200) {
  X <- matrix(NA_real_, nrow = m, ncol = K)
  W <- extraDistr::rbern(m, p_sig)
  D_exp <- diag(rchisq(n = m, df = chi_df)) # mean expression value of each gene
  for (i in 1:m) X[i, ] <- sample(W[i] * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - W[i]) * extraDistr::rdirichlet(1, rep(1, K)))
  X <- D_exp %*% X
  list(X = X, W = W)
}

## TODO: simulate when q = 2, 3, 4 (1 discrete + 1 continuous, 1 discrete + 2 continuous, 2 discrete + 2 continuous)
## TODO: simulate continuous functions only when q = 1.
## allow errors (err = TRUE)
## build a new model with an intercept
dSVA_model_sim_intercept <- function(m, n, K, q = 1:4, p_sig = 0.5, lambda = 3, gamma = 2, chi_df = 200, err = FALSE, first_effect = c("bin", "con", "small"), p_w2 = 0.5) {
  
  ## generate signature matrices
  sig_ls <- get_signatures(m, n, K, p_sig, lambda, chi_df)
  X <- sig_ls$X
  W <- sig_ls$W
  
  ## generate proportions
  P_star <- matrix(nrow = K, ncol = n)
  for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = rep(1, K))
  
  ## generating the regression part and the latent variable part
  Y_reg <- X %*% P_star
  
  ## generate the categorical and continuous latent variables in the D matrix
  # Z <- matrix(nrow = m, ncol = q)
  D <- matrix(nrow = q, ncol = n)
  D[1, ] <- c(rep(0, floor(n/2)), rep(1, ceiling(n/2)))
  if (q >= 2) {
    ## add another continuous feature (purely positive)
    D[2, ] <- rchisq(n = n, df = 30)
    # D[2, ] <- abs(rnorm(n = n, mean = 5, sd = 1))
  }  
  if (q >= 3) {
    ## add another continuous feature (could be negative)
    # D[3, ] <- 
  }
  if (q == 4) {
    ## add another discrete feature
    # D[4, ] <- 
  }
  ## choose half of the 0, 1 group to be up regulated
  
  ## encode the latent effects of the binary groups directly instead of using Z * D
  if (first_effect == "small") {
    Y_lat <- get_Y_lat_small_effects(m, n, chi_df, gamma, W)
  } else if (first_effect == "bin") {
    Y_lat <- get_Y_lat_binary(m, n, chi_df, gamma, W)
  } else if (first_effect == "con") {
    Y_lat <- get_Y_lat_continuous(m, n, chi_df, gamma, W, p_w2)
  }
  
  if (q >= 2) {
    Y_lat2 <- matrix(0, nrow = m, ncol = n)
    W2 <- extraDistr::rbern(m, 0.5)
    for (i in 1:n) {
      # Y_lat2[, i] <- W2 * (W * rchisq(m, df = gamma * chi_df/4) + (1 - W) * rchisq(m, df = chi_df/4)) * D[2, i]
      # Y_lat2[, i] <- W2 * (W * rnorm(m, gamma * chi_df, sqrt(chi_df/3)) + (1 - W) * rnorm(m, chi_df, sqrt(chi_df/3))) * D[2, i]
      Y_lat2[, i] <- W2 * (W * rnorm(m, gamma * chi_df/4, sqrt(gamma * chi_df/3)) + (1 - W) * rnorm(m, chi_df/4, sqrt(chi_df/3))) * D[2, i]
    }
    Y_lat <- Y_lat + Y_lat2
  }
  # if (q >= 3) {
  #   Y_lat3 <- matrix(0, nrow = m, ncol = n)
  #   W3 <- extraDistr::rbern(m, 0.4)
  # }
  
  ## generate the measurement error
  Y <- Y_reg + Y_lat 
  
  if (err) {
    E <- matrix(0, nrow = m, ncol = n)
    for (i in 1:m) {
      E[i, ] <- rnorm(n, mean(Y_reg[i, ]), sqrt(mean(Y_reg[i, ])/5))
    }
    Y <- Y + E
  } else {
    E <- NA
  }
  
  ## get the proportion of Y being negative
  p_neg <- sum(Y < 0)/(m * n)
  Y[Y < 0] <- 0
  
  
  ## gather the results into a list
  ls <- list(X = X,
             # Z = Z, # Z is not important!
             P_star = P_star,
             D = D,
             Y = Y,
             E = E,
             Y_lat = Y_lat,
             p_neg = p_neg)
  return(ls)
}

get_Y_lat_binary <- function(m, n, chi_df, gamma, W) {
  Y_lat <- matrix(0, m, n)
  
  ## for the group where D = 0, select m/4 genes that increase in expression
  m1 <- ceiling(m/4)
  n1 <- floor(n/2)
  W1 <- W[seq(m1)]
  Y_lat[seq(m1), seq(n1)] <- W1 * rchisq(n = m1 * n1, df = gamma * chi_df) + (1 - W1) * rchisq(n = m1 * n1, df = chi_df)
  
  ## for the group where D = 1, do the same
  m2 <- floor(m/4)
  n2 <- ceiling(n/2)
  W2 <- W[seq(m - m2 + 1, m)]
  Y_lat[seq(m - m2 + 1, m), seq(n - n2 + 1, n)] <- W2 * rchisq(n = m2 * n2, df = gamma * chi_df) + (1 - W2) * rchisq(n = m2 * n2, df = chi_df)
  
  ## return Y_lat
  Y_lat
}

get_Y_lat_small_effects <- function(m, n, chi_df, gamma, W) {
  Y_lat <- matrix(0, m, n)
  
  ## for the group where D = 0, select m/4 genes that increase in expression
  m1 <- ceiling(m/4)
  n1 <- floor(n/2)
  W1 <- W[seq(m1)]
  Y_lat[seq(m1), seq(n1)] <- W1 * rnorm(n = m1 * n1, mean = gamma * chi_df, sd = sqrt(chi_df/3)) + (1 - W1) * rnorm(n = m1 * n1, mean = chi_df, sd = sqrt(chi_df/3))
  
  ## for the group where D = 1, do the same
  m2 <- floor(m/4)
  n2 <- ceiling(n/2)
  W2 <- W[seq(m - m2 + 1, m)]
  Y_lat[seq(m - m2 + 1, m), seq(n - n2 + 1, n)] <- W2 * rnorm(n = m2 * n2, mean = gamma * chi_df, sd = sqrt(chi_df/3)) + (1 - W2) * rnorm(n = m2 * n2, mean = chi_df, sd = sqrt(chi_df/3))
  
  ## return Y_lat
  Y_lat
}

get_Y_lat_continuous <- function(m, n, chi_df, gamma, W, p_w2 = 0.5) {
  Y_lat2 <- matrix(0, nrow = m, ncol = n)
  W2 <- extraDistr::rbern(m, p_w2)
  for (i in 1:n) {
    Y_lat2[, i] <- W2 * (W * rnorm(m, gamma * chi_df/4, sqrt(gamma * chi_df/3)) + (1 - W) * rnorm(m, chi_df/4, sqrt(chi_df/3))) * D[2, i]
  }
  Y_lat2
}
