## function for creating disturbances
disturb <- function(alpha) {
  K <- length(alpha)
  alpha2 <- vector(mode = "double", length = K)
  alpha2[seq(K - 1)] <- alpha[seq(K - 1)] + runif(K - 1, 0, mean(alpha)/(K - 1))
  alpha2[K] <- alpha[K] - runif(1, 0, mean(alpha))
  alpha2
}

softmax <- function(x) exp(x)/sum(exp(x))
row_center <- function(X) t(apply(X, 1, function(x) x - mean(x)))


## function to check assumption 5
check_a5 <- function(X, Y) {
  # the K x q covariance matrix
  mat_out <- cov_mat <- matrix(NA, nrow = nrow(X), ncol = nrow(Y))
  for (k in 1:nrow(X)) {
    for (l in 1:nrow(Y)) {
      cov_mat[k, l] <- cov(X[k, ], Y[l, ])
    }
  }
  
  # the q x 1 variance vector
  Y_var <- apply(Y, 1, var)
  
  # calculate the ratios
  for (l in 1:nrow(Y)) mat_out[, l] <- cov_mat[, l]/Y_var[l]
  
  # return the matrix of ratios
  mat_out
}


# For the case-control structure ------------------------------------------


## Using the hierarchical exponential generative model for P
# n = 10
n_seq <- c(10, 100, 1000, 10000)
ratio_mat <- matrix(NA, nrow = 100, ncol = length(n_seq))
colnames(ratio_mat) <- as.character(n_seq)
for (n in n_seq) {
  for (b in 1:100) {
    i <- which(n == n_seq)
    Alpha <- Alpha2 <- t(extraDistr::rdirichlet(n, 1:5))
    Alpha2[, seq(n/2 + 1, n)] <- apply(Alpha2[, seq(n/2 + 1, n)], 2, disturb)
    P_star <- apply(Alpha, 2, softmax)
    P1 <- apply(Alpha2, 2, softmax)
    D <- P1 - P_star
    P_star_tilde <- row_center(P_star)
    D_tilde <- row_center(D)
    
    ratio_mat[b, i] <- max(abs(check_a5(P_star_tilde, D_tilde))) 
  }
}

boxplot(ratio_mat) # it did show improvement!


# For the binary confounding ----------------------------------------------



# For the continuous confounding ------------------------------------------

n_seq <- c(10, 100, 1000, 10000)
ratio_mat <- matrix(NA, nrow = 100, ncol = length(n_seq))
for (n in n_seq) {
  for (b in 1:100) {
    i <- which(n == n_seq)
    P_star <- t(extraDistr::rdirichlet(n = n, alpha = rep(1, 5)))
    D <- matrix(rchisq(n, df = 200/5), nrow = 1)
    ratio_mat[b, i] <- max(abs(check_a5(P_star, D))) 
  }
}
boxplot(ratio_mat, main = "Continuous Latent Structure") # also shows improvement!


# cor_mat <- cor(t(P_star_tilde), t(D_tilde))
# var_mat <- cor(t(D_tilde), t(D_tilde))
# delt <- cor_mat/var_mat



# P1[, seq(n/2 + 1, n)] <- apply(P1[, seq(n/2 + 1, n)], function(x) disturb(x))
