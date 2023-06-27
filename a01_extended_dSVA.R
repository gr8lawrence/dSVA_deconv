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

dsva_ext <- function(Y, X, q) {
  
  ## step 1: obtain the canonical model residual 
  B_star_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, sum = 1)$x})
  M_x <- X %*% solve(t(X) %*% X) %*% t(X)
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

# mean(abs(B - P_hat))
# mean(abs(B - B_star_hat))

