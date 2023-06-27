## simulation based on the dSVA paper. 
set.seed(100)
library(extraDistr)
library(MASS)

## algorithm parameters
m <- 1000
n <- 20
K <- 5
p <- 2
p_sig <- 0.5
p_class <- 0.5
rho <- 0.5
lambda <- 5

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

Y <- X %*% P_star + Z %*% D_star + E  