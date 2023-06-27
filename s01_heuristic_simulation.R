## simulation based on the dSVA paper. 
set.seed(100)
library(extraDistr)

## algorithm parameters
m <- 1000
n <- 20
K <- 5
p <- 2
p_sig <- 0.5
p_class <- 0.5
lambda <- 5

## generate signature matrices
## (1) no signature (2) 50% signature (3) all signatures
X <- matrix(nrow = m, ncol = K)
D_exp <- diag(rchisq(n = m, df = 200)) # mean expression value of each gene
for (i in 1:m) {
  ii <- rbern(1, p_sig)
  X[i, ] <- sample(ii * rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - ii) * rdirichlet(1, rep(1, K)))
}
X <- D_exp %*% X 

## generate the categorical and continuous variables


