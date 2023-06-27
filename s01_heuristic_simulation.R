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
# S <- 5

true_data <- dSVA_deconv_sim(m, n, K, p, p_sig, p_class, rho, lambda)
P_dSVA <- dsva_ext(Y = true_data$Y, X = true_data$X, q = p)
P_nnls <- apply(true_data$Y, 2, function(y) {lsei::nnls(a = true_data$X, b = y)$x})
P_pnnls <- apply(true_data$Y, 2, function(y) {lsei::pnnls(a = true_data$X, b = y, sum = 1)$x})

