## this page is used to understand the lsei package
# source("s_sources.R")
# set.seed(100)
# n <- 20
# m <- 1000
# K <- 5
# q <- 1
# err <- TRUE
# p <- 0.75
# lambda <- 5
# gamma <- 4
# n_sv <- 6 # the number of sv we want to plot
# first_effect <- "bin"
# true_data <- dSVA_model_sim_intercept(m, n, K, q, p, lambda, gamma, err = err, first_effect = first_effect)
# 
# ## use lsei on the
# Theta <- true_data$X
# X <- model.matrix(~1 + Theta)
# P_hat <- apply(true_data$Y, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(Theta), sum = 1)$x})
# 
# ## example given in lsei
# a = matrix(rnorm(40), nrow = 10)
# b = drop(a %*% c(0,1,-1,1)) + rnorm(10)
# p1 <- lsei::pnnls(a, b, k=2, sum=1)$x
# r <- b - a[, 3:4] %*% p1[3:4] 
# # mat <- cbind(r, a[, 1:2])
# # colnames(mat) <- c("r", "a1", "a2")
# lm1 <- lm(r ~ a[, 1:2] - 1)
# all.equal(coef(lm1), p1[1:2])


## this section checks what happens for continuous confounding
Y <- true_data$Y
Y_true <- true_data$X %*% true_data$P_star + true_data$E 
all.equal(Y - Y_true, true_data$Y_lat) # passed
X <- model.matrix(~ 1 + true_data$X)
P_known <- apply(Y_true, 2, function(y) {lsei::pnnls(a = X, b = y , sum = 1)$x})
P_pnnls <- NNLS_ext(Y = Y, Theta = true_data$X, alg = "pnnls", centralized_residual = FALSE)
R_known <- Y - X %*% P_known
R_pnnls <- Y - X %*% P_pnnls

mean(R_known)
mean(R_pnnls)


