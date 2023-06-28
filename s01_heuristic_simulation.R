## simulation based on the dSVA paper. 
set.seed(100)
source("s02_simulation_functions.R")
library(extraDistr)
library(MASS)
library(dplyr)
library(tidyr)

## algorithm parameters
# m <- 10
# n <- 200
# K <- 5
# p <- 2
# p_sig <- 0.5
# p_class <- 0.5
# rho <- 0.5
# lambda <- 5
# mu <- 1 # confounding variable means
# # S <- 5
# B <- 100
# 
# result_matrix <- matrix(nrow = B, ncol = 8)
# 
# ## run the simulation B times
# for (b in 1:B) {
#   true_data <- dSVA_deconv_sim(m, n, K, p, p_sig, p_class, rho, lambda, mu)
#   P_dSVA <- dsva_ext(Y = true_data$Y, X = true_data$X, q = p)
#   P_nnls <- apply(true_data$Y, 2, function(y) {lsei::nnls(a = true_data$X, b = y)$x})
#   P_pnnls <- apply(true_data$Y, 2, function(y) {lsei::pnnls(a = true_data$X, b = y, sum = 1)$x})
#   
#   ## if we know both X and Z
#   X_new = cbind(true_data$Z, true_data$X)
#   P_hat <- apply(true_data$Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = p, sum = 1)$x})[-seq(p), ]
#   
#   cor_dSVA <- mean(my_cor(true_data$P_star, P_dSVA))
#   cor_nnls <- mean(my_cor(true_data$P_star, P_nnls))
#   cor_pnnls <- mean(my_cor(true_data$P_star, P_pnnls))
#   cor_best <- mean(my_cor(true_data$P_star, P_hat))
#   
#   mse_dSVA <- my_mse(true_data$P_star, P_dSVA)
#   mse_nnls <- my_mse(true_data$P_star, P_nnls)
#   mse_pnnls <- my_mse(true_data$P_star, P_pnnls)
#   mse_best <- my_mse(true_data$P_star, P_hat)
#   
#   result_matrix[b, ] <- c(cor_dSVA, cor_nnls, cor_pnnls, cor_best, mse_dSVA, mse_nnls, mse_pnnls, mse_best)
# }
# 
# result_df <- as_tibble(cbind(1:B, result_matrix))
# colnames(result_df) <- c("b", "cor_dSVA", "cor_nnls", "cor_pnnls", "cor_known", "mse_dSVA", "mse_nnls", "mse_pnnls", "mse_known")
# 
# par(mfrow = c(1, 2))
# boxplot(x = result_df$cor_dSVA, result_df$cor_nnls, result_df$cor_pnnls, result_df$cor_known, xlab = "Method", ylab = "Cor")
# axis(1, at = c("1", "2", "3", "4"), labels = c("dSVA", "NNLS", "PNNLS", "Known"))
# boxplot(result_df$mse_dSVA, result_df$mse_nnls, result_df$mse_pnnls, result_df$mse_known, xlab = "Method", ylab = "MSE")
# axis(1, at = c("1", "2", "3", "4"), labels = c("dSVA", "NNLS", "PNNLS", "Known"))

## plot the results
# result_df2 <- pivot_longer(data = result_df, cols = c("cor_dSVA", "cor_nnls", "cor_pnnls"), names_to = "cor", values_to = "cor_method")
# result_df3 <- pivot_longer(data = result_df2, cols = c("mse_dSVA", "mse_nnls", "mse_pnnls"), names_to = "mse", values_to = "mse_method")

## heuristic simulations for proportions

n <- 20
m <- 1000
K <- 5
p <- 2
p_sig <- 0.5
p_class <- 0.5
rho <- 0.5
lambda <- 5
mu <- 1 # confounding variable means
# S <- 5
B <- 100

result_matrix <- matrix(nrow = B, ncol = 6)
for (b in 1:B) {
  true_data <- dSVA_deconv_sim_people(m, n, K, p, p_sig, p_class, rho, lambda, mu)
  X_dSVA <- dsva_ext2(Y = true_data$Y, X = true_data$P, q = p)
  X_nnls <- apply(true_data$Y, 2, function(y) {lsei::nnls(a = true_data$P, b = y)$x})
  # X_pnnls <- apply(true_data$Y, 2, function(y) {lsei::pnnls(a = true_data$P, b = y, sum = 1)$x})
  
  ## if we know both X and Z
  P_new = cbind(true_data$Z, true_data$P)
  X_hat <- apply(true_data$Y, 2, function(y) {lsei::pnnls(a = P_new, b = y, k = p, sum = NULL)$x})[-seq(p), ]
  
  cor_dSVA <- mean(my_cor(t(true_data$X_star), t(X_dSVA)))
  cor_nnls <- mean(my_cor(t(true_data$X_star), t(X_nnls)))
  # cor_pnnls <- mean(my_cor(true_data$P_star, P_pnnls))
  cor_best <-  mean(my_cor(t(true_data$X_star), t(X_hat)))
  
  mse_dSVA <- my_mse(true_data$X_star, X_dSVA)
  mse_nnls <- my_mse(true_data$X_star, X_nnls)
  # mse_pnnls <- my_mse(true_data$P_star, P_pnnls)
  mse_best <- my_mse(true_data$X_star, X_hat)
  
  result_matrix[b, ] <- c(cor_dSVA, cor_nnls, cor_best, mse_dSVA, mse_nnls, mse_best)
}

result_df <- as_tibble(cbind(1:B, result_matrix))
colnames(result_df) <- c("b", "cor_dSVA", "cor_nnls", "cor_known", "mse_dSVA", "mse_nnls", "mse_known")

par(mfrow = c(1, 2))
boxplot(x = result_df$cor_dSVA, result_df$cor_nnls, result_df$cor_known, xlab = "Method", ylab = "Cor")
axis(1, at = c("1", "2", "3"), labels = c("dSVA", "NNLS", "Known"))
boxplot(result_df$mse_dSVA, result_df$mse_nnls, result_df$mse_known, xlab = "Method", ylab = "MSE")
axis(1, at = c("1", "2", "3"), labels = c("dSVA", "NNLS", "Known"))
