## source files and libraries
source("s_sources.R")
set.seed(100)

## TODO: write a script to run for different values of p_sig/gamma_seq with only methods with intercepts
## set the parameters
n <- 20
m <- 1000
K <- 5
q <- 2
err <- TRUE
p_sig <- c(0.25, 0.5, 0.75)
lambda <- 5
# gamma <- 3
gamma_seq <- c(1/4, 1/2, 1, 2, 4)

## partition the canvas
# pdf("plots/simulation_with_binary_latent_factor_q_1_diff_p.pdf")
# par(mfrow = c(3, 1))
result_matrix <- matrix(ncol = 6)
## simulation functions
for (p in p_sig) {
  for (gamma in gamma_seq) {
    for (b in 1:B) {
      ## create a list to hold all the estimated proportions
      P_hat_ls <- list()
      true_data <- dSVA_model_sim_intercept(m, n, K, q, p, lambda, gamma, err = err, small_effect = FALSE)
      
      ## dSVA
      P_hat_ls$dSVA <- dsva_for_sim(Y = true_data$Y, Theta = true_data$X, n_comp = q)
      
      ## NNLS
      P_hat_ls$nnls <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "nnls", 
                                centralized_residual = FALSE)
      ## constrained NNLS
      P_hat_ls$pnnls <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "pnnls", 
                                 centralized_residual = FALSE)
      
      ## if we know both X and Z
      P_hat_ls$known <- apply(true_data$Y - true_data$Y_lat, 2, function(y) {lsei::pnnls(a = true_data$X, b = y , sum = 1)$x})
      
      ## mean sample-wise correlations
      cor_ls <- lapply(P_hat_ls, my_cor, P2 = true_data$P_star)
      
      ## mean CCC
      ccc_ls <- lapply(P_hat_ls, my_ccc, P2 = true_data$P_star)
      
      ## mean squared errors
      mse_ls <- lapply(P_hat_ls, my_mse, P2 = true_data$P_star)
      
      ## organize the results
      res_mat <- cbind(rbind(unlist(cor_ls), unlist(ccc_ls), unlist(mse_ls)), p, gamma)
      rownames(res_mat) <- c("cor", "ccc", "mse")
      result_matrix <- rbind(result_matrix, res_mat)
    }
  }
}
result_matrix <- result_matrix[-1, ]
result_df <- as_tibble(result_matrix, rownames = "metric")
result_df_long <- pivot_longer(result_df,
                               c("dSVA", "nnls", "pnnls", "known"),
                               names_to = "method",
                               values_to = "value")
result_df_long$method <- factor(result_df_long$method, levels = c("dSVA", "nnls", "pnnls", "known"))
p1 <- ggplot(result_df_long %>% filter(metric == "cor"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = paste("q =", q), x = "Method", y = "Pearson's correlation", caption = "Rows: p_sig; columns: gamma") +
  facet_grid(p ~ gamma) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45))

p2 <- ggplot(result_df_long %>% filter(metric == "ccc"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = paste("q =", q), x = "Method", y = "Concordance correlation coefficient", caption = "Rows: p_sig; columns: gamma") +
  facet_grid(p ~ gamma) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45))

p3 <- ggplot(result_df_long %>% filter(metric == "mse"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = paste("q =", q), x = "Method", y = "Mean squared error", caption = "Rows: p_sig; columns: gamma") +
  facet_grid(p ~ gamma) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45))

## use ggplot to plot the results
pdf("plots/simulation_with_binary_latent_factor_q_2_diff_p.pdf")
print(p1)
print(p2)
print(p3)
dev.off()
