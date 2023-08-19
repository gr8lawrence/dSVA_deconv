## source files and libraries
source("s_sources.R")
set.seed(100)

## TODO: write code to perform the simulation of Bingxin

n <- 20
m <- 1000
K <- 5
p_sig <- c(0.25, 0.5, 0.75)
chi_df <- 200
d_seq <- seq(0, 0.2, 0.05)
n_sv <- 6
B <- 20
err <- TRUE
result_matrix <- matrix(ncol = 7)
## first figure out how many latent factors there are
## supposedly 19
q = 19
## the main simulation body
for (p in p_sig) {
  for (d in d_seq) {
    for (b in 1:B) {
      P_hat_ls <- list() # create a list to hold results
      true_data <- get_Y_bz(m, n, K, p, chi_df, d, err)
      
      ## dSVA
      P_hat_ls$dSVA <- dsva_for_sim(Y = true_data$Y, Theta = true_data$X1, n_comp = 19)
      
      ## limma + PNNLS
      # P_hat_ls$limma_pnnls <- limma_ext(Y = true_data$Y, Theta = true_data$X1, batch = true_data$D[1, ])
      
      ## RUVr + PNNLS
      P_hat_ls$RUVr_pnnls <- RUVr_ext(Y = true_data$Y, Theta = true_data$X1, n_comp = q)
      
      ## NNLS
      P_hat_ls$nnls <- NNLS_ext(Y = true_data$Y, Theta = true_data$X1, alg = "nnls", 
                                centralized_residual = FALSE)
      ## constrained NNLS
      P_hat_ls$pnnls <- NNLS_ext(Y = true_data$Y, Theta = true_data$X1, alg = "pnnls", 
                                 centralized_residual = FALSE)
      
      ## if we know the real X instead of X1
      P_hat_ls$known <- apply(true_data$Y, 2, function(y) {lsei::pnnls(a = true_data$X, b = y , sum = 1)$x})
      
      # print("Here")
      
      ## mean sample-wise correlations
      cor_ls <- lapply(P_hat_ls, my_cor, P2 = true_data$P_star)
      
      ## mean CCC
      ccc_ls <- lapply(P_hat_ls, my_ccc, P2 = true_data$P_star)
      
      ## mean squared errors
      mse_ls <- lapply(P_hat_ls, my_mse, P2 = true_data$P_star)
      
      ## organize the results
      res_mat <- cbind(rbind(unlist(cor_ls), unlist(ccc_ls), unlist(mse_ls)), p, d)
      rownames(res_mat) <- c("cor", "ccc", "mse")
      result_matrix <- rbind(result_matrix, res_mat)
    }
  }
}

result_matrix <- result_matrix[-1, ]
result_df <- as_tibble(result_matrix, rownames = "metric")
result_df_long <- pivot_longer(result_df,
                               c("dSVA", "RUVr_pnnls", "nnls", "pnnls", "known"),
                               names_to = "method",
                               values_to = "value")
result_df_long$method <- factor(result_df_long$method, levels = c("dSVA", "RUVr_pnnls", "nnls", "pnnls", "known"))
p1 <- ggplot(result_df_long %>% filter(metric == "cor"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = paste("q =", q), x = "Method", y = "Pearson's correlation", caption = "Rows: p_sig; columns: gamma") +
  facet_grid(p ~ d) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))

p2 <- ggplot(result_df_long %>% filter(metric == "ccc"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = paste("q =", q), x = "Method", y = "Concordance correlation coefficient", caption = "Rows: p_sig; columns: gamma") +
  facet_grid(p ~ d) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))

p3 <- ggplot(result_df_long %>% filter(metric == "mse"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = paste("q =", q), x = "Method", y = "Mean squared error", caption = "Rows: p_sig; columns: gamma") +
  facet_grid(p ~ d) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))

## use ggplot to plot the results
pdf("plots/benchmark_sim_q_1_dropouts.pdf")
print(p1)
print(p2)
print(p3)
dev.off()

