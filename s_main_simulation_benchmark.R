## source files and libraries
source("s_sources.R")
set.seed(100)

## TODO: write code for recording specific singular value distributions of true latent factors and residuals for q = 1 and q = 2.
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
n_sv <- 6 # the number of sv we want to plot
first_effect <- "miss"
second_effect <- "bin"
B <- 100

## parameters for debugging
# p_sig <- 0.75
# gamma <- 4

## partition the canvas
# pdf("plots/simulation_with_continuous_latent_factor_q_1_diff_p.pdf")
# par(mfrow = c(3, 1))

## list to hold bi-plots
# bi_Y_ls <- list()
# bi_Y_lat_ls <- list()
# bi_R_ls <- list()

result_matrix <- matrix(ncol = 8)
sv_results <- tibble(gamma = NaN, p = NaN, rank = NA, d_Y = NaN, d_Y_lat = NaN, d_R = NaN)
## simulation functions
for (p in p_sig) {
  for (gamma in gamma_seq) {
    # p = 0.25
    # gamma = 1/4
    for (b in 1:B) {
      # if (b == 1) {
      #   bi_ls <- get_bi_plot(true_data = true_data)
      #   bi_ls$bi_Y + labs(title = "Y")
      #   print(bi_ls$bi_Y_lat + labs(title = "Latent Factors"))
      #   print(bi_ls$bi_R + labs(title = "Fitted Residuals"))
      # }
      
      ## create a list to hold all the estimated proportions
      P_hat_ls <- list()
      true_data <- dSVA_model_sim_intercept(m, n, K, q, p, lambda, gamma, err = err, first_effect = first_effect, second_effect = second_effect)
      
      ## dSVA
      P_hat_ls$dSVA <- dsva_for_sim(Y = true_data$Y, Theta = true_data$X, n_comp = ifelse(first_effect == "me", K - 1, q))
      
      ## limma + PNNLS
      P_hat_ls$limma_pnnls <- limma_ext(Y = true_data$Y, Theta = true_data$X, batch = true_data$D[1, ])
      
      ## RUVr + PNNLS
      P_hat_ls$RUVr_pnnls <- RUVr_ext(Y = true_data$Y, Theta = true_data$X, n_comp = q)
      
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
      mae_ls <- lapply(P_hat_ls, my_mae, P2 = true_data$P_star)
      
      ## organize the results
      res_mat <- cbind(rbind(unlist(cor_ls), unlist(ccc_ls), unlist(mae_ls)), p, gamma)
      rownames(res_mat) <- c("cor", "ccc", "mae")
      result_matrix <- rbind(result_matrix, res_mat)
      
      ## for recording the distribution of eigenvalues
      mat <- get_sing_vals(true_data = true_data, n_sv = n_sv)
      mat2 <- cbind(gamma, p, 1:n_sv, t(mat))
      colnames(mat2) <- colnames(sv_results)
      sv_results <- rbind(sv_results, mat2)
    }
  }
}
result_matrix <- result_matrix[-1, ]
result_df <- as_tibble(result_matrix, rownames = "metric")
result_df_long <- pivot_longer(result_df,
                               c("dSVA", "limma_pnnls", "RUVr_pnnls", "nnls", "pnnls", "known"),
                               names_to = "method",
                               values_to = "value")
result_df_long$method <- factor(result_df_long$method, levels = c("dSVA", "limma_pnnls", "RUVr_pnnls", "nnls", "pnnls", "known"))
p1 <- ggplot(result_df_long %>% filter(metric == "cor"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = paste("q =", q), x = "Method", y = "Pearson's correlation", caption = "Rows: p_sig; columns: gamma") +
  facet_grid(p ~ gamma) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))

p2 <- ggplot(result_df_long %>% filter(metric == "ccc"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = paste("q =", q), x = "Method", y = "Concordance correlation coefficient", caption = "Rows: p_sig; columns: gamma") +
  facet_grid(p ~ gamma) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))

p3 <- ggplot(result_df_long %>% filter(metric == "mae"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = paste("q =", q), x = "Method", y = "Mean absolute error", caption = "Rows: p_sig; columns: gamma") +
  facet_grid(p ~ gamma) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))

## use ggplot to plot the results
if (q == 1) {
  pdf(paste0("plots/benchmark_sim_q_", q, "_", first_effect, "_v5.pdf"))
} else if (q == 2) {
  pdf(paste0("plots/benchmark_sim_q_", q, "_", first_effect, "_", second_effect, "_v5.pdf"))
}
print(p1)
print(p2)
print(p3)
dev.off()

## now plot the sv_s
sv_results <- sv_results[-1, ]
sv_results$rank <- factor(sv_results$rank, levels = as.character(1:n_sv))
sv_long <- pivot_longer(sv_results, cols = 4:6, names_to = "type", values_to = "sv")
sv_long$type <- factor(sv_long$type, 
                       levels = c("d_Y", "d_Y_lat", "d_R"),
                       labels = c("bulk (true + latent)", "latent", "fitted residuals"))

sv_p1 <- ggplot(sv_long, aes(x = rank, y = log10(sv), color = type)) +
  geom_boxplot() +
  facet_grid(p ~ gamma) +
  labs(title = paste("q =", q), x = "Rank", y = "log10(eigenvalue)") +
  theme(legend.position = "bottom")

sv_p2 <- ggplot(sv_long %>% filter(type == "bulk (true + latent)"), aes(x = rank, y = log10(sv))) +
  geom_boxplot() +
  facet_grid(p ~ gamma) +
  labs(title = "Sample-wise PCA on Y", x = "Rank", y = "log10(eigenvalue)") +
  theme(legend.position = "bottom")

sv_p3 <- ggplot(sv_long %>% filter(type == "latent"), aes(x = rank, y = log10(sv))) +
  geom_boxplot() +
  facet_grid(p ~ gamma) +
  labs(title = "Sample-wise PCA on Latent Factors", x = "Rank", y = "log10(eigenvalue)") +
  theme(legend.position = "bottom")

sv_p4 <- ggplot(sv_long %>% filter(type == "fitted residuals"), aes(x = rank, y = log10(sv))) +
  geom_boxplot() +
  facet_grid(p ~ gamma) +
  labs(title = "Sample-wise PCA on Fitted Residuals", x = "Rank", y = "log10(eigenvalue)") +
  theme(legend.position = "bottom")

if (q == 1) {
  pdf(paste0("plots/benchmark_sim_sample_PCA_q_", q, "_", first_effect, ".pdf"))
} else if (q == 2) {
  pdf(paste0("plots/benchmark_sim_sample_PCA_q_", q, "_", first_effect, "_", second_effect, ".pdf"))
}

print(sv_p1)
print(sv_p2)
print(sv_p3)
print(sv_p4)
dev.off()

# svd(Y)$d[1]/sum(svd(Y)$d)
# (svd(Y)$d[1])^2/sum(svd(Y)$d^2)


# pdf("plots/benchmark_pca_biplot_q_1_p=0.75_gamma=4.pdf")
# bi_ls <- get_bi_plot(true_data, orient = "gene")
# print(bi_ls$bi_Y)
# print(bi_ls$bi_Y_lat)
# print(bi_ls$bi_R)
# 
# bi_ls <- get_bi_plot(true_data, orient = "sample")
# print(bi_ls$bi_Y)
# print(bi_ls$bi_Y_lat)
# print(bi_ls$bi_R)
# dev.off()

## let's test out some PCA functions
# A = matrix(1:12, 3, 4)
# M_j = matrix(1/3, 3, 3)
# A_cent = (diag(1, 3, 3) - M_j) %*% A
# 
# prcomp(A)
# prcomp(A_cent) # same eigenvalues
