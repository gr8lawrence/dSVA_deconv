## source files and libraries
source("s_sources.R")
set.seed(100)

## run simulations when the proportion of genes affected in Theta2 is changed

## command line arguments (q and first_effect)
# cmd_args <- commandArgs(TRUE)
# q <- as.numeric(cmd_args[1])
# first_effect <- cmd_args[2]

n <- 20
m <- 500
K <- 4
q <- 1
B <- 10

# n <- 40
# m <- 2000
# K <- 10
# q <- 1
err <- TRUE
p_sig <- c(0.1, 0.25, 0.5)
# p_change <- c(0, 0.01, 0.05, 0.1, 0.25)
p_change <- c(0, 0.25, 0.5, 0.75, 1)
lambda <- 5
n_sv <- 6 # the number of sv we want to plot
# B <- 100
plot_sv <- FALSE # whether to plot the singular/eigenvalues
mod <- NULL # ComBat covariate

## set up a color palette
my_palette <- c("#ffd380", "#ff8531", "#ff6361", "#bc5090", "#8a508f", "#2c4875")

## parameters for debugging
# p = 0.25
# gamma = 1/4
# p_sig <- 0.75

## partition the canvas
# pdf("plots/simulation_with_continuous_latent_factor_q_1_diff_p.pdf")
# par(mfrow = c(3, 1))

## list to hold bi-plots
# bi_Y_ls <- list()
# bi_Y_lat_ls <- list()
# bi_R_ls <- list()

## make title
title = paste0("Cell Type with Missing Gene Signatures", " (q = ", q, ")")
fig_title = paste0(title, " (q = ", q, ")")
fig_subtitle = "pNNLS Solver"

## specify batches
batch <- c(rep(0, floor(n/2)), rep(1, ceiling(n/2)))
result_matrix <- matrix(ncol = 8) # ncol = 8 as of 01/19/2024

## initiate matrices for estimated q_hat
q_hat_mat <- matrix(NaN, nrow = length(p_sig), ncol = length(p_change))
rownames(q_hat_mat) <- p_sig
colnames(q_hat_mat) <- p_change

if (plot_sv) sv_results <- tibble(gamma = NaN, p = NaN, rank = NA, d_Y = NaN, d_Y_lat = NaN, d_R = NaN)
## simulation functions
for (p in p_sig) {
  for (p_c in p_change) {
    for (b in 1:B) {
      # p = .1
      # b = 1
 
      ## create a list to hold all the estimated proportions
      P_hat_ls <- list()
      true_data <- dSVA_model_sim_pert(m, n, K, q, p, lambda, err = err, p_pert = p_pert)
      
      ## estimate n_comp once for each setting
      if (b == 1) dSVA_n_comp <- estimate_n_comp(Y = true_data$Y, Theta = true_data$X, method = "be", B = 49, seed = 100)
      
      ## dSVA + PNNLS
      # P_hat_ls$dSVA <- dsva_for_sim(Y = true_data$Y, Theta = true_data$X, n_comp = ifelse(first_effect == "me", K - 1, q))
      P_hat_ls$dSVA_pnnls <- dsva_for_sim(Y = true_data$Y, Theta = true_data$X, n_comp = dSVA_n_comp, alg = "pnnls", solver = "lsei") # using the estimated components

      ## limma + PNNLS
      P_hat_ls$limma_pnnls <- limma_ext(Y = true_data$Y, Theta = true_data$X, batch = batch, alg = "pnnls")
      
      ## change the covariates if q = 2 and first_effect == con
      # if (q == 2 & first_effect == "con") mod <- matrix(true_data$D[1, ], n, 1)
      
      ## ComBat + PNNLS
      P_hat_ls$ComBat_pnnls <- combat_ext(Y = true_data$Y, Theta = true_data$X, batch = batch, combat_seq = FALSE, alg = "pnnls", mod = mod)
      
      ## Combat_seq + PNNLS
      P_hat_ls$ComBat_seq_pnnls <- combat_ext(Y = true_data$Y, Theta = true_data$X, batch = batch, combat_seq = TRUE, alg = "pnnls", covar_mod = mod)
      
      ## PNNLS (sum-to-one NNLS)
      P_hat_ls$pnnls <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "pnnls", centralized_residual = FALSE)
      
      ## if we know the underlying truths
      # Theta_full <- cbind(true_data$X, true_data$Z)
      # P_hat_ls$known <- NNLS_ext(Y = true_data$Y, Theta = Theta_full, alg = "pnnls", centralized_residual = FALSE)
      
      P_hat_ls$known <- get_p_known(true_data, "cc", TRUE, "pnnls")
       
      # if (first_effect == "cc") {
      #   P_hat_ls$known <- cbind(
      #     apply(true_data$Y[, true_data$D[1, ] == 0], 2, function(y) {lsei::pnnls(a = true_data$X, b = y , sum = 1)$x}),
      #     apply(true_data$Y[, true_data$D[1, ] == 1], 2, function(y) {lsei::pnnls(a = true_data$X + true_data$X2, b = y , sum = 1)$x})
      #   )
      # } else {
      #   P_hat_ls$known <- apply(true_data$X %*% true_data$P_star + true_data$E, 2, function(y) {lsei::pnnls(a = true_data$X, b = y , sum = 1)$x})
      # }

      ## mean sample-wise correlations
      cor_ls <- lapply(P_hat_ls, my_cor, P2 = true_data$P_star)
      
      ## mean CCC
      ccc_ls <- lapply(P_hat_ls, my_ccc, P2 = true_data$P_star)
      
      ## mean squared errors
      mae_ls <- lapply(P_hat_ls, my_mae, P2 = true_data$P_star)
      
      ## organize the results
      res_mat <- cbind(rbind(unlist(cor_ls), unlist(ccc_ls), unlist(mae_ls)), p, p_c)
      rownames(res_mat) <- c("cor", "ccc", "mae")
      result_matrix <- rbind(result_matrix, res_mat)
      
      ## for recording the distribution of eigenvalues
      if (plot_sv) {
        mat <- get_sing_vals(true_data = true_data, n_sv = n_sv)
        mat2 <- cbind(gamma, p, 1:n_sv, t(mat))
        colnames(mat2) <- colnames(sv_results)
        sv_results <- rbind(sv_results, mat2)
      }
    }
  }
}
result_matrix <- result_matrix[-1, ]
result_df <- as_tibble(result_matrix, rownames = "metric")
result_df_long <- pivot_longer(result_df,
                               c("dSVA_pnnls", "limma_pnnls", "ComBat_pnnls", "ComBat_seq_pnnls", "pnnls", "known"),
                               names_to = "method",
                               values_to = "value")

## labeling the results
result_df_long$method <- factor(result_df_long$method, 
                                levels = c("dSVA_pnnls", "limma_pnnls", "ComBat_pnnls", "ComBat_seq_pnnls", "pnnls", "known"),
                                labels = c("dSVA", "limma", "ComBat", "ComBat-seq", "No adjustment", "No latent effects"))
                                # labels = c("dSVA + pNNLS", "limma + pNNLS", "ComBat + pNNLS", "ComBat-seq + pNNLS", "pNNLS", "Known Z + pNNLS"))

result_df_long$p2 <- factor(result_df_long$p,
                            levels = as.character(p_sig),
                            labels = as.expression(parse(text = paste("p[sig]","==", p_sig, sep="")))
)
result_df_long$p_c2 <- factor(result_df_long$p_c,
                                levels = as.character(p_change),
                                labels = as.expression(parse(text = paste("p[pert]","==", p_change, sep=""))))

p1 <- ggplot(result_df_long %>% filter(metric == "cor"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = fig_title,
       x = "Method", y = "Pearson's correlation") +
  facet_grid(p2 ~ p_c2, labeller = label_parsed, scales = "free") +
  theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_manual(values = my_palette, name = "Method") 

p2 <- ggplot(result_df_long %>% filter(metric == "ccc"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = fig_title,
       x = "Method", y = "Concordance correlation coefficient") +
  facet_grid(p2 ~ p_c2, labeller = label_parsed, scales = "free") +
  theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_manual(values = my_palette, name = "Method")

p3 <- ggplot(result_df_long %>% filter(metric == "mae"), aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  labs(title = fig_title,
       x = "Method", y = "Mean absolute error") +
  facet_grid(p2 ~ p_c2, labeller = label_parsed, scales = "free") +
  theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_manual(values = my_palette, name = "Method")

## use ggplot to plot the results

pdf(paste0("plots/cc_COVAR_sim_q_", q, "_", "cc_20240229_pNNLS_be.pdf"))
print(p1)
print(p2)
print(p3)
dev.off()

## print the q_hat_mat
print(q_hat_mat)
