## source files and libraries
source("s_sources.R")
set.seed(100)

## set the parameters
n <- 20
m <- 1000
K <- 5
p <- 1
p_sig <- 0.5
lambda <- 5
# gamma <- 3
gamma_seq <- 1:3

## partition the canvas
pdf("plots/simulation_with_binary_latent_factor_q_1_no_error.pdf")
par(mfrow = c(3, 1))
result_matrix <- matrix(ncol = 11)
## simulation functions
for (gamma in gamma_seq) {
  for (b in 1:B) {
    ## create a list to hold all the estimated proportions
    P_hat_ls <- list()
    true_data <- dSVA_model_sim_intercept(m, n, K, p, p_sig, lambda, gamma)
    
    ## dSVA
    P_hat_ls$dSVA_no_int <- dsva_for_sim(Y = true_data$Y, Theta = true_data$X, n_comp = 1,
                                         intercept = FALSE)
    P_hat_ls$dSVA <- dsva_for_sim(Y = true_data$Y, Theta = true_data$X, n_comp = 1)
    
    ## NNLS
    P_hat_ls$nnls_no_int <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "nnls", 
                                     intercept = FALSE,
                                     centralized_residual = FALSE)
    P_hat_ls$nnls <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "nnls", 
                              centralized_residual = FALSE)
    P_hat_ls$nnls_no_int_ext <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "nnls", 
                                         intercept = FALSE)
    P_hat_ls$nnls_ext <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "nnls")
    
    ## constrained NNLS
    P_hat_ls$pnnls_no_int <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "pnnls", 
                                      intercept = FALSE,
                                      centralized_residual = FALSE)
    P_hat_ls$pnnls <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "pnnls", 
                               centralized_residual = FALSE)
    P_hat_ls$pnnls_no_int_ext <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "pnnls", 
                                          intercept = FALSE)
    P_hat_ls$pnnls_ext <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "pnnls")
    
    ## if we know both X and Z
    P_hat_ls$known <- apply(true_data$Y - true_data$Y_lat, 2, function(y) {lsei::pnnls(a = true_data$X, b = y , sum = 1)$x})
    
    ## mean sample-wise correlations
    cor_ls <- lapply(P_hat_ls, my_cor, P2 = true_data$P_star)
    
    ## mean CCC
    ccc_ls <- lapply(P_hat_ls, my_ccc, P2 = true_data$P_star)
  
    ## mean squared errors
    mse_ls <- lapply(P_hat_ls, my_mse, P2 = true_data$P_star)
    
    ## organize the results
    res_mat <- rbind(unlist(cor_ls), unlist(ccc_ls), unlist(mse_ls))
    rownames(res_mat) <- c("cor", "ccc", "mse")
    result_matrix <- rbind(result_matrix, res_mat)
  }
  result_matrix <- result_matrix[-1, ]
  result_df <- as_tibble(result_matrix, rownames = "metric")
  result_df_long <- pivot_longer(result_df,
                                 -metric,
                                 names_to = "method",
                                 values_to = "value")
  result_df_long$method <- factor(result_df_long$method, levels = colnames(result_matrix))
  boxplot(value ~ method,
          data = result_df_long, 
          subset = metric == "cor",
          main = paste0("gamma = ", gamma), xlab = "Method", ylab = "Pearson's correlation",
          col = c("violet", "orange", "azure", "lightpink"))
  #axis(1, at = as.character(seq(nlevels(result_df_long$method))), labels = c("dSVA", "NNLS", "PNNLS", "Known"))
  
  boxplot(value ~ method,
          data = result_df_long, 
          subset = metric == "ccc",
          main = paste0("gamma = ", gamma), xlab = "Method", ylab = "Concordance correlation coefficient",
          col = c("violet", "orange", "azure", "lightpink"))
  #axis(1, at = as.character(seq(nlevels(result_df_long$method))), labels = c("dSVA", "NNLS", "PNNLS", "Known"))
  
  boxplot(value ~ method,
          data = result_df_long, 
          subset = metric == "ccc",
          main = paste0("gamma = ", gamma), xlab = "Method", ylab = "Mean squared error",
          col = c("violet", "orange", "azure", "lightpink"))
  #axis(1, at = as.character(seq(nlevels(result_df_long$method))), labels = c("dSVA", "NNLS", "PNNLS", "Known"))
  # hist(q_hats, 
  #      # main = "Estimated Number of Latent Factors by Horn's PA",
  #      xlab = "Estimated q")
}
dev.off()