library(tidyverse)
library(gridExtra)
source("a_dSVA_functions.R")
source("f_customized_package_functions.R")

## CYTOF data
cytof_ma <- readRDS("../dSVA_datasets/microarray_CyTOF.rds")
# ciber_sc <- readRDS("../dSVA_datasets/ciber_GEP.rds")

ciber_sc <- readRDS("../dSVA_datasets/5mwn_GEP.rds")

d_genes <- intersect(rownames(cytof_ma), rownames(ciber_sc))
ciber_sc <- ciber_sc[d_genes, ]
ciber_sc_filtered <- ciber_sc[which(rowSums(ciber_sc) > 10), ]

## build the patient data
col_data <- readxl::read_xlsx("../dSVA_datasets/CyTOF.xlsx", skip = 1) %>% 
  dplyr::select(FileName, Group) %>% 
  dplyr::filter(FileName %in% colnames(cytof_ma)) %>% 
  dplyr::arrange(FileName) %>% 
  dplyr::rename(subject = FileName) %>% 
  mutate(subject = as.character(subject))

select_gene <- function(X, n_top = 200) {
  gene_list <- vector()
  X_sort <- t(apply(X, 1, sort, decreasing = TRUE))
  ratio <- X_sort[, 1]/X_sort[, 2] 
  max_ind <- apply(X, 1, function(x) which(x == max(x)))
  
  ## find the top marker genes
  for (i in 1:ncol(X)) {
    sub_ratio <- sort(ratio[which(max_ind == i)], decreasing = TRUE)
    gene_list <- c(gene_list, names(sub_ratio[1:n_top]))
  }
  
  X[gene_list, ]
}
ciber_sc_select <- select_gene(ciber_sc_filtered)
cytof_select <- cytof_ma[rownames(ciber_sc_select), ]

cytof_select <- cytof_select[, order(colnames(cytof_select))]

Y <- cytof_select
X <- model.matrix(~ 1 + ciber_sc_select)
B <- solve(t(X) %*% (X)) %*% t(X) %*% Y  
R <- Y - X %*% B 
pca <- prcomp(t(R), scale. = TRUE, center = TRUE)

## the elbow plot
plot(pca$sdev[1:20]^2, main = "Elbow Plot \nEigenvalues of R (CyTOF)", ylab = "Eigenvalue",
     log = "y")

autoplot(pca, data = col_data, colour = "Group")

## perform deconvolution
plot_P <- function(P) {
  P_wide <- as_tibble(P, rownames = "cell_type")
  P_long <- pivot_longer(data = P_wide,
                         cols = -cell_type,
                         names_to = "subject",
                         values_to = "proportion")
  P_all <- left_join(P_long, col_data, by = "subject")
  p <- ggplot(P_all, aes(x = cell_type, y = proportion, col = Group, group_by = Group)) +
    geom_boxplot()
  p
}

perform_deconv <- function(bulk_mat, sig_mat, q_method, Title = NULL) {
  # bulk_mat = cytof_select
  # sig_mat = ciber_sc_select
  # q_method = "be"
  if (q_method == "no") {
    q_hat <- 0
  } else {
    q_hat <- estimate_n_comp(Y = bulk_mat, Theta = sig_mat, intercept = TRUE, method = q_method, B = 49, seed = 100) # 2
    print(q_hat)
  }
  
  P_nnls <- dsva_for_sim(Y = bulk_mat, 
                         Theta = sig_mat, 
                         n_comp = q_hat, 
                         intercept = TRUE, 
                         alg = "nnls", 
                         solver = "lsei")
  
  rownames(P_nnls) <- colnames(sig_mat)
  colnames(P_nnls) <- colnames(bulk_mat)
  p <- plot_P(P_nnls) + labs(title = Title)
  p
}

p_be_un <- perform_deconv(cytof_select, ciber_sc_select, "be", "dSVA + NNLS \nBE: q_hat = 4")
p_cu_un <- perform_deconv(cytof_select, ciber_sc_select, "cutoff", "dSVA + NNLS \nCutoff: q_hat = 1")
p_no_un <- perform_deconv(cytof_select, ciber_sc_select, "no", "dSVA + NNLS \nNNLS: q_hat = 0")


p_be_un
p_cu_un
p_no_un

## what happens after dSVA correction
# svd_R <- svd(R)
# U_q <- as.matrix(svd_R$u[, 1]) 
# Psi_hat <- rbind(svd_R$d[1] %*% t(svd_R$v)[1, ])
# U_q <- svd_R$u[, 1]
# Psi_hat <- svd_R$d[1] * t(svd_R$v)[1, ]
# Psi_hat1 <- Psi_hat - mean(Psi_hat)
# C <- Psi_hat1 %*% t(Psi_hat1) 
# if (1 / C[1, 1] < 1e-15) message("q_hat = 1 and the inverse of Psi_hat %*% D_jn %*% t(Psi_hat) is smaller than 1e-15.")
# Gamma_hat <- U_q + (X %*% B %*% (Psi_hat1)) / C[1, 1]
# X_dSVA <- cbind(X, Gamma_hat)
# B2 <- solve(t(X_dSVA) %*% (X_dSVA)) %*% t(X_dSVA) %*% Y  
# R2 <- Y - X_dSVA %*% B2 
# pca <- prcomp(t(R2), scale. = TRUE, center = TRUE)
# autoplot(pca, data = col_data, colour = "Group")
