library(MuSiC)
library(tidyverse)
# library(ggpubr)
library(gridExtra)
# library(devtools)
source("a_dSVA_functions.R")
source("f_customized_package_functions.R")

## read the bulk data
fadista_bulk <- readRDS("../dSVA_datasets/bulk_Fadista.rds")
fadista_bulk <- fadista_bulk[, !is.na(fadista_bulk$hba1c)]
fadista_bulk$group <- ifelse(fadista_bulk$hba1c > 5.7, "T2D", "Healthy") # use 5.7 as in Fan et al.

## make the column data
col_data <- tibble(subject = colnames(fadista_bulk),
                   group = fadista_bulk$group,
                   hba1c = fadista_bulk$hba1c,
                   age = fadista_bulk$age,
                   bmi = fadista_bulk$bmi,
                   tissue = fadista_bulk$tissue,
                   gender = fadista_bulk$gender)

## read the healthy sc data
cell_types <- c("acinar", "alpha", "beta", "delta", "ductal", "gamma")
healthy_sc <- readRDS("../dSVA_datasets/sc_healthy_Segerstolpe.rds")
healthy_sc <- healthy_sc[, healthy_sc$cellType %in% cell_types]
healthy_sc1 <- readRDS("../dSVA_datasets/EMTABsce_healthy.rds")

## plot the distribution of the cell-type expressions (in log)


## see the intersecting genes
d_genes <- intersect(rownames(healthy_sc), rownames(fadista_bulk))

## make a simple signature matrix (averaging over cell types)
avg_exp <- matrix(NA, nrow(healthy_sc), length(cell_types))
rownames(avg_exp) <- rownames(healthy_sc)
colnames(avg_exp) <- cell_types
for (i in seq(1, length(cell_types))) {
  exp_sub <- healthy_sc[, healthy_sc$cellType == cell_types[i]]
  avg_exp[, i] <- rowMeans(exprs(exp_sub))
}
avg_exp <- avg_exp[d_genes, ]
avg_exp <- avg_exp[rowSums(avg_exp) > 20, ] # 10554 features

# n_top = 1200
## choose the gene with the top variation
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
  # X_norm <- t(apply(X, 1, function(x) x/sqrt(sum(x^2))))
  # rowvars <- MatrixGenerics::rowVars(X_norm, useNames = TRUE)
  # top_names <- names(sort(-rowvars))[seq(n_top)]
  # X_top <- X[top_names, ]
  # X_top
}
healthy_sig_mat <- select_gene(X = avg_exp)
sig_mat <- healthy_sig_mat
bulk_mat <- exprs(fadista_bulk)[rownames(sig_mat), ]
# bulk_mat <- t(apply(bulk_mat, 1, function(y) y/sum(y)))
# all(rownames(sig_mat) == rownames(bulk_mat)) # TRUE

## cell mRNA amount
mRNA_amt <- vector(length = 6, mode = "double")
names(mRNA_amt) <- cell_types
for (i in seq(1, length(cell_types))) {
  exp_sub <- healthy_sc[, healthy_sc$cellType == cell_types[i]]
  mRNA_amt[i] <- mean(colSums(exprs(exp_sub)))
}
mRNA_amt_norm <- mRNA_amt / sum(mRNA_amt)
sig_mat2 <- apply(sig_mat, 2, function(x) x/sum(x)) %*% diag(mRNA_amt)
colnames(sig_mat2) <- colnames(sig_mat)
# sig_mat <- t(apply(sig_mat %*% diag(1/mRNA_amt), 1, function(x) x/sum(x)))

## plot the PCA on the residuals of the analysis
Y <- bulk_mat
X <- model.matrix(~1 + sig_mat)
B <- solve(t(X) %*% (X)) %*% t(X) %*% Y  
R <- Y - X %*% B 
pca <- prcomp(t(R), scale. = TRUE, center = TRUE)

plot(pca$sdev[1:20]^2, main = "Elbow Plot \nEigenvalues of R (Segerstolpe)", ylab = "Eigenvalue",
     log = "y")

autoplot(pca, data = col_data, colour = "group") +
  labs(title = "By Status", caption = "(T2D: Hba1c > 5.7)")

autoplot(pca, data = col_data, colour = "hba1c") +
  labs(title = "By Hba1c Value") +
  scale_color_viridis_b()

autoplot(pca, data = col_data, colour = "age") +
  labs(title = "By Age") +
  scale_color_viridis_b()

autoplot(pca, data = col_data, colour = "bmi") +
  labs(title = "By BMI") +
  scale_color_viridis_b()

autoplot(pca, data = col_data, colour = "gender") +
  labs(title = "By Gender") 


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
# 
# pca <- prcomp(t(R2), scale. = TRUE, center = TRUE)
# autoplot(pca, data = col_data, colour = "group") +
#   labs(title = "By Status", caption = "(T2D: Hba1c > 5.7)")
# 
# autoplot(pca, data = col_data, colour = "hba1c") +
#   labs(title = "By Hba1c Value") +
#   scale_color_viridis_b()
# 
# autoplot(pca, data = col_data, colour = "age") +
#   labs(title = "By Age") +
#   scale_color_viridis_b()
# 
# autoplot(pca, data = col_data, colour = "bmi") +
#   labs(title = "By BMI") +
#   scale_color_viridis_b()
# 
# autoplot(pca, data = col_data, colour = "gender") +
#   labs(title = "By Gender")

## plot the Ps
plot_P <- function(P) {
  P_wide <- as_tibble(P, rownames = "cell_type")
  P_long <- pivot_longer(data = P_wide,
                         cols = -cell_type,
                         names_to = "subject",
                         values_to = "proportion")
  P_all <- left_join(P_long, col_data, by = "subject")
  p <- ggplot(P_all, aes(x = cell_type, y = proportion, col = group, group_by = group)) +
    geom_boxplot()
  p
}

## perform dSVA-based deconvolution
## use the uncorrected signature matrix
perform_deconv <- function(bulk_mat, sig_mat, q_method, Title = NULL) {
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
  
  P_nnls <- apply(P_nnls, 2, function(x) x/sum(x))
  p <- plot_P(P_nnls)
  P_beta <- P_nnls["beta", ] 
  beta_df <- tibble(subject = names(P_beta),
                    prop = P_beta)
  all_df <- beta_df %>% 
    left_join(col_data, by = "subject")
  lm_obj <- lm(prop ~ hba1c + age + bmi + gender, data = all_df)
  q <- ggplot(all_df, aes(x = hba1c, y = prop, col = group)) +
    geom_point() +
    geom_abline(slope = lm_obj$coefficients[2],
                intercept = lm_obj$coefficients[1]) +
    annotate("text", x = 7, y = 0.2,
             label = paste0("p-value = ", round(anova(lm_obj)$`Pr(>F)`[1], 3))) +
    labs(x = "HbA1c",
         y = "Proportion of beta cells")
  p_all <- ggarrange(p, q, nrow = 1, ncol = 2, common.legend = TRUE,
                     legend = "bottom")
  p_all
}

p_be_un <- perform_deconv(bulk_mat, sig_mat, "be", "dSVA + NNLS \n(Uncorrected; BE q_hat = 6)")
p_cu_un <- perform_deconv(bulk_mat, sig_mat, "cutoff", "dSVA + NNLS \n(Uncorrected; CU q_hat = 1)")
p_no_un <- perform_deconv(bulk_mat, sig_mat, "no", "dSVA + NNLS \n(Uncorrected; No q_hat = 0)")


p_be_un
p_cu_un
p_no_un


p_be_co <- perform_deconv(bulk_mat, sig_mat2, "be", "dSVA + NNLS \n(Uncorrected; BE q_hat = 6)")
p_cu_co <- perform_deconv(bulk_mat, sig_mat2, "cutoff", "dSVA + NNLS \n(Uncorrected; CU q_hat = 1)")
p_no_co <- perform_deconv(bulk_mat, sig_mat2, "no", "dSVA + NNLS \n(Uncorrected; CU q_hat = 0)")

p_be_co
p_cu_co
p_no_co
# Title = "dSVA + NNLS \n(Uncorrected; BE q_hat = 6)"
## plot the correlation btw Hb1Ac and P_beta



# P_nnls <- dsva_for_sim(Y = bulk_mat, Theta = sig_mat, n_comp = 0, alg = "nnls", solver = "lsei")
## just use NNLS

# bulk.control.mtx <- exprs(fadista_bulk)[, fadista_bulk$group == "Healthy"]
# bulk.case.mtx <- exprs(fadista_bulk)[, fadista_bulk$group == "T2D"]


# est = music2_prop_t_statistics(bulk.control.mtx = bulk.control.mtx, 
#                                bulk.case.mtx = bulk.case.mtx, 
#                                sc.sce = healthy_sc1, 
#                                clusters = 'cellType', 
#                                samples = 'sampleID', 
#                                select.ct = c('acinar','alpha','beta','delta','ductal','gamma'), 
#                                n_resample=20, sample_prop=0.5,cutoff_c=0.05,cutoff_r=0.01)
# 
# est.prop = est$Est.prop


## Bulk RNA-seq data
# benchmark.eset = readRDS("../dSVA_datasets/bulk-eset.rds")
# benchmark.eset
# bulk.control.mtx = exprs(benchmark.eset)[, benchmark.eset$group == 'healthy']
# bulk.case.mtx = exprs(benchmark.eset)[, benchmark.eset$group == 't2d']
# set.seed(1234)
# est = my_music2_prop_t_statistics(bulk.control.mtx = bulk.control.mtx,
#                                   bulk.case.mtx = bulk.case.mtx,
#                                   sc.sce = healthy_sc1,
#                                   clusters = 'cellType',
#                                   samples = 'sampleID',
#                                   select.ct = c('acinar','alpha','beta','delta','ductal','gamma'),
#                                   n_resample=20, 
#                                   sample_prop=0.5, 
#                                   cutoff_c=0.05, 
#                                   cutoff_r=0.01)
# 
# est.prop <- est$Est.prop
# saveRDS(est.prop, "../dSVA_datasets/music_2_est.rds")
# col_data$sample_id <- as.character(fadista_bulk$sampleID)
# 
# P_wide <- as_tibble(t(est.prop), rownames = "cell_type")
# P_long <- pivot_longer(data = P_wide,
#                        cols = -cell_type,
#                        names_to = "sample_id",
#                        values_to = "proportion")
# P_all <- left_join(P_long, col_data, by = "sample_id")
# ggplot(P_all, aes(x = cell_type, y = proportion, col = group, group_by = group)) +
#   geom_boxplot()

## just run music as music2 is slow
# Est.prop = music_prop(bulk.mtx = exprs(fadista_bulk),
#                       sc.sce = healthy_sc1,
#                       clusters = 'cellType',
#                       samples = 'sampleID',
#                       select.ct = c('alpha', 'beta', 'delta', 'gamma', 'acinar', 'ductal'),
#                       verbose = F)
# P2 <- Est.prop$Est.prop.weighted
# P3 <- Est.prop$Est.prop.allgene
# 
# 
# P_wide <- as_tibble(t(P3), rownames = "cell_type")
# P_long <- pivot_longer(data = P_wide,
#                        cols = -cell_type,
#                        names_to = "subject",
#                        values_to = "proportion")
# P_all2 <- left_join(P_long, col_data, by = "subject")
# ggplot(P_all2, aes(x = cell_type, y = proportion, col = group, group_by = group)) +
#   geom_boxplot() +
#   labs(title = "MuSiC")
# Est.prop$Weight.gene
# Est.prop$Var.prop
