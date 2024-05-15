library(MuSiC)
library(tidyverse)
library(ggfortify)
# library(devtools)
source("a_dSVA_functions.R")
source("f_customized_package_functions.R")

## read the bulk data
fadista_bulk <- readRDS("../dSVA_datasets/bulk_Fadista.rds")
fadista_bulk <- fadista_bulk[, !is.na(fadista_bulk$hba1c)]
fadista_bulk$group <- ifelse(fadista_bulk$hba1c > 5.7, "T2D", "Healthy") # use 5.7 as in Fan et al.

## make the column data
col_data <- tibble(subject = colnames(fadista_bulk),
                   group = fadista_bulk$group)

## read the healthy sc data (Baron et al.)
cell_types <- c("acinar", "alpha", "beta", "delta", "ductal", "gamma")
healthy_sc <- readRDS("../dSVA_datasets/baron_count_raw.rds")
sample_orig <- str_split_i(colnames(healthy_sc), "_", 1)
healthy_sc_cell_types <- readRDS("../dSVA_datasets/baron_cell_types.rds")
colnames(healthy_sc) <- healthy_sc_cell_types

ind_keep <- colnames(healthy_sc) %in% cell_types
healthy_sc <- healthy_sc[, ind_keep]
sample_proc <- sample_orig[ind_keep]

## create a column data of original counts
col_data_sc <- tibble(subject = sample_proc,
                      cell_type = colnames(healthy_sc))
col_data_summary <- col_data_sc %>% 
  group_by(subject, cell_type) %>% 
  summarise(n = n()) %>% 
  ungroup(cell_type) %>% 
  mutate(prop = n/sum(n)) %>% 
  ungroup()

# col_data_summary %>% 
#   filter(subject == "human1") %>% 
#   select(prop) %>% 
#   unlist() %>% 
#   sum() # checking the proportion ranges

## summarizing the healthy ct proportions
p <- col_data_summary %>% 
  ggplot(aes(x = cell_type, y = prop, col = cell_type)) +
  geom_boxplot() +
  geom_jitter() +
  labs(title = "Major Cell-type Proportions in Healthy Subjects (Baron et al.)") 

ggsave("real_data/plots/baron_sc_prop_20231030.pdf", p)

## TODO: plot the distribution of the cell-type expressions (in log)

## see the intersecting genes
d_genes <- intersect(rownames(healthy_sc), rownames(fadista_bulk)) # 17336 features

## make a simple signature matrix (averaging over cell types)
avg_exp <- matrix(NA, nrow(healthy_sc), length(cell_types))
rownames(avg_exp) <- rownames(healthy_sc)
colnames(avg_exp) <- cell_types
for (i in seq(1, length(cell_types))) {
  exp_sub <- healthy_sc[, colnames(healthy_sc) == cell_types[i]]
  avg_exp[, i] <- rowMeans(exp_sub)
}
avg_exp <- avg_exp[d_genes, ]
avg_exp <- avg_exp[rowSums(avg_exp) > 1, ] # 4735 features

## choose the marker genes
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
healthy_sig_mat <- select_gene(X = avg_exp) # 1200 features
bulk_mat <- Biobase::exprs(fadista_bulk)[rownames(healthy_sig_mat), ]

## cell mRNA amount
mRNA_amt <- vector(length = 6, mode = "double")
names(mRNA_amt) <- cell_types
for (i in seq(1, length(cell_types))) {
  exp_sub <- healthy_sc[, colnames(healthy_sc) == cell_types[i]]
  mRNA_amt[i] <- mean(colSums(exp_sub))
}
mRNA_amt_norm <- mRNA_amt / sum(mRNA_amt)
sig_mat <- apply(healthy_sig_mat, 2, function(x) x/sum(x)) %*% diag(mRNA_amt_norm)  
colnames(sig_mat) <- colnames(healthy_sig_mat)
cor(sig_mat)

## plot the PCA on the residuals of the analysis
Y <- bulk_mat
X <- sig_mat
B <- solve(t(X) %*% (X)) %*% t(X) %*% Y  
R <- Y - X %*% B 
pca <- prcomp(t(R), scale. = TRUE, center = TRUE)
autoplot(pca, data = col_data, colour = "group")


## deconvolution with dSVA
q_hat <- estimate_n_comp(Y = bulk_mat, 
                         Theta = sig_mat, 
                         intercept = FALSE,
                         method = "cutoff",
                         B = 49, 
                         seed = 100) 
# q_hat2 <- estimate_n_comp(Y = bulk_mat, Theta = sig_mat, intercept = TRUE, method = "be", B = 49, seed = 100) # 2

# q_hat = 0
P_nnls <- dsva_for_sim(Y = bulk_mat, 
                       Theta = sig_mat, 
                       n_comp = q_hat, 
                       intercept = FALSE, 
                       alg = "nnls", 
                       solver = "lsei")

rownames(P_nnls) <- colnames(sig_mat)
colnames(P_nnls) <- colnames(bulk_mat)

P_nnls <- apply(P_nnls, 2, function(x) x/sum(x))

## plot the Ps

P_wide <- as_tibble(P_nnls, rownames = "cell_type")
P_long <- pivot_longer(data = P_wide,
                       cols = -cell_type,
                       names_to = "subject",
                       values_to = "proportion")
P_all <- left_join(P_long, col_data, by = "subject")
ggplot(P_all, aes(x = cell_type, y = proportion, col = group, group_by = group)) +
  geom_boxplot() + 
  labs(title = "dSVA + NNLS")
ggsave("real_data/plots/baron_dSVA_prop_20231030.pdf")

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
Est.prop = music_prop(bulk.mtx = exprs(fadista_bulk),
                      sc.sce = healthy_sc1,
                      clusters = 'cellType',
                      samples = 'sampleID',
                      select.ct = c('alpha', 'beta', 'delta', 'gamma', 'acinar', 'ductal'),
                      verbose = F)
P2 <- Est.prop$Est.prop.weighted
P3 <- Est.prop$Est.prop.allgene


P_wide <- as_tibble(t(P3), rownames = "cell_type")
P_long <- pivot_longer(data = P_wide,
                       cols = -cell_type,
                       names_to = "subject",
                       values_to = "proportion")
P_all2 <- left_join(P_long, col_data, by = "subject")
ggplot(P_all2, aes(x = cell_type, y = proportion, col = group, group_by = group)) +
  geom_boxplot() +
  labs(title = "MuSiC")
# Est.prop$Weight.gene
# Est.prop$Var.prop
