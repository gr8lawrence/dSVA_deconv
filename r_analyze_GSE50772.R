library(GEOquery)
library(GEOmetadb)
library(tidyverse)
library(ggfortify)
library(ggrepel)
library(gprofiler2)
library(ggthemes)
source("a_dSVA_functions.R")
source("r_real_data_funs.R")

# Get the raw data --------------------------------------------------------

# gse50772 <-  getGEO('GSE50772',GSEMatrix=TRUE)
# sm <- gse50772$GSE50772_series_matrix.txt.gz
# tibble(
#   sample = sm$geo_accession,
#   disease_status = sm$`disease status:ch1`,
# ) %>% 
#   write_csv("../dSVA_datasets/GSE50772_col.csv")
#
# 
# gse50772_tbl <-  getGEO("GSE50772",GSEMatrix=FALSE)
# # gse50772_tbl @gsms$GSM1228860 @ dataTable @table
# 
# raw_df_ls <- lapply(gse50772_tbl @ gsms, function(x) as_tibble(x @ dataTable @ table) )
# raw_df <- plyr::join_all(raw_df_ls, by = "ID_REF", type = "left")
# colnames(raw_df) <- c("probe", names(raw_df_ls))
# # head(raw_df_ls$GSM1228861)
# raw_df_final <- as_tibble(raw_df)
# write_csv(raw_df_final, "../dSVA_datasets/GSE50772_bulk.csv")

# Explore and pre-process data --------------------------------------------
bulk_df <- read_csv("../dSVA_datasets/GSE50772_bulk.csv")
bulk_col <- read_csv("../dSVA_datasets/GSE50772_col.csv")

bulk_mat_raw <- data.matrix(bulk_df[, -1])
rownames(bulk_mat_raw) <- bulk_df$probe


pca <- prcomp(t(bulk_mat_raw), scale. = TRUE, center = TRUE)
# pca$sdev
plot(pca$sdev^2/sum(pca$sdev^2), main = "GSE50772 Variation Explained")

autoplot(pca, data = bulk_col, colour = "disease_status") +
  labs(title = "PCA Plot on SLE Bulk Expression",
       caption = "GEO: GSE50722; platform: Affymetrix microarray") +
  geom_label_repel(aes(label = sample)) +
  theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Disease status", values = c("#ffd380", "#8a508f"))  


# Process the data to be annotated ----------------------------------------

## bulk
probe_name_df <- read_csv("../dSVA_datasets/Affy_probe_name.csv")
probe_convert_df <- tibble(probe = bulk_df$probe) %>% 
  left_join(probe_name_df, by = "probe") %>% 
  filter(!is.na(name))
bulk_raw_subset <- bulk_mat_raw[probe_convert_df$probe, ]
rownames(bulk_raw_subset) <- probe_convert_df$name
saveRDS(bulk_raw_subset, "../dSVA_datasets/GSE50772_subset_annotated.rds")

## signature
## ABIS-microarray
sig_mat <- read.delim("../dSVA_datasets/sigmatrixMicro.txt")
intersect_genes <- intersect(rownames(bulk_raw_subset), rownames(sig_mat))
Y_deconv <- bulk_raw_subset[intersect_genes, ]
Theta_deconv <- data.matrix(sig_mat[intersect_genes, ])


# Helper functions --------------------------------------------------------

get_long_form2 <- function(P_hat, bulk_cols) {
  P_df <- as_tibble(t(P_hat), rownames = "sample")
  P_long <- P_df %>% 
    pivot_longer(-1, 
                 names_to = "Cell_type",
                 values_to = "Proportion") %>% 
    left_join(bulk_cols, by = "sample")
  P_long
}

get_l2_CPM <- function(x, offset = FALSE) {
  log2(x/sum(x) * 10^6 + offset)
}

plot_result <- function(P_raw, P_adj, lab1 = "Unadjusted", lab2 = "Adjusted (dSVA)") {
  df_raw <- get_long_form2(P_raw, bulk_col) %>% mutate(adj = lab1)
  df_adj <- get_long_form2(P_adj, bulk_col) %>% mutate(adj = lab2)
  df_all <- rbind(df_raw, df_adj)
  p <- df_all %>% 
    ggplot(aes(x = Cell_type, y = Proportion, fill = disease_status)) +
    geom_boxplot() +
    theme_base(base_size = 12, base_family = "Times") +
    facet_wrap(~adj) +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_fill_manual(name = "Disease status", values = c("#ffd380", "#8a508f"))  
  p
}

get_R_pca <- function(Y, Theta, intercept = TRUE, bulk_col2 = bulk_col) {
  # Theta = Theta_deconv
  # intercept = TRUE
  # bulk_col2 = bulk_col[bulk_col$sample %in% colnames(Y_sub),]
  # Y = Y_sub[, bulk_col2$sample]
  
  ## add an intercept 
  if (intercept) {
    X <- model.matrix(~1 + Theta)
  } else {
    X <- Theta
  }
  
  ## perform first pass regression
  Bhat <- solve(t(X) %*% X) %*% t(X) %*% Y  
  R <- Y - X %*% Bhat
  
  ## plot pca
  pca2 <- prcomp(t(R), scale. = TRUE, center = TRUE)
  # pca$sdev
  plot(pca$sdev^2/sum(pca$sdev^2), main = "Bulk Residuals PCA Variation Explained",
       sub = "On a patient subset",
       ylab = "Variance explained by PC", xlab = "PC")
  
  p2 <- autoplot(pca2, data = bulk_col2, colour = "disease_status") +
    labs(title = "PCA Plot on SLE Bulk Residuals",
         caption = "GEO: GSE50722; platform: Affymetrix microarray") +
    geom_label_repel(aes(label = sample)) +
    ggthemes::theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_manual(name = "Disease status", values = c("#ffd380", "#8a508f")) 
  p2
}

# Deconvolution -----------------------------------------------------------

## get residual PCA
p_res <- get_R_pca(Y_deconv, Theta_deconv)

## raw count
P_ls <- run_dSVA(Y_deconv, Theta_deconv, intercept = TRUE, alg = "nnls", solver = "lsei", method = "none")
P_ls2 <- run_dSVA(Y_deconv, Theta_deconv, intercept = TRUE, alg = "nnls", solver = "lsei", method = "cutoff")

cor(Theta_deconv, P_ls2$Gamma_hat)

# P_nnls4 <- run_dSVA(Y_deconv, Theta_deconv, intercept = FALSE, alg = "nnls", solver = "lsei", q = 2)
# 
# P_nnls3 <- run_dSVA(Y_deconv, Theta_deconv, intercept = FALSE, alg = "nnls", solver = "lsei", method = "be", B = 19)
# P_nnls5 <- run_dSVA(Y_deconv, Theta_deconv, intercept = FALSE, alg = "nnls", solver = "lsei", method = "tw")
P_nnls <- P_ls$P_hat
P_nnls2 <- P_ls2$P_hat

plot_result(P_nnls, P_nnls2)
ggsave("plots/SLE/deconv_all_results.pdf")

# get_long_form2(P_nnls, bulk_col) %>% 
#   ggplot(aes(x = Cell_type, y = Proportion, fill = disease_status)) +
#   geom_boxplot()
# get_long_form2(P_nnls2, bulk_col) %>% 
#   ggplot(aes(x = Cell_type, y = Proportion, fill = disease_status)) +
#   geom_boxplot()
# get_long_form2(P_nnls5, bulk_col) %>% 
#   ggplot(aes(x = Cell_type, y = Proportion, fill = disease_status)) +
#   geom_boxplot()

## log2CPM

Y_log_cpm <- apply(Y_deconv, 2, get_l2_CPM, offset = TRUE)
Theta_log_cpm <- apply(Theta_deconv, 2, get_l2_CPM, offset = TRUE)

P_nnls6 <- run_dSVA(Y_log_cpm, Theta_log_cpm, intercept = TRUE, alg = "nnls", solver = "lsei", method = "none")
P_nnls7 <- run_dSVA(Y_log_cpm, Theta_log_cpm, intercept = TRUE, alg = "nnls", solver = "lsei", method = "cutoff")


get_long_form2(P_nnls6, bulk_col) %>% 
  ggplot(aes(x = Cell_type, y = Proportion, fill = disease_status)) +
  geom_boxplot()
get_long_form2(P_nnls7, bulk_col) %>% 
  ggplot(aes(x = Cell_type, y = Proportion, fill = disease_status)) +
  geom_boxplot()


# Deconvolution on a subset -----------------------------------------------

## get patient labels
control_id <- bulk_col$sample[bulk_col$disease_status == "Control"]
sle_id <- bulk_col$sample[bulk_col$disease_status == "SLE"]

## get the loadings of PC1 for SLE patients
lds <- pca$x[sle_id, "PC1"]
sle_id_keep <- names(lds[lds >= lds["GSM1228894"]])
all_id_keep <- c(control_id, sle_id_keep)

## perform deconvolution
Y_sub <- Y_deconv[, all_id_keep]
P_nnls_sub_1 <- run_dSVA(Y_sub, Theta_deconv, intercept = TRUE, alg = "nnls", solver = "lsei", method = "none")
P_nnls_sub_2 <- run_dSVA(Y_sub, Theta_deconv, intercept = TRUE, alg = "nnls", solver = "lsei", method = "cutoff")
plot_result(P_nnls_sub_1, P_nnls_sub_2)
ggsave("plots/SLE/deconv_sub_results.pdf")

P_nnls_sub_3 <- run_dSVA(Y_sub, Theta_deconv, intercept = TRUE, alg = "nnls", solver = "lsei", method = "be", B = 19)
P_nnls_sub_4 <- run_dSVA(Y_sub, Theta_deconv, intercept = TRUE, alg = "nnls", solver = "lsei", method = "tw")
plot_result(P_nnls_sub_3, P_nnls_sub_4, lab1 = "Adjusted (dSVA - be)", lab2 = "Adjusted (dSVA - tw)")
ggsave("plots/SLE/deconv_sub_results_be_tw.pdf")
