library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)

sc_DT <- fread("../dSVA_datasets/sc-dendritic/sc_signature_case-1024-lognormalized.txt")
unique(colnames(sc_DT))

table(colnames(sc_DT))

## plot the total expression
total_exp <- colSums(sc_DT[, -1])
total_exp_df <- tibble(cell_type = names(total_exp),
                       total_exp = total_exp)
total_exp_df$cell_type <- factor(total_exp_df$cell_type, levels = unique(colnames(sc_DT))
[-1])
total_exp_df %>% 
  ggplot(aes(x = cell_type, y = total_exp, col = cell_type)) +
  geom_boxplot() +
  scale_y_log10()
saveRDS(total_exp_df, "../dSVA_datasets/Chuwen_sc_total_exp.rds")

## gather the mean expressions
T_means <- sc_DT[, grepl("T cells", colnames(sc_DT)), with = FALSE] %>% rowMeans()
B_means <- sc_DT[, grepl("B cells", colnames(sc_DT)), with = FALSE] %>% rowMeans()
NK_means <- sc_DT[, grepl("NK cells", colnames(sc_DT)), with = FALSE] %>% rowMeans()
Mono_means <- sc_DT[, grepl("Monocytes", colnames(sc_DT)), with = FALSE] %>% rowMeans()
DC_means <- sc_DT[, grepl("Dendritic", colnames(sc_DT)), with = FALSE] %>% rowMeans()

sc_sig_full <- cbind(T_means, B_means, NK_means, Mono_means, DC_means)
rownames(sc_sig_full) <- sc_DT$genesymbol
colnames(sc_sig_full) <- c("T", "B", "NK", "Mono", "DC")

## save the object
saveRDS(sc_sig_full, "../dSVA_datasets/Chuwen_sc_sig_full.rds")

## Re-do it and only keep the middle 15-85 percentiles
get_middle_mean <- function(df, low = .15, up = .85) {
  # df = sc_DT[, grepl("T cells", colnames(sc_DT)), with = FALSE]
  # low = .15
  # up = .85
  dm <- data.matrix(df)
  col_sums <- colSums(dm)
  keep_inds <- which(col_sums >= quantile(col_sums, low) & col_sums <= quantile(col_sums, up))
  return(rowMeans(dm[, keep_inds]))
}

T_means_mid <- sc_DT[, grepl("T cells", colnames(sc_DT)), with = FALSE] %>% 
  get_middle_mean()
B_means_mid <- sc_DT[, grepl("B cells", colnames(sc_DT)), with = FALSE] %>% 
  get_middle_mean()
NK_means_mid <- sc_DT[, grepl("NK cells", colnames(sc_DT)), with = FALSE] %>% 
  get_middle_mean()
Mono_means_mid <- sc_DT[, grepl("Monocytes", colnames(sc_DT)), with = FALSE] %>% 
  get_middle_mean()
DC_means_mid <- sc_DT[, grepl("Dendritic", colnames(sc_DT)), with = FALSE] %>% 
  get_middle_mean()

sc_sig_mid <- cbind(T_means_mid, B_means_mid, NK_means_mid, Mono_means_mid, DC_means_mid)
rownames(sc_sig_mid) <- sc_DT$genesymbol
colnames(sc_sig_mid) <- c("T", "B", "NK", "Mono", "DC")

## save the object
saveRDS(sc_sig_mid, "../dSVA_datasets/Chuwen_sc_sig_mid.rds")

### version 3: use the sum of cellular expression
get_middle_sum <- function(df, low = .15, up = .85, seed = 1000, n_cell = 500) {
  set.seed(seed)
  dm <- data.matrix(df)
  col_sums <- colSums(dm)
  keep_inds <- which(col_sums >= quantile(col_sums, low) & col_sums <= quantile(col_sums, up))
  dm2 <- dm[, keep_inds]
  samp_inds <- sample(1:ncol(dm2), n_cell)
  return(rowSums(dm2[, samp_inds]))
}

T_sums_mid <- sc_DT[, grepl("T cells", colnames(sc_DT)), with = FALSE] %>% 
  get_middle_sum()
B_sums_mid <- sc_DT[, grepl("B cells", colnames(sc_DT)), with = FALSE] %>% 
  get_middle_sum()
NK_sums_mid <- sc_DT[, grepl("NK cells", colnames(sc_DT)), with = FALSE] %>% 
  get_middle_sum()
Mono_sums_mid <- sc_DT[, grepl("Monocytes", colnames(sc_DT)), with = FALSE] %>% 
  get_middle_sum()
DC_sums_mid <- sc_DT[, grepl("Dendritic", colnames(sc_DT)), with = FALSE] %>% 
  get_middle_sum()

sc_sig_sum_mid <- cbind(T_sums_mid, B_sums_mid, NK_sums_mid, Mono_sums_mid, DC_sums_mid)
rownames(sc_sig_sum_mid) <- sc_DT$genesymbol
colnames(sc_sig_sum_mid) <- c("T", "B", "NK", "Mono", "DC")

sum(rowSums(sc_sig_sum_mid) > 5) # 15,836 genes with total expression > 5 in 2,500 cells 

saveRDS(sc_sig_sum_mid, "../dSVA_datasets/Chuwen_sc_sig_sum_mid_500_cells.rds")

## make another signature matrix using only diseased cells
sc_DT_controls <- fread("../dSVA_datasets/sc-dendritic/sc_signature_control-1024-lognormalized.txt")
unique(colnames(sc_DT_controls))

table(colnames(sc_DT_controls))

T_sums_mid <- sc_DT_controls[, grepl("T cells", colnames(sc_DT_controls)), with = FALSE] %>% 
  get_middle_sum(n_cell = 250)
B_sums_mid <- sc_DT_controls[, grepl("B cells", colnames(sc_DT_controls)), with = FALSE] %>% 
  get_middle_sum(n_cell = 250)
NK_sums_mid <- sc_DT_controls[, grepl("NK cells", colnames(sc_DT_controls)), with = FALSE] %>% 
  get_middle_sum(n_cell = 250)
Mono_sums_mid <- sc_DT_controls[, grepl("Monocytes", colnames(sc_DT_controls)), with = FALSE] %>% 
  get_middle_sum(n_cell = 250)
DC_sums_mid <- sc_DT_controls[, grepl("Dendritic", colnames(sc_DT_controls)), with = FALSE] %>% 
  get_middle_sum(n_cell = 250)

sc_sig_sum_mid <- cbind(T_sums_mid, B_sums_mid, NK_sums_mid, Mono_sums_mid, DC_sums_mid)
rownames(sc_sig_sum_mid) <- sc_DT_controls$genesymbol
colnames(sc_sig_sum_mid) <- c("T", "B", "NK", "Mono", "DC")

sum(rowSums(sc_sig_sum_mid) > 5) # 13,781 genes with total expression > 5 in 2,500 cells 

saveRDS(sc_sig_sum_mid, "../dSVA_datasets/Chuwen_sc_sig_sum_mid_500_cells_control.rds")
