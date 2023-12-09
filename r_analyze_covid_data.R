library(tidyverse)
library(ggpubr)
source("a_dSVA_functions.R")
source("f_customized_package_functions.R")


# Chuwen Data -------------------------------------------------------------


controls_df <- read_delim("../dSVA_datasets/chuwen_bulk_controls.txt")
cases_df <- read_delim("../dSVA_datasets/chuwen_bulk_cases.txt")
sig_df <- read_delim("../dSVA_datasets/chuwen_sig_controls.txt")

## process data
cons <- data.matrix(controls_df[, -1])
cases <- data.matrix(cases_df[, -1])
both <- cbind(cons, cases)
rownames(both) <- controls_df$genesymbol
covid_cols <- readRDS("../dSVA_datasets/covid_col.rds")

sig <- data.matrix(sig_df[, -1])
rownames(sig) <- sig_df$NAME
genes <- intersect(rownames(both), rownames(sig))

Y <- both[genes, ]
Theta <- sig[genes, ]

q_hat <- estimate_n_comp(Y = Y, Theta = Theta, intercept = TRUE, method = "be", B = 49, seed = 100) # 5
P_nnls <- dsva_for_sim(Y = Y, 
                       Theta = Theta, 
                       n_comp = q_hat, 
                       intercept = TRUE, 
                       alg = "nnls", 
                       solver = "lsei")

rownames(P_nnls) <- colnames(Theta)
colnames(P_nnls) <- colnames(Y)

P_nnls <- apply(P_nnls, 2, function(x) x/sum(x))

P_nnls2 <- dsva_for_sim(Y = Y, 
                       Theta = Theta, 
                       n_comp = 1, 
                       intercept = TRUE, 
                       alg = "nnls", 
                       solver = "lsei")

rownames(P_nnls2) <- colnames(Theta)
colnames(P_nnls2) <- colnames(Y)

P_nnls2 <- apply(P_nnls2, 2, function(x) x/sum(x))


# GSE64655 Signature Matrix -----------------------------------------------


## Now, let's analyze it using a different reference (GSE 64655)
bulk_full1 <- readRDS("../dSVA_datasets/covid_bulk_count_processed.rds")
# sig_full <- readRDS("../dSVA_datasets/GSE64655_sorted_count_processed.rds")
sig_full <- readRDS("../dSVA_datasets/GSE64655_sorted_count_0d_processed.rds")

## remove the two suspect samples 
bulk_full <- bulk_full1[, !grepl("S155", colnames(bulk_full1)) & !grepl("S179", colnames(bulk_full1))]

sig_full <- sig_full[, 1:4] # remove DC cells

bulk_cols <- readRDS("../dSVA_datasets/covid_col.rds")
bulk_cols$severity <- factor(bulk_cols$severity, levels = c("Healthy",
                                                            "Convalescent",
                                                            "Moderate",
                                                            "Severe",
                                                            "ICU"))
sig_full2 <- sig_full[rowSums(sig_full) > 10, ]

## find overlap genes
overlap_genes <- intersect(rownames(bulk_full), rownames(sig_full2)) # 15,466 genes
bulk_sub1 <- bulk_full[overlap_genes, ]
sig_sub1 <- sig_full2[overlap_genes, ]

## use ComBat-seq to remove batch effects
# all_sub <- cbind(bulk_sub1, sig_sub1)
# batches <- c(rep(1, 32), rep(2, 4))
# all_combat <- sva::ComBat_seq(all_sub, batch = batches)
# bulk_sub1 <- all_combat[, 1:32]
# sig_sub1 <- all_combat[, -(1:32)]

# saveRDS(bulk_sub1, "../dSVA_datasets/bulk_sub_GSE64655.rds")
# saveRDS(sig_sub1, "../dSVA_datasets/sig_sub_GSE64655.rds")

# saveRDS(bulk_sub1, "../dSVA_datasets/bulk_sub_ComBat_GSE64655.rds")
# saveRDS(sig_sub1, "../dSVA_datasets/sig_sub_ComBat_GSE64655.rds")

sig_sub1[grepl("CD3", rownames(sig_sub1)), ] # high on T cell
sig_sub1[grepl("CD4", rownames(sig_sub1)), ] # high on monocyte (not normal)
sig_sub1[grepl("CD8", rownames(sig_sub1)), ] # high on T cell
sig_sub1[grepl("CD14", rownames(sig_sub1)), ] # high on monocyte

sig_sub1[grepl("CD38", rownames(sig_sub1)), ] # high on NK (should be on B?)
sig_sub1[grepl("CD19", rownames(sig_sub1)), ] # high on B
sig_sub1[grepl("CD20", rownames(sig_sub1)), ] # high on B


## find the marker genes (100 per cell type)
rat_df <- tibble(gene = rownames(sig_sub1),
                 cell_type = apply(sig_sub1, 1, function(x) colnames(sig_sub1)[which(x == max(x))[1]]))
max_val <- apply(sig_sub1, 1, max)
sec_val <- apply(sig_sub1, 1, function(x) sort(x, decreasing = TRUE)[2])
rat_df <- rat_df %>% 
  mutate(ratio = max_val/sec_val)
marker_genes <- rat_df %>% 
  group_by(cell_type) %>% 
  arrange(desc(ratio)) %>% 
  filter(row_number() <= 300) %>% 
  ungroup %>% 
  dplyr::select(gene) %>% 
  unlist

bulk_marker <- bulk_sub1[marker_genes, ]
sig_marker <- sig_sub1[marker_genes, ]
cor(sig_marker)
# sig_marker2 <- sig_marker[, !grepl("DC", colnames(sig_marker))]

## run dSVA deconv
q_hat2 <- estimate_n_comp(Y = bulk_marker, Theta = sig_marker, intercept = TRUE, method = "be", B = 49, seed = 100)
# q_hat2 <- estimate_n_comp(Y = bulk_marker, Theta = sig_marker, intercept = TRUE, method = "cutoff")

# q_hat2 <- estimate_n_comp(Y = bulk_marker, Theta = sig_marker, intercept = FALSE, method = "be", B = 49, seed = 100)

P_nnls2 <- dsva_for_sim(Y = bulk_marker, 
                        Theta = sig_marker, 
                        n_comp = q_hat2, 
                        intercept = TRUE, 
                        alg = "pnnls", 
                        solver = "lsei") # using pnnls the results get interesting

rownames(P_nnls2) <- colnames(sig_marker)
colnames(P_nnls2) <- colnames(bulk_marker)

P_df <- as_tibble(t(P_nnls2), rownames = "subject_id")
P_long <- P_df %>% 
  pivot_longer(-subject_id, 
               names_to = "Cell_type",
               values_to = "Proportion") %>% 
  left_join(bulk_cols, by = "subject_id")
P_long %>% 
  ggplot(aes(x = Cell_type, y = Proportion, col = disease_state)) +
  geom_boxplot() +
  labs(title = "m = 1600 (Day 0 Reference, no DC, \n2 Sample Removed)", x = "Cell type") +
  ggpubr::theme_pubr()

P_long %>% 
  filter(Cell_type == "T") %>% 
  ggplot(aes(x = severity, y = Proportion, col = severity)) +
  geom_boxplot(alpha = 0.4) +
  geom_jitter() +
  labs(title = "T-Cell Proportion (m = 1600, Day 0 Reference, No DC, \n2 Sample Removed)", x = "Disease Severity") +
  scale_color_viridis_d() +
  ggpubr::theme_pubr() 

P_long %>% 
  filter(Cell_type == "Mono") %>% 
  ggplot(aes(x = severity, y = Proportion, col = severity)) +
  geom_boxplot(alpha = 0.4) +
  geom_jitter() +
  labs(title = "Monocyte Proportion (m = 1600, Day 0 Reference, No DC, \n2 Sample Removed)", x = "Disease Severity") +
  scale_color_viridis_d() +
  ggpubr::theme_pubr()

## compare healthy and covid-19 patients before and after the correction
P_null <- dsva_for_sim(Y = bulk_marker, 
                       Theta = sig_marker, 
                       n_comp = 0, 
                       intercept = TRUE, 
                       alg = "pnnls", 
                       solver = "lsei") # using pnnls the results get interesting

P_adjust <- dsva_for_sim(Y = bulk_marker, 
                         Theta = sig_marker, 
                         n_comp = q_hat2, 
                         intercept = TRUE, 
                         alg = "pnnls", 
                         solver = "lsei")

rownames(P_null) <- rownames(P_adjust) <- colnames(sig_marker)
colnames(P_null) <- colnames(P_adjust) <- colnames(bulk_marker)

saveRDS(P_adjust, "../dSVA_datasets/dSVA_Padj_q_3.rds")
saveRDS(P_null, "../dSVA_datasets/dSVA_Pnull_q_0.rds")

get_long_form <- function(P_hat) {
  P_df <- as_tibble(t(P_hat), rownames = "subject_id")
  P_long <- P_df %>% 
    pivot_longer(-subject_id, 
                 names_to = "Cell_type",
                 values_to = "Proportion") %>% 
    left_join(bulk_cols, by = "subject_id")
  P_long
}

P_null_long <- get_long_form(P_null) %>% 
  rename(Proportion_unadjusted = Proportion)
P_adjust_long <- get_long_form(P_adjust) %>% 
  rename(Proportion_adjusted = Proportion)

P_all <- tibble(subject_id = P_null_long$subject_id,
                prop_unadj = P_null_long$Proportion_unadjusted,
                prop_adj = P_adjust_long$Proportion_adjusted,
                severity = P_null_long$severity,
                disease_state = P_null_long$disease_state,
                cell_type = P_null_long$Cell_type)


# P_all <- P_null_long %>% 
#   left_join(P_adjust_long %>% dplyr::select(subject_id, Proportion_adjusted), by = "subject_id")

ggplot(P_all, aes(x = prop_unadj, y = prop_adj, col = severity)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_fixed(ratio = 1) +
  labs(x = "Unadjusted proportion (q = 0)", y = paste0("Adjusted proportion (q = ", q_hat2, ")")) +
  facet_wrap(~cell_type) +
  ggpubr::theme_pubr()



## can't know lymphopenia bc not whole blood
# P_long %>% 
#   filter(Cell_type == "T" | Cell_type == "B" | Cell_type == "NK") %>%
#   group_by(subject_id) %>% 
#   summarise(Lymph_prop = sum(Proportion)) %>% 
#   left_join(col_data, by = "subject_id") %>% 
#   ggplot(aes(x = disease_state, y = Lymph_prop, col = disease_state)) +
#   geom_boxplot() +
#   labs(title = "m = 50", x = "Cell type")


# Chuwen's matrix ---------------------------------------------------------

## use the controls' scRNA data from Chuwen to do this again
# chuwen_full <- readRDS("../dSVA_datasets/Chuwen_sc_sig_full.rds") # signature using all cells
# chuwen_full <- readRDS("../dSVA_datasets/Chuwen_sc_sig_mid.rds") # signature using cells whose total counts are in the middle
chuwen_full <- readRDS("../dSVA_datasets/Chuwen_sc_sig_sum_mid_500_cells.rds") # signature using sum of expression from cells whose total counts are in the middle
# chuwen_full <- readRDS("../dSVA_datasets/Chuwen_sc_sig_sum_mid_500_cells_control.rds") # CONTROLS: healthy sample signature using sum of expression from cells whose total counts are in the middle

chuwen_full <- chuwen_full[, 1:4] # remove DCs

chuwen_full2 <- chuwen_full[rowSums(chuwen_full) > 5, ]
overlap_genes2 <- intersect(rownames(bulk_full), rownames(chuwen_full2)) # 1,161 genes
bulk_sub2 <- bulk_full[overlap_genes2, ]
sig_sub2 <- chuwen_full2[overlap_genes2, ]

## check the signature matrix further
dim(sig_sub2)
dim(sig_full2)

cor(sig_sub2[grepl("CD4", rownames(sig_sub2)), "T"], sig_sub2[grepl("CD4", rownames(sig_sub2)), "NK"])

sig_sub2[grepl("CD3", rownames(sig_sub2)), ] # high on T cell
sig_sub2[grepl("CD4", rownames(sig_sub2)), ] # high on monocyte (not normal)
sig_sub2[grepl("CD8", rownames(sig_sub2)), ] # high on T cell
sig_sub2[grepl("CD14", rownames(sig_sub2)), ] # high on monocyte

sig_sub2[grepl("CD38", rownames(sig_sub2)), ] # high on NK in the healthy reference (should be on B?) - High on B in the COVID-19 reference!
sig_sub2[grepl("CD19", rownames(sig_sub2)), ] # high on B 


## find the marker genes (100 per cell type)
rat_df2 <- tibble(gene = rownames(sig_sub2),
                 cell_type = apply(sig_sub2, 1, function(x) colnames(sig_sub2)[which(x == max(x))[1]]))
max_val2 <- apply(sig_sub2, 1, max)
sec_val2 <- apply(sig_sub2, 1, function(x) sort(x, decreasing = TRUE)[2])
rat_df2 <- rat_df2 %>% 
  mutate(ratio = max_val2/sec_val2)
marker_genes2 <- rat_df2 %>%
  group_by(cell_type) %>%
  arrange(desc(ratio)) %>%
  filter(row_number() <= 200) %>%
  # filter(ratio >= 7) %>% 
  ungroup %>%
  dplyr::select(gene) %>%
  unlist

# rat_df2 %>%
#   group_by(cell_type) %>%
#   arrange(desc(ratio)) %>%
#   filter(ratio >= 7) %>%
#   group_by(cell_type) %>%
#   summarise(n = n())

# marker_genes3 <- rat_df2 %>% 
#   group_by(cell_type) %>% 
#   arrange(desc(ratio)) %>% 
#   filter(ratio >= 4) %>% 
#   ungroup %>% 
#   dplyr::select(gene) %>% 
#   unlist

bulk_marker2 <- bulk_sub2[marker_genes2, ]
sig_marker2 <- sig_sub2[marker_genes2, ]
cor(sig_marker2)
# sig_marker2 <- sig_marker[, !grepl("DC", colnames(sig_marker))]

## run dSVA deconv
q_hat3 <- estimate_n_comp(Y = bulk_marker2, Theta = sig_marker2, intercept = TRUE, method = "be", B = 49, seed = 100) 
# q_hat2 <- estimate_n_comp(Y = bulk_marker2, Theta = sig_marker2, intercept = TRUE, method = "cutoff")

# q_hat2 <- estimate_n_comp(Y = bulk_marker, Theta = sig_marker, intercept = FALSE, method = "be", B = 49, seed = 100)

P_nnls3 <- dsva_for_sim(Y = bulk_marker2, 
                        Theta = sig_marker2, 
                        n_comp = q_hat3, 
                        intercept = TRUE, 
                        alg = "nnls", 
                        solver = "lsei") # using pnnls the results get interesting

rownames(P_nnls3) <- colnames(sig_marker2)
colnames(P_nnls3) <- colnames(bulk_marker2)

P_df2 <- as_tibble(t(P_nnls3), rownames = "subject_id")
P_long2 <- P_df2 %>% 
  pivot_longer(-subject_id, 
               names_to = "Cell_type",
               values_to = "Proportion") %>% 
  left_join(bulk_cols, by = "subject_id")
P_long2 %>% 
  ggplot(aes(x = Cell_type, y = Proportion, col = disease_state)) +
  geom_boxplot() +
  ggpubr::theme_pubr() +
  labs(title = "m = 1600 (COVID-19 signature from Chuwen), Remove 2 Samples \n Using Middle Cells as References", x = "Cell type")


P_long2$severity <- factor(P_long2$severity, levels = c("Healthy",
                                                        "Convalescent",
                                                        "Moderate",
                                                        "Severe",
                                                        "ICU"))
P_long2 %>% 
  filter(Cell_type == "T") %>% 
  ggplot(aes(x = severity, y = Proportion, col = severity)) +
  geom_boxplot(alpha = 0.4) +
  geom_jitter() +
  labs(title = "T-cell Proportion, m = 1600 (COVID-19 signature from Chuwen),  \n Remove 2 Samples \n Using Middle Cells as References", x = "Disease Severity") +
  scale_color_viridis_d() +
  ggpubr::theme_pubr() 

P_long2 %>% 
  filter(Cell_type == "Mono") %>% 
  ggplot(aes(x = severity, y = Proportion, col = severity)) +
  geom_boxplot(alpha = 0.4) +
  geom_jitter() +
  labs(title = "Monocyte Proportion, m = 1600 (COVID-19 signature from Chuwen),  \n Remove 2 Samples \n Using Middle Cells as References", x = "Disease Severity") +
  scale_color_viridis_d() +
  ggpubr::theme_pubr() 

# P_long2 %>% 
#   filter(Cell_type == "T") %>% 
#   ggplot(aes(x = post_infection, y = Proportion, col = severity)) +
#   geom_jitter() +
#   labs(title = "Selected ratio >= 4", x = "Disease Severity") +
#   scale_color_viridis_d() +
#   ggpubr::theme_pubr() 

## draw the comparison plots
P_null <- dsva_for_sim(Y = bulk_marker2, 
                       Theta = sig_marker2, 
                       n_comp = 0, 
                       intercept = TRUE, 
                       alg = "nnls", 
                       solver = "lsei") # using pnnls the results get interesting

P_adjust <- dsva_for_sim(Y = bulk_marker2, 
                         Theta = sig_marker2, 
                         n_comp = q_hat3, 
                         intercept = TRUE, 
                         alg = "nnls", 
                         solver = "lsei")

rownames(P_null) <- rownames(P_adjust) <- colnames(sig_marker2)
colnames(P_null) <- colnames(P_adjust) <- colnames(bulk_marker2)

P_null_long <- get_long_form(P_null) %>% 
  rename(Proportion_unadjusted = Proportion)
P_adjust_long <- get_long_form(P_adjust) %>% 
  rename(Proportion_adjusted = Proportion)

P_all <- tibble(subject_id = P_null_long$subject_id,
                prop_unadj = P_null_long$Proportion_unadjusted,
                prop_adj = P_adjust_long$Proportion_adjusted,
                severity = P_null_long$severity,
                disease_state = P_null_long$disease_state,
                cell_type = P_null_long$Cell_type)

ggplot(P_all, aes(x = prop_unadj, y = prop_adj, col = severity)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_fixed(ratio = 1) +
  labs(x = "Unadjusted proportion (q = 0)", y = paste0("Adjusted proportion (q = ", q_hat3, ")")) +
  facet_wrap(~cell_type) +
  ggpubr::theme_pubr()


# P_long %>% 
#   filter(Cell_type == "B") %>% 
#   ggplot(aes(x = severity, y = Proportion, col = severity)) +
#   geom_boxplot(alpha = 0.4) +
#   geom_jitter() +
#   labs(title = "B-Cell Proportion (m = 500)", x = "Disease Severity") +
#   scale_color_viridis_d() +
#   ggpubr::theme_pubr()