library(tidyverse)
library(ggpubr)
source("a_dSVA_functions.R")
source("r_real_data_funs.R")
# source("f_customized_package_functions.R")

# bulk_full <- readRDS("../dSVA_datasets/covid19_bulk_zhang.rds")
bulk_full <- readRDS("../dSVA_datasets/covid19_bulk_zhang_ComBat.rds")
sig_full <- readRDS("../dSVA_datasets/GSE64655_sorted_count_0d_processed.rds")
bulk_cols <- readRDS("../dSVA_datasets/covid19_col_zhang.rds")

## add the pseudo-block membership to each subject
bulk_cols$pseudo_block <- c(
  rep("1", 24),
  rep("2", 33),
  rep("1", 13)
)
outlier_subjects <- c("Mi_2", "Mo_13", "RP_4", "RP_7", "RP_11", "Re_4", "As_9")
bulk_cols$outlier <- "No"
bulk_cols$outlier[bulk_cols$subject_id %in% outlier_subjects] <- "Yes"
bulk_cols$designation <- paste(bulk_cols$pseudo_block, bulk_cols$outlier, sep = "_")

# sig_full <- sig_full[, 1:4] # remove DC cells
sig_full2 <- sig_full[rowSums(sig_full) > 10, ]

## find overlap genes
overlap_genes <- intersect(rownames(bulk_full), rownames(sig_full2)) # 15,466 genes
bulk_sub1 <- bulk_full[overlap_genes, ]
sig_sub1 <- sig_full2[overlap_genes, ]

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
  filter(row_number() <= 500) %>% 
  ungroup %>% 
  dplyr::select(gene) %>% 
  unlist

bulk_marker <- bulk_sub1[marker_genes, ]
sig_marker <- sig_sub1[marker_genes, ]
cor(sig_marker)
# sig_marker2 <- sig_marker[, !grepl("DC", colnames(sig_marker))]

## run dSVA deconv
# q_hat2 <- estimate_n_comp(Y = bulk_marker, Theta = sig_marker, intercept = TRUE, method = "be", B = 49, seed = 100)
# P_nnls2 <- dsva_for_sim(Y = bulk_marker, 
#                         Theta = sig_marker, 
#                         n_comp = q_hat2, 
#                         intercept = TRUE, 
#                         alg = "nnls", 
#                         solver = "lsei") # using pnnls the results get interesting
# rownames(P_nnls2) <- colnames(sig_marker)
# colnames(P_nnls2) <- colnames(bulk_marker)

P_nnls2 <- run_dSVA(bulk_mat = bulk_marker, sig_mat = sig_marker)

# P_df <- as_tibble(t(P_nnls2), rownames = "subject_id")
# P_long <- P_df %>% 
#   pivot_longer(-subject_id, 
#                names_to = "Cell_type",
#                values_to = "Proportion") %>% 
#   left_join(bulk_cols, by = "subject_id")
P_long <- get_long_form(P_nnls2)

## estimated cell type proportions by disease state group
P_long %>% 
  ggplot(aes(x = Cell_type, y = Proportion, col = immune_response)) +
  geom_boxplot() +
  labs(title = "Estimated cell type proportions", x = "Cell type",
       y = "Estimated proportion",
       caption = paste0("m = ", nrow(sig_marker), "; ", "q_hat = ", q_hat2, ".")) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Disease state")

## estimated cell type proportions by disease state group
P_long %>% 
  ggplot(aes(x = Cell_type, y = Proportion, col = immune_response)) +
  geom_boxplot() +
  geom_jitter(aes(shape = designation)) +
  labs(title = "Estimated cell type proportions", x = "Cell type",
       y = "Estimated proportion",
       caption = paste0("m = ", nrow(sig_marker), "; ", "q_hat = ", q_hat2, ".")) +
  theme_minimal() +
  ggsci::scale_color_tron(name = "Disease state")

## compare healthy and covid-19 patients before and after the correction
P_null <- dsva_for_sim(Y = bulk_marker, 
                       Theta = sig_marker, 
                       n_comp = 0, 
                       intercept = TRUE, 
                       alg = "nnls", 
                       solver = "lsei") # using pnnls the results get interesting

# q_hat2 <- estimate_n_comp(Y = bulk_marker, Theta = sig_marker, intercept = TRUE, method = "be", B = 49, seed = 100)
# P_adjust <- dsva_for_sim(Y = bulk_marker, 
#                          Theta = sig_marker, 
#                          n_comp = q_hat2, 
#                          intercept = TRUE, 
#                          alg = "nnls", 
#                          solver = "lsei")

rownames(P_null) <- colnames(sig_marker)
colnames(P_null) <- colnames(bulk_marker)

P_null_long <- get_long_form(P_null) %>% 
  rename(Proportion_unadjusted = Proportion)
P_adjust_long <- get_long_form(P_nnls2) %>% 
  rename(Proportion_adjusted = Proportion)

P_all <- tibble(subject_id = P_null_long$subject_id,
                prop_unadj = P_null_long$Proportion_unadjusted,
                prop_adj = P_adjust_long$Proportion_adjusted,
                immune_response = P_null_long$immune_response,
                disease_state = P_null_long$disease_state,
                cell_type = P_null_long$Cell_type,
                pseudo_block = P_null_long$pseudo_block,
                designation = P_null_long$designation)

## see which proportions are adjusted more
## by disease state
ggplot(P_all, aes(x = prop_unadj, y = prop_adj, col = immune_response, shape = designation)) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_fixed(ratio = 1) +
  labs(x = "Unadjusted proportion (q = 0)", y = paste0("Adjusted proportion (q = ", q_hat2, ")"),
       caption = paste0("m = ", nrow(sig_marker), "; ", "q_hat = ", q_hat2, ".")) +
  facet_wrap(~cell_type) +
  ggsci::scale_color_tron(name = "Disease state") +
  scale_shape_discrete(name = "Designation") +
  theme_minimal()

## by pseudo-block
ggplot(P_all, aes(x = prop_unadj, y = prop_adj, col = pseudo_block, shape = designation)) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_fixed(ratio = 1) +
  labs(x = "Unadjusted proportion (q = 0)", y = paste0("Adjusted proportion (q = ", q_hat2, ")"),
       caption = paste0("m = ", nrow(sig_marker), "; ", "q_hat = ", q_hat2, ".")) +
  facet_wrap(~cell_type) +
  ggsci::scale_color_tron(name = "Pseudo-Block") +
  scale_shape_discrete(name = "Designation") +
  theme_minimal()

P_lymph_null <- rbind(colSums(P_null[1:3, ]), P_null[-(1:3), ])
P_lymph_adjust <- rbind(colSums(P_adjust[1:3, ]), P_adjust[-(1:3), ])
rownames(P_lymph_null) <- rownames(P_lymph_adjust) <- c("Lymphocyte", colnames(sig_marker)[-(1:3)])

P_null_long2 <- get_long_form(P_lymph_null) %>% 
  rename(Proportion_unadjusted = Proportion)
P_adjust_long2 <- get_long_form(P_lymph_adjust) %>% 
  rename(Proportion_adjusted = Proportion)

P_all2 <- tibble(subject_id = P_null_long2$subject_id,
                prop_unadj = P_null_long2$Proportion_unadjusted,
                prop_adj = P_adjust_long2$Proportion_adjusted,
                immune_response = P_null_long2$immune_response,
                disease_state = P_null_long2$disease_state,
                cell_type = P_null_long2$Cell_type)

ggplot(P_all2, aes(x = prop_unadj, y = prop_adj, col = immune_response)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_fixed(ratio = 1) +
  labs(x = "Unadjusted proportion (q = 0)", y = paste0("Adjusted proportion (q = ", q_hat2, ")")) +
  facet_wrap(~cell_type) +
  ggpubr::theme_pubr()

## lymphocyte vs monocyte
P_df <- as_tibble(t(P_lymph_adjust), rownames = "subject_id")
P_long <- P_df %>% 
  pivot_longer(-subject_id, 
               names_to = "Cell_type",
               values_to = "Proportion") %>% 
  left_join(bulk_cols, by = "subject_id")
P_long %>% 
  ggplot(aes(x = Cell_type, y = Proportion, col = immune_response)) +
  geom_boxplot() +
  labs(title = "Estimated cell type proportions", x = "Cell type",
       y = "Estimated proportion",
       caption = paste0("m = ", nrow(sig_marker), "; ", "q_hat = ", q_hat2, ".")) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Disease state")


