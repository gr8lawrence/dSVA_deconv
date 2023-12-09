## use combat-seq to adjust the counts
library(sva)
bulk_full1 <- readRDS("../dSVA_datasets/covid_bulk_count_processed.rds")
sig_full <- readRDS("../dSVA_datasets/GSE64655_sorted_count_0d_processed.rds")
bulk_cols <- readRDS("../dSVA_datasets/covid_col.rds")

## remove the two suspect samples 
bulk_full <- bulk_full1[, !grepl("S155", colnames(bulk_full1)) & !grepl("S179", colnames(bulk_full1))]
bulk_full2 <- bulk_full[rowSums(bulk_full) > 0, ]

bulk_cols2 <- bulk_cols %>% 
  filter(!grepl("S155", subject_id) & !grepl("S179", subject_id))
bulk_cols2 %>% 
  group_by(disease_state) %>% 
  summarise(n = n())
# all(colnames(bulk_full) == bulk_cols2$subject_id) # TRUE
batch <- c(rep(1, 15), rep(2, 17)) # group convalescent to COVID-19
bulk_adjust <- ComBat_seq(counts = bulk_full, batch = batch)
bulk_adjust <- bulk_adjust[rowSums(bulk_adjust) > 0, ]


## the PCs for the original bulk after removing the two samples
pca <- prcomp(t(bulk_full2), scale. = TRUE, center = TRUE)
samp_id <- str_split_i(colnames(bulk_full2), "_", 1)
plot(pca$sdev[1:20]^2, main = "Elbow Plot \nEigenvalues of Raw Bulk Matrix (COVID-19 Study)", ylab = "Eigenvalue",
     log = "y")
autoplot(pca, data = bulk_cols2, colour = "disease_state") +
  labs(title = "PC of COVID-19 Bulk Samples", subtitle = "By Disease States after Removing 2 samples") +
  geom_label_repel(aes(label = samp_id)) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Disease state")

autoplot(pca, data = bulk_cols2, colour = "severity") +
  labs(title = "PC of COVID-19 Bulk Samples", subtitle = "By Disease Severity after Removing 2 samples") +
  geom_label_repel(aes(label = samp_id)) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Disease severity")


## Let's check the altered bulk data
pca <- prcomp(t(bulk_adjust), scale. = TRUE, center = TRUE)
samp_id <- str_split_i(colnames(bulk_adjust), "_", 1)
plot(pca$sdev[1:20]^2, main = "Elbow Plot \nEigenvalues of Raw Bulk Matrix (COVID-19 Study)", ylab = "Eigenvalue",
     log = "y")
autoplot(pca, data = bulk_cols2, colour = "disease_state") +
  labs(title = "PC of COVID-19 Bulk Samples", subtitle = "By Disease States after Removing 2 samples and ComBat-Seq") +
  geom_label_repel(aes(label = samp_id)) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Disease state")

autoplot(pca, data = bulk_cols2, colour = "severity") +
  labs(title = "PC of COVID-19 Bulk Samples", subtitle = "By Disease Severity after Removing 2 samples and ComBat-Seq") +
  geom_label_repel(aes(label = samp_id)) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Disease severity")

## find overlap genes
sig_full <- sig_full[, 1:4] # remove DC cells
overlap_genes <- intersect(rownames(bulk_adjust), rownames(sig_full2)) # 15,466 genes
bulk_sub1 <- bulk_adjust[overlap_genes, ]
sig_sub1 <- sig_full2[overlap_genes, ]


bulk_cols$severity <- factor(bulk_cols$severity, levels = c("Healthy",
                                                            "Convalescent",
                                                            "Moderate",
                                                            "Severe",
                                                            "ICU"))
sig_full2 <- sig_full[rowSums(sig_full) > 10, ]

## find overlap genes
overlap_genes <- intersect(rownames(bulk_adjust), rownames(sig_full2)) # 15,466 genes
bulk_sub1 <- bulk_adjust[overlap_genes, ]
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
  filter(row_number() <= 200) %>% 
  ungroup %>% 
  dplyr::select(gene) %>% 
  unlist

bulk_marker <- bulk_sub1[marker_genes, ]
sig_marker <- sig_sub1[marker_genes, ]
cor(sig_marker)

P_nnls2 <- dsva_for_sim(Y = bulk_marker, 
                        Theta = sig_marker, 
                        n_comp = 0, 
                        intercept = TRUE, 
                        alg = "pnnls", 
                        solver = "lsei")
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
  labs(title = "m = 800 (Day 0 Reference, no DC, \n2 Sample Removed, ComBat)", x = "Cell type") +
  ggpubr::theme_pubr()

saveRDS(P_nnls2, "../dSVA_datasets/ComBat_P_hat.rds")
P_null <- readRDS("../dSVA_datasets/dSVA_Pnull_q_0.rds")
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
P_adjust_long <- get_long_form(P_nnls2) %>% 
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
  labs(x = "Unadjusted proportion (q = 0)", y = "Adjusted proportion (ComBat)") +
  facet_wrap(~cell_type) +
  ggpubr::theme_pubr()


