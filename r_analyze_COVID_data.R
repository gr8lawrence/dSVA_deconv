library(tidyverse)
library(ggpubr)
library(ggfortify)
source("a_dSVA_functions.R")
source("r_real_data_funs.R")
# source("f_customized_package_functions.R")

# Chuwen Data -------------------------------------------------------------

# controls_df <- read_delim("../dSVA_datasets/chuwen_bulk_controls.txt")
# cases_df <- read_delim("../dSVA_datasets/chuwen_bulk_cases.txt")
# sig_df <- read_delim("../dSVA_datasets/chuwen_sig_controls.txt")
# 
# ## process data
# cons <- data.matrix(controls_df[, -1])
# cases <- data.matrix(cases_df[, -1])
# both <- cbind(cons, cases)
# rownames(both) <- controls_df$genesymbol
# covid_cols <- readRDS("../dSVA_datasets/covid_col.rds")
# 
# sig <- data.matrix(sig_df[, -1])
# rownames(sig) <- sig_df$NAME
# genes <- intersect(rownames(both), rownames(sig))
# 
# Y <- both[genes, ]
# Theta <- sig[genes, ]
# 
# q_hat <- estimate_n_comp(Y = Y, Theta = Theta, intercept = TRUE, method = "be", B = 49, seed = 100) # 5
# P_nnls <- dsva_for_sim(Y = Y, 
#                        Theta = Theta, 
#                        n_comp = q_hat, 
#                        intercept = TRUE, 
#                        alg = "nnls", 
#                        solver = "lsei")
# 
# rownames(P_nnls) <- colnames(Theta)
# colnames(P_nnls) <- colnames(Y)
# 
# P_nnls <- apply(P_nnls, 2, function(x) x/sum(x))
# 
# P_nnls2 <- dsva_for_sim(Y = Y, 
#                        Theta = Theta, 
#                        n_comp = 1, 
#                        intercept = TRUE, 
#                        alg = "nnls", 
#                        solver = "lsei")
# 
# rownames(P_nnls2) <- colnames(Theta)
# colnames(P_nnls2) <- colnames(Y)
# 
# P_nnls2 <- apply(P_nnls2, 2, function(x) x/sum(x))


# ABIS-seq Signature Matrix -----------------------------------------------

## Now, let's analyze it using a different reference (GSE 64655 or ABIS-seq)
bulk_full1 <- readRDS("../dSVA_datasets/covid_bulk_count_processed.rds")
# sig_full <- readRDS("../dSVA_datasets/GSE64655_sorted_count_processed.rds")
# sig_full <- readRDS("../dSVA_datasets/GSE64655_sorted_count_0d_processed.rds")
sig_full <- read.table("../dSVA_datasets/sigmatrixRNAseq.txt", header = TRUE, sep = "\t") %>% 
  data.matrix()

## remove the two suspect samples 
bulk_full <- bulk_full1[, !grepl("S155", colnames(bulk_full1)) & !grepl("S179", colnames(bulk_full1))]

# sig_full <- sig_full[, 1:4] # remove DC cells

bulk_cols <- readRDS("../dSVA_datasets/covid_col.rds")
bulk_cols$severity <- factor(bulk_cols$severity, levels = c("ICU",
                                                            "Severe",
                                                            "Moderate",
                                                            "Convalescent",
                                                            "Healthy"))
bulk_cols$disease_state_bin <- case_when(
  bulk_cols$disease_state == "Healthy" | bulk_cols$disease_state == "Convalescent" ~ "Healthy/Convalescent",
  .default = "COVID-19"
)
bulk_cols$post_infection3 <- replace_na(bulk_cols$post_infection, 0)

sig_full2 <- sig_full[rowSums(sig_full) > 10, ]

## find overlap genes
overlap_genes <- intersect(rownames(bulk_full), rownames(sig_full2)) # 15,466 genes
bulk_sub1 <- bulk_full[overlap_genes, ]
sig_sub1 <- sig_full2[overlap_genes, ]

saveRDS(overlap_genes, "../dSVA_datasets/covid_deconv_genes.rds")

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
# rat_df <- tibble(gene = rownames(sig_sub1),
#                  cell_type = apply(sig_sub1, 1, function(x) colnames(sig_sub1)[which(x == max(x))[1]]))
# max_val <- apply(sig_sub1, 1, max)
# sec_val <- apply(sig_sub1, 1, function(x) sort(x, decreasing = TRUE)[2])
# rat_df <- rat_df %>% 
#   mutate(ratio = max_val/sec_val)
# marker_genes <- rat_df %>% 
#   group_by(cell_type) %>% 
#   arrange(desc(ratio)) %>% 
#   filter(row_number() <= 300) %>% 
#   ungroup %>% 
#   dplyr::select(gene) %>% 
#   unlist

bulk_marker <- bulk_sub1[marker_genes, ]
sig_marker <- sig_sub1[marker_genes, ]
cor(sig_marker)
# sig_marker2 <- sig_marker[, !grepl("DC", colnames(sig_marker))]


# Check residuals ---------------------------------------------------------

  # Theta = Theta_deconv
  # intercept = TRUE
  # bulk_col2 = bulk_col[bulk_col$sample %in% colnames(Y_sub),]
  # Y = Y_sub[, bulk_col2$sample]

## check what the PC1 and PC2 of R is
Y <- bulk_marker
Theta <- sig_marker
intercept <- TRUE
bulk_col2 <- bulk_cols[!grepl("S155", bulk_cols$subject_id) & !grepl("S179", bulk_cols$subject_id), ]
bulk_col2$sample <- str_split_i(bulk_col2$subject_id, "_", 1)
  
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
plot(pca2$sdev^2/sum(pca2$sdev^2), main = "Bulk Residuals PCA Variation Explained",
     sub = "Using all cell types in a signature matrix",
     ylab = "Variance explained by PC", xlab = "PC") %>% print()

p2 <- autoplot(pca2, data = bulk_col2, colour = "disease_state") +
  labs(title = "PCA Plot on COVID-19 Bulk Residuals") +
  # geom_label_repel(aes(label = sample)) +
  ggthemes::theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Disease state", values = c("#ffd380", "#8a508f", "#2c4875")) 
p2
ggsave("plots/COVID_19/Residual_PC_Plot_Full_Theta.pdf")

# See if the residuals are correlated with Monocyte.c signatures ----------

mono_c_sig <- sig_marker[, 1]
cor(R, mono_c_sig)
# cor(t(R), mono_c_sig)

## What if we remove the monocyte C?
sig_marker2 <- sig_marker[, -1]
Theta <- sig_marker2
intercept <- TRUE

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
pca2b <- prcomp(t(R), scale. = TRUE, center = TRUE)
# pca$sdev
plot(pca2b$sdev^2/sum(pca2b$sdev^2), main = "Bulk Residuals PCA Variation Explained",
     sub = "Without classical monocytes in the signature matrix",
     ylab = "Variance explained by PC", xlab = "PC") %>% print()

p3 <- autoplot(pca2b, data = bulk_col2, colour = "disease_state") +
  labs(title = "PCA Plot on COVID-19 Bulk Residuals", subtitle = "Without classical monocytes in the signature matrix") +
  # geom_label_repel(aes(label = sample)) +
  ggthemes::theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Disease state", values = c("#ffd380", "#8a508f", "#2c4875")) 
p3
ggsave("plots/COVID_19/Residual_PC_Plot_Theta_Minus_MC.pdf")


# get_R_pca(bulk_marker, sig_marker, bulk_col2 = bulk_col2)

## check if PC1 loadings correlate with the proportions of monocytes

pc1_val <- pca2b$x[, "PC1"]
nnls_mono_c <- P_nnls["Monocytes.C", ]
# plot(pc1_val, nnls_mono_c)
disease_states <- bulk_col2$disease_state
names(disease_states) <- bulk_col2$subject_id
colors <- ifelse(disease_states == "Healthy", "#ffd380","#8a508f")
pdf("plots/COVID_19/NNLS_Mono_C_vs_PC1.pdf")
plot(pc1_val, nnls_mono_c, col = colors, 
     # main = "NNLS C. Monocyte Props vs. PC 1 Coordinates",
     sub = "Without classical monocytes in the signature matrix",
     xlab = "PC1 of R", ylab = "NNLS classical monocyte proportions")
dev.off()

# Run dSVA ----------------------------------------------------------------

## run dSVA deconv
q_hat2 <- estimate_n_comp(Y = bulk_marker, Theta = data.matrix(sig_marker), intercept = TRUE, method = "cutoff", B = 49, seed = 100)
# q_hat2 <- 0
# q_hat2 <- estimate_n_comp(Y = bulk_marker, Theta = sig_marker, intercept = TRUE, method = "cutoff")
# q_hat2 <- estimate_n_comp(Y = bulk_marker, Theta = sig_marker, intercept = FALSE, method = "be", B = 49, seed = 100)
P_ls <- run_dSVA(bulk_marker, sig_marker, intercept = TRUE, exclude = NULL, alg = "nnls", solver = "lsei", method = "none")
P_ls2 <- run_dSVA(bulk_marker, sig_marker, intercept = TRUE, exclude = NULL, alg = "nnls", solver = "lsei", method = "cutoff")

## see different version with exclude = 1
# P_ls3 <- run_dSVA(bulk_marker, sig_marker, intercept = TRUE, exclude = 1, alg = "nnls", solver = "lsei", method = "cutoff")
P_ls3 <- run_dSVA(bulk_marker, sig_marker, intercept = TRUE, exclude = 1, alg = "nnls", solver = "lsei", method = "be")

P_ls4 <- run_dSVA(bulk_marker, sig_marker, intercept = TRUE, exclude = 1, alg = "nnls", solver = "lsei", method = "none")

P_ls5 <- run_dSVA(bulk_marker, sig_marker, intercept = TRUE, exclude = 1, alg = "nnls", solver = "lsei", method = "cutoff")

P_ls6 <- run_dSVA(bulk_marker, sig_marker, intercept = TRUE, exclude = 1, alg = "nnls", solver = "lsei", q = 2)



## SSR comparison between the corrected model and the uncorrected model
P_ls$SSR/P_ls3$SSR

P_ls$Gamma_hat

cor(sig_marker , P_ls2$Gamma_hat)
cor(sig_marker , P_ls3$Gamma_hat)
cor(sig_marker , P_ls5$Gamma_hat)


# Plot the results --------------------------------------------------------
P_nnls <- P_ls4$P_hat
P_nnls2 <- P_ls3$P_hat
P_nnls3 <- P_ls5$P_hat
P_nnls4 <- P_ls6$P_hat

cor(sig_marker, P_ls3$Gamma_hat)

P_long_nnls <- get_long_form(P_nnls) %>% 
  mutate(adjust = "Unadjusted")
P_long_nnls2 <- get_long_form(P_nnls2) %>% 
  mutate(adjust = "Adjusted (dSVA)")
P_long_nnls3 <- get_long_form(P_nnls3) %>% 
  mutate(adjust = "Adjusted (dSVA)")
P_long_nnls4 <- get_long_form(P_nnls4) %>% 
  mutate(adjust = "Adjusted (dSVA)")
P_long <- rbind(P_long_nnls, P_long_nnls2)
P_long$disease_state_bin <- factor(P_long$disease_state_bin, levels = c("Healthy/Convalescent", "COVID-19"))


P_long_nnls %>% 
  ggplot(aes(x = Cell_type, y = Proportion, fill = disease_state)) +
  geom_boxplot() +
  labs(title = "Deconvoluted PBMC proportions", 
       subtitle = "Unadjusted",
       x = "Cell type",
       caption = paste("Using ABIS-seq as reference; q =", P_ls$q_hat)) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f", "#2c4875")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background = element_blank(),
        legend.position = c(0.01, 0.98), 
        legend.justification = c(0.01, 0.98),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5)) +
  ggsci::scale_color_tron(name = "Disease state")

P_long_nnls2 %>% 
  ggplot(aes(x = Cell_type, y = Proportion, fill = disease_state)) +
  geom_boxplot() +
  labs(title = "Deconvoluted PBMC proportions", 
       subtitle = "dSVA-Adjusted (BE/TW Method)",
       x = "Cell type",
       caption = paste("Using ABIS-seq as reference; q =", P_ls3$q_hat)) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f", "#2c4875")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background = element_blank(),
        legend.position = c(0.01, 0.98), 
        legend.justification = c(0.01, 0.98),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5)) +
  ggsci::scale_color_tron(name = "Disease state")

P_long_nnls3 %>% 
  ggplot(aes(x = Cell_type, y = Proportion, fill = disease_state)) +
  geom_boxplot() +
  labs(title = "Deconvoluted PBMC proportions", 
       subtitle = "dSVA-Adjusted (Cutoff Method)",
       x = "Cell type",
       caption = paste("Using ABIS-seq as reference; q =", P_ls5$q_hat)) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f", "#2c4875")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background = element_blank(),
        legend.position = c(0.01, 0.98), 
        legend.justification = c(0.01, 0.98),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5)) +
  ggsci::scale_color_tron(name = "Disease state")

P_long_nnls4 %>% 
  ggplot(aes(x = Cell_type, y = Proportion, fill = disease_state)) +
  geom_boxplot() +
  labs(title = "Deconvoluted PBMC proportions", 
       subtitle = "dSVA-Adjusted (q = 2)",
       x = "Cell type",
       caption = paste("Using ABIS-seq as reference; q =", 2)) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f", "#2c4875")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background = element_blank(),
        legend.position = c(0.01, 0.98), 
        legend.justification = c(0.01, 0.98),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5)) +
  ggsci::scale_color_tron(name = "Disease state")

# P_df <- as_tibble(t(P_nnls2), rownames = "subject_id")
# P_long <- P_df %>% 
#   pivot_longer(-subject_id, 
#                names_to = "Cell_type",
#                values_to = "Proportion") %>% 
#   left_join(bulk_cols, by = "subject_id")
# P_long %>%
#   ggplot(aes(x = Cell_type, y = Proportion, col = disease_state)) +
#   geom_boxplot() +
#   labs(title = "m = 1600 (Day 0 Reference, no DC, \n2 Sample Removed)", x = "Cell type") +
#   ggpubr::theme_pubr() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# P_long %>%
#   ggplot(aes(x = Cell_type, y = Proportion, fill = disease_state_bin)) +
#   geom_boxplot() +
#   facet_grid(adjust ~ .) +
#   labs(title = "Deconvoluted PBMC proportions", x = "Cell type",
#        caption = paste("Using ABIS-seq as reference; q =", q_hat2)) +
#   theme_linedraw() +
#   scale_fill_manual(name = "Disease State", values = c("#D5896F", "#70A288")) +
#   # ggsci::scale_color_tron(name = "Disease State") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#         legend.position = c(0.01, 0.98), 
#         legend.justification = c(0.01, 0.98),
#         legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
#         panel.grid = element_line(color = "grey"))
        # legend.position = "bottom",
        # legend.box = "horizontal") 

col_vals <- c("#ffd380", "#ff8531", "#ff6361", "#bc5090", "#8a508f", "#2c4875")

P_long %>%
  ggplot(aes(x = Cell_type, y = Proportion, fill = disease_state_bin)) +
  geom_boxplot() +
  facet_grid(adjust ~ .) +
  labs(title = "Deconvoluted PBMC proportions", x = "Cell type",
       caption = paste("Using ABIS-seq as reference; q =", P_ls3$q_hat)) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background = element_blank(),
        legend.position = c(0.01, 0.98), 
        legend.justification = c(0.01, 0.98),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5)) +
  ggsci::scale_color_tron(name = "Disease state")

ggsave("plots/COVID_19/COVID_19_main_prop_figure.pdf")

# p_pDC <- P_long %>%
#   filter(Cell_type == "pDCs") %>% 
#   ggplot(aes(x = disease_state_bin, y = Proportion, fill = disease_state_bin)) +
#   geom_boxplot() +
#   facet_grid(.~adjust) +
#   labs(title = "pDC Abundance", x = "Cell type"
#        # ,caption = paste("Using ABIS-seq as reference;", "q =", q_hat2)
#        ) +
#   theme_linedraw() +
#   theme(axis.text.x = element_blank(),
#         panel.grid = element_line(color = "grey")) +
#   scale_fill_manual(name = "Disease State", values = c("#D5896F", "#70A288")) 
  # ggsci::scale_color_tron(name = "Disease State") 
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p_pDC1 <- P_long %>%
  filter(Cell_type == "pDCs") %>% 
  ggplot(aes(x = disease_state_bin, y = Proportion, fill = disease_state_bin)) +
  geom_boxplot() +
  facet_grid(.~adjust) +
  labs(title = "pDC Abundance", x = "Disease state"
       # ,caption = paste("Using ABIS-seq as reference;", "q =", q_hat2)
  ) +
  theme_base(base_size = 12, base_family = "Times") +
  theme(axis.text.x = element_blank(),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5),
        legend.position = "bottom",
        plot.background=element_blank()) +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f"))  
p_pDC <- p_pDC1 + ggpubr::stat_compare_means(method = "wilcox.test")
# Change method
# p + stat_compare_means(method = "t.test")

# p_plasma <- P_long %>%
#   filter(Cell_type == "Plasmablasts") %>% 
#   ggplot(aes(x = disease_state_bin, y = Proportion, fill = disease_state_bin)) +
#   geom_boxplot() +
#   facet_grid(.~adjust) +
#   labs(title = "Plasmablast Abundance", x = "Cell type"
#        # ,caption = paste("Using ABIS-seq as reference; q_hat =", q_hat2)
#        ) +
#   theme_linedraw() +
#   theme(axis.text.x = element_blank(),
#         panel.grid = element_line(color = "grey")) +
#   scale_fill_manual(name = "Disease State", values = c("#D5896F", "#70A288")) 
  # ggsci::scale_color_tron(name = "Disease State") 

p_plasma1 <- P_long %>%
  filter(Cell_type == "Plasmablasts") %>% 
  ggplot(aes(x = disease_state_bin, y = Proportion, fill = disease_state_bin)) +
  geom_boxplot() +
  facet_grid(.~adjust) +
  labs(title = "Plasmablast Abundance", x = "Disease state"
       # ,caption = paste("Using ABIS-seq as reference; q_hat =", q_hat2)
  ) +
  theme_base(base_size = 12, base_family = "Times") +
  theme(axis.text.x = element_blank(),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5),
        legend.position = "bottom",
        plot.background=element_blank()) +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f")) 
p_plasma <- p_plasma1 + ggpubr::stat_compare_means(method = "wilcox.test")


# p_CD8 <- P_long %>%
#   filter(grepl("CD8", Cell_type)) %>%
#   group_by(subject_id) %>%
#   summarise(Proportion = sum(Proportion)) %>%
#   right_join(P_long %>% select(-Proportion) %>% distinct()) %>%
#   ggplot(aes(x = disease_state_bin, y = Proportion, fill = disease_state_bin)) +
#   geom_boxplot() +
#   facet_grid(.~adjust) +
#   labs(title = "Total CD8 Cell Abundance", x = "Cell type"
#        # ,caption = paste("Using ABIS-seq as reference; q_hat =", q_hat2)
#        ) +
#   theme_linedraw() +
#   theme(axis.text.x = element_blank(),
#         panel.grid = element_line(color = "grey")) +
#   scale_fill_manual(name = "Disease State", values = c("#D5896F", "#70A288")) 
#   # ggsci::scale_color_tron(name = "Disease State") 

p_CD81 <- P_long %>%
  filter(grepl("CD8", Cell_type)) %>%
  group_by(subject_id) %>%
  summarise(Proportion = sum(Proportion)) %>%
  right_join(P_long %>% dplyr::select(-Proportion) %>% distinct()) %>%
  ggplot(aes(x = disease_state_bin, y = Proportion, fill = disease_state_bin)) +
  geom_boxplot() +
  facet_grid(.~adjust) +
  labs(title = "Total CD8 Cell Abundance", x = "Disease state"
       # ,caption = paste("Using ABIS-seq as reference; q_hat =", q_hat2)
  ) +
  theme_base(base_size = 12, base_family = "Times") +
  theme(axis.text.x = element_blank(),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5),
        legend.position = "bottom",
        plot.background=element_blank()) +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f")) 
p_CD8 <- p_CD81 + ggpubr::stat_compare_means(method = "wilcox.test")


# p_CD8_post_infection <- P_long %>%
#   filter(grepl("CD8", Cell_type)) %>%
#   group_by(subject_id) %>%
#   summarise(Proportion = sum(Proportion)) %>%
#   right_join(P_long %>% select(-Proportion) %>% distinct()) %>%
#   ggplot(aes(x = post_infection3, y = Proportion, col = disease_state_bin)) +
#   geom_point() +
#   facet_grid(adjust~.) +
#   labs(title = "Total CD8 Cell Abundance", x = "Days post infection",
#        caption = paste("Using ABIS-seq as reference; q =", q_hat2)) +
#   theme_linedraw() +
#   theme(panel.grid = element_line(color = "grey")) +
#   scale_color_manual(name = "Disease State", values = c("#D5896F", "#70A288"))
#   # ggsci::scale_color_tron(name = "Disease State") 

p_CD8_post_infection <- P_long %>%
  filter(grepl("CD8", Cell_type)) %>%
  group_by(subject_id) %>%
  summarise(Proportion = sum(Proportion)) %>%
  right_join(P_long %>% dplyr::select(-Proportion) %>% distinct()) %>%
  ggplot(aes(x = post_infection3, y = Proportion, col = disease_state_bin)) +
  geom_point() +
  facet_grid(adjust~.) +
  labs(title = "Total CD8 Cell Abundance", x = "Days post infection",
       caption = paste("Using ABIS-seq as reference; q =", q_hat2)) +
  theme_base(base_size = 12, base_family = "Times") +
  theme(axis.text.x = element_blank(),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5),
        legend.position = "bottom",
        plot.background=element_blank()) +
  scale_color_manual(name = "Disease State", values = c("#ffd380", "#8a508f"))

ggpubr::ggarrange(p_pDC, p_plasma, p_CD8, p_CD8_post_infection, common.legend = TRUE)


ggsave("plots/COVID_19/COVID_19_cell_subset_prop_figure.pdf")


# P_long %>%
#   filter(grepl("CD8", Cell_type)) %>%
#   filter(post_infection3 < 30) %>% 
#   group_by(subject_id) %>%
#   summarise(Proportion = sum(Proportion)) %>%
#   right_join(P_long %>% select(-Proportion) %>% distinct()) %>%
#   ungroup() %>% 
#   ggplot(aes(x = post_infection3, y = Proportion)) +
#   geom_point(aes(col = severity)) +
#   geom_smooth() +
#   facet_grid(adjust~.) +
#   labs(title = "Total CD8 Cell Abundance",
#        x = "Days post infection",
#        caption = paste("Using ABIS-seq as reference; q_hat =", q_hat2, "; convalescent patient filtered")) +
#   theme_light() +
#   ggsci::scale_color_tron(name = "Severity") 


P_long %>%
  filter(grepl("CD8", Cell_type)) %>%
  filter(post_infection3 < 30) %>% 
  group_by(subject_id) %>%
  summarise(Proportion = sum(Proportion)) %>%
  right_join(P_long %>% select(-Proportion) %>% distinct()) %>%
  ungroup() %>% 
  ggplot(aes(x = post_infection3, y = Proportion)) +
  geom_point(aes(col = severity)) +
  xlim(0, 25) +
  geom_smooth() +
  facet_grid(adjust~.) +
  labs(title = "Total CD8 Cell Abundance",
       x = "Days post infection",
       caption = paste("Using ABIS-seq as reference; q_hat =", q_hat2, "; convalescent patient filtered")) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_color_manual(name = "Severity", values = col_vals[5:1]) +
  theme(legend.position = "bottom",
        plot.background=element_blank())

ggsave("plots/COVID_19/COVID_19_CD8_abundance_trend.pdf")


# P_long %>% 
#   filter(Cell_type == "T") %>% 
#   ggplot(aes(x = severity, y = Proportion, col = severity)) +
#   geom_boxplot(alpha = 0.4) +
#   geom_jitter() +
#   labs(title = "T-Cell Proportion (m = 1600, Day 0 Reference, No DC, \n2 Sample Removed)", x = "Disease Severity") +
#   scale_color_viridis_d() +
#   ggpubr::theme_pubr() 
# 
# P_long %>% 
#   filter(Cell_type == "Mono") %>% 
#   ggplot(aes(x = severity, y = Proportion, col = severity)) +
#   geom_boxplot(alpha = 0.4) +
#   geom_jitter() +
#   labs(title = "Monocyte Proportion (m = 1600, Day 0 Reference, No DC, \n2 Sample Removed)", x = "Disease Severity") +
#   scale_color_viridis_d() +
#   ggpubr::theme_pubr()

## compare healthy and covid-19 patients before and after the correction
P_null <- dsva_for_sim(Y = bulk_marker, 
                       Theta = sig_marker, 
                       n_comp = 0, 
                       intercept = TRUE, 
                       alg = "nnls", 
                       solver = "lsei") # using pnnls the results get interesting

P_adjust <- dsva_for_sim(Y = bulk_marker, 
                         Theta = sig_marker, 
                         n_comp = q_hat2, 
                         intercept = TRUE, 
                         alg = "nnls", 
                         solver = "lsei")

rownames(P_null) <- rownames(P_adjust) <- colnames(sig_marker)
colnames(P_null) <- colnames(P_adjust) <- colnames(bulk_marker)

saveRDS(P_adjust, "../dSVA_datasets/dSVA_Padj_q_2_ABIS.rds")
saveRDS(P_null, "../dSVA_datasets/dSVA_Pnull_q_0_ABIS.rds")

P_null_long <- get_long_form(P_null) %>% 
  dplyr::rename(Proportion_unadjusted = Proportion)
P_adjust_long <- get_long_form(P_adjust) %>% 
  dplyr::rename(Proportion_adjusted = Proportion)

P_all <- tibble(subject_id = P_null_long$subject_id,
                prop_unadj = P_null_long$Proportion_unadjusted,
                prop_adj = P_adjust_long$Proportion_adjusted,
                severity = P_null_long$severity,
                disease_state = P_null_long$disease_state,
                disease_state_bin = P_null_long$disease_state_bin,
                cell_type = P_null_long$Cell_type)


# P_all <- P_null_long %>% 
#   left_join(P_adjust_long %>% dplyr::select(subject_id, Proportion_adjusted), by = "subject_id")

ggplot(P_all, aes(x = prop_unadj, y = prop_adj, col = disease_state_bin)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_fixed(ratio = 1) +
  labs(x = "Unadjusted proportion (q = 0)", y = paste0("Adjusted proportion (q = ", q_hat2, ")")) +
  facet_wrap(~cell_type) +
  scale_color_manual(name = "Disease State", values = c("#8a508f", "#ffd380")) +
  # ggsci::scale_color_tron(name = "Disease state") +
  # theme_linedraw() +
  theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        # panel.grid = element_line(color = "grey"),
        plot.background=element_blank()) 
ggsave("plots/COVID_19/cell_subset_comparison.pdf")

ggplot(P_all, aes(x = prop_unadj, y = prop_adj, col = severity)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_fixed(ratio = 1) +
  labs(x = "Unadjusted proportion (q = 0)", y = paste0("Adjusted proportion (q = ", q_hat2, ")")) +
  facet_wrap(~cell_type) +
  # ggsci::scale_color_frontiers(name = "Severity")+
  theme_base(base_size = 12, base_family = "Times") +
  scale_color_manual(name = "Severity", values = col_vals[5:1]) +
  # theme_linedraw() +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        plot.background=element_blank()) 
  # ggpubr::theme_pubr()

## plot the total lymphocyte amount
lymph_ct_ind <- grepl("T", rownames(P_nnls2)) | grepl("B", rownames(P_nnls2)) | grepl("NK", rownames(P_nnls2))
P_lymph <- P_nnls2[lymph_ct_ind, ]
P_other <- P_nnls2[!lymph_ct_ind, ]
P_aggreg <- rbind(colSums(P_lymph), colSums(P_other))
rownames(P_aggreg) <- c("Lymphocyte", "Other")
P_aggreg_long <- get_long_form(P_aggreg)
P_aggreg_long %>% 
  ggplot(aes(x = Cell_type, y = Proportion, col = disease_state_bin)) + 
  geom_boxplot() +
  labs( caption = paste("Using ABIS-seq as reference; q_hat =", q_hat2))

## Plot the T, B, NK, Monocytes
P_nnls2
# get_all_prop <- function(P_nnls2) {
#   T_ct_ind <- lymph_ct_ind <- grepl("T", rownames(P_nnls2))
#   B_ct_ind <- lymph_ct_ind <- grepl("B", rownames(P_nnls2))
#   NK_ct_ind <- lymph_ct_ind <- grepl("NK", rownames(P_nnls2))
#   # Mono_ct_ind <- lymph_ct_ind <- grepl("Mono", rownames(P_nnls2))
#   P_T <- P_nnls2[T_ct_ind, ]
#   P_B <- P_nnls2[B_ct_ind, ]
#   P_NK <- P_nnls2[NK_ct_ind, ]
#   # P_Mono <- P_nnls2[Mono_ct_ind, ]
#   all_prop <- rbind(colSums(P_T), colSums(P_B), P_NK, colSums(P_Mono))
#   rownames(all_prop) <- c("T cells", "B cells", "NK cells", "Monocytes")
#   return(all_prop)
# }

## function to combine deconvolved cell types into a large category
get_all_prop <- function(P_nnls2) {
  T_ct_ind <- lymph_ct_ind <- grepl("T", rownames(P_nnls2))
  B_ct_ind <- lymph_ct_ind <- grepl("B", rownames(P_nnls2))
  NK_ct_ind <- lymph_ct_ind <- grepl("NK", rownames(P_nnls2))
  # Mono_ct_ind <- lymph_ct_ind <- grepl("Mono", rownames(P_nnls2))
  P_T <- P_nnls2[T_ct_ind, ]
  P_B <- P_nnls2[B_ct_ind, ]
  P_NK <- P_nnls2[NK_ct_ind, ]
  # P_Mono <- P_nnls2[Mono_ct_ind, ]
  all_prop <- rbind(colSums(P_T), colSums(P_B), P_NK)
  rownames(all_prop) <- c("T cells", "B cells", "NK cells")
  return(all_prop)
}

## Plot the deconvolution results in Chuwen's methods against true single-cell data
## create a metadata
meta_data <- tibble(
  sample = paste0("cov", c("01", "02", "03", "04", "07", "08", "09", "10", "11", "12", "17", "18")),
  disease = c(rep("COVID-19", 4), rep("Control", 3), rep("COVID-19", 3), rep("Control", 2)),
  subject_id = c("S150", "S153", "S147", "S152", "S071", "S063", "S070", "S155", "S175", "S157", "S066", "S062")
)

## read Chuwen's data
case_true <- read_csv("../dSVA_datasets/case-celltype.csv") %>% 
  select(-1) %>% 
  dplyr::rename(cell_type = y,
         true_prop = prop) %>%
  left_join(meta_data, by = "sample")
control_true <- read_csv("../dSVA_datasets/control-celltype.csv") %>%   
  select(-1) %>% 
  dplyr::rename(cell_type = y,
         true_prop = prop) %>% 
  left_join(meta_data, by = "sample")
all_true <- rbind(case_true, control_true)


## enumerate all results here
all_prop <- get_all_prop(P_nnls2) # q = 3
all_prop2 <- get_all_prop(P_nnls3) # q = 1
all_prop3 <- get_all_prop(P_nnls4) # q = 2
# all_prop_n <- get_all_prop(P_nnls) # q = 0

# all_prop4 <- get_all_prop(P_nnls5)


# all_prop <- apply(all_prop, 2, function(x) x/sum(x))


all_est <- all_prop3 %>% 
  as_tibble(rownames = "cell_type") %>% 
  pivot_longer(-cell_type, 
               names_to = "subject_id", 
               values_to = "est_prop") %>% 
  mutate(subject_id = str_split_i(subject_id, "_", 1))

## combine the truths and the estimates
all_df <- all_true %>% 
  left_join(all_est, by = c("subject_id", "cell_type"))
all_df$cell_type <- factor(all_df$cell_type, levels = c("T cells", "B cells", "NK cells", "Monocytes"))

## calculate summary statistics
summary_df <- all_df %>% 
  group_by(disease) %>% 
  summarise(pearson = cor(true_prop, est_prop, use = "complete.obs"),
            spearman = cor(true_prop, est_prop, "spearman", use = "complete.obs"),
            mae = mean(abs(true_prop - est_prop), na.rm = TRUE)) %>% 
  ungroup()
all_df2 <- all_df %>% 
  left_join(summary_df, by = "disease")

ggplot(all_df2, aes(x = true_prop, y = est_prop, col = cell_type, shape = cell_type)) +
  geom_point(size = 1.5) +
  # geom_text(aes(label = paste0("Pearson's correlation = ", pearson)), x = 0.6, y = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, col = "blue") + 
  scale_color_manual(name = "Cell type", values = c("#D72638", "#3F88C5", "#F49D37", "#140F2D")) +
  scale_shape_discrete(name = "Cell type") +
  coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_wrap(~disease) +
  theme_linedraw() +
  labs(x = "Single cell proportion", y = "Estimated proportion") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical") 

## can't know lymphopenia bc not whole blood
# P_long %>% 
#   filter(Cell_type == "T" | Cell_type == "B" | Cell_type == "NK") %>%
#   group_by(subject_id) %>% 
#   summarise(Lymph_prop = sum(Proportion)) %>% 
#   left_join(col_data, by = "subject_id") %>% 
#   ggplot(aes(x = disease_state, y = Lymph_prop, col = disease_state)) +
#   geom_boxplot() +
#   labs(title = "m = 50", x = "Cell type")

all_prop2 <- get_all_prop(P_nnls)


all_est <- all_prop2 %>% 
  as_tibble(rownames = "cell_type") %>% 
  pivot_longer(-cell_type, 
               names_to = "subject_id", 
               values_to = "est_prop") %>% 
  mutate(subject_id = str_split_i(subject_id, "_", 1))

## combine the truths and the estimates
all_df <- all_true %>% 
  left_join(all_est, by = c("subject_id", "cell_type"))
all_df$cell_type <- factor(all_df$cell_type, levels = c("T cells", "B cells", "NK cells", "Monocytes"))

## calculate summary statistics
summary_df <- all_df %>% 
  group_by(disease) %>% 
  summarise(pearson = cor(true_prop, est_prop, use = "complete.obs"),
            spearman = cor(true_prop, est_prop, "spearman", use = "complete.obs"),
            mae = mean(abs(true_prop - est_prop), na.rm = TRUE)) %>% 
  ungroup()
all_df2 <- all_df %>% 
  left_join(summary_df, by = "disease")

ggplot(all_df2, aes(x = true_prop, y = est_prop, col = cell_type, shape = cell_type)) +
  geom_point(size = 1.5) +
  # geom_text(aes(label = paste0("Pearson's correlation = ", pearson)), x = 0.6, y = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, col = "blue") + 
  scale_color_manual(name = "Cell type", values = c("#D72638", "#3F88C5", "#F49D37", "#140F2D")) +
  scale_shape_discrete(name = "Cell type") +
  coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_wrap(~disease) +
  theme_linedraw() +
  labs(x = "Single cell proportion", y = "Estimated proportion") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical") 


# Incorporate Gamma_2 into estimation -------------------------------------

# P_ls6$Gamma_hat[, 2] %>% max

Gamma2 <- P_ls6$Gamma_hat[, 2]

X <- cbind(Gamma2, 1, sig_marker)
colnames(X)[2] <- "intercept"

B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(sig_marker))$x})
P_ext <- B_hat[-seq(2), ]
P_ext_hat <- apply(P_ext, 2, function(x) x/sum(x))
rownames(P_ext_hat) <- colnames(sig_marker)

P_long_ext <- get_long_form(P_ext_hat) %>% 
  mutate(adjust = "Adjusted (second comp)")

P_long_ext$disease_state_bin <- factor(P_long_ext$disease_state_bin, levels = c("Healthy/Convalescent", "COVID-19"))

P_long_ext %>% 
  ggplot(aes(x = Cell_type, y = Proportion, fill = disease_state_bin)) +
  geom_boxplot() +
  labs(title = "Deconvoluted PBMC proportions", 
       subtitle = "Incorporate Gamma2 from dSVA",
       x = "Cell type",
       caption = paste("Using ABIS-seq as reference; q =", 1)) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background = element_blank(),
        legend.position = "bottom",
        # legend.position = c(0.70, 0.98),
        # legend.direction = "horizontal",
        # legend.justification = c(0.01, 0.98),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5)) +
  ggsci::scale_color_tron(name = "Disease state")
ggsave("plots/COVID_19/dSVA_Gamma2_estimate_20240228.pdf")

P_long_nnls <- get_long_form(P_ls$P_hat) %>% 
  mutate(adjust = "Unadjusted")

P_long_nnls %>%
  ggplot(aes(x = Cell_type, y = Proportion, fill = disease_state_bin)) +
  geom_boxplot() +
  labs(title = "Deconvoluted PBMC proportions",
       subtitle = "NNLS Estimates",
       x = "Cell type",
       caption = paste("Using ABIS-seq as reference; q =", 0)) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background = element_blank(),
        legend.position = "bottom",
        # legend.position = c(0.70, 0.98),
        # legend.direction = "horizontal",
        # legend.justification = c(0.01, 0.98),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5)) +
  ggsci::scale_color_tron(name = "Disease state")
ggsave("plots/COVID_19/NNLS_estimate_20240228.pdf")

## Contrast all proportions
P_null_long <- get_long_form(P_ls$P_hat) %>% 
  dplyr::rename(Proportion_unadjusted = Proportion)
P_adjust_long <- get_long_form(P_ext_hat) %>% 
  dplyr::rename(Proportion_adjusted = Proportion)

P_all <- tibble(subject_id = P_null_long$subject_id,
                prop_unadj = P_null_long$Proportion_unadjusted,
                prop_adj = P_adjust_long$Proportion_adjusted,
                severity = P_null_long$severity,
                disease_state = P_null_long$disease_state,
                disease_state_bin = P_null_long$disease_state_bin,
                cell_type = P_null_long$Cell_type)

ggplot(P_all, aes(x = prop_unadj, y = prop_adj, col = disease_state_bin)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2, col = "navy", alpha = .8) +
  # coord_fixed(ratio = 1) +
  labs(x = "Unadjusted proportion (q = 0)", y = paste0("Adjusted proportion (incorporating gamma2)")) +
  facet_wrap(~cell_type, scales = "free") +
  scale_color_manual(name = "Disease State", values = c("#8a508f", "#ffd380")) +
  # ggsci::scale_color_tron(name = "Disease state") +
  # theme_linedraw() +
  theme_base(base_size = 10, base_family = "Times") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom",
        legend.box = "horizontal",
        # panel.grid = element_line(color = "grey"),
        plot.background=element_blank()) 
ggsave("plots/COVID_19/cell_subset_comparison_disease_new.pdf")

ggplot(P_all, aes(x = prop_unadj, y = prop_adj, col = severity)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, col = "navy", alpha = .8) +
  # coord_fixed(ratio = 1) +
  labs(x = "Unadjusted proportion (q = 0)", y = paste0("Adjusted proportion (incorporating gamma2)")) +
  facet_wrap(~cell_type, scales = "free") +
  # ggsci::scale_color_frontiers(name = "Severity")+
  theme_base(base_size = 10, base_family = "Times") +
  scale_color_manual(name = "Severity", values = col_vals[5:1]) +
  # theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom",
        legend.box = "horizontal",
        plot.background=element_blank()) 
ggsave("plots/COVID_19/cell_subset_comparison_sev_new.pdf")


## combine the truths and the estimates
# all_df <- all_true %>% 
#   left_join(all_est, by = c("subject_id", "cell_type"))
# all_df$cell_type <- factor(all_df$cell_type, levels = c("T cells", "B cells", "NK cells", "Monocytes"))

## calculate summary statistics
# summary_df <- all_df %>% 
#   group_by(disease) %>% 
#   summarise(pearson = cor(true_prop, est_prop, use = "complete.obs"),
#             spearman = cor(true_prop, est_prop, "spearman", use = "complete.obs"),
#             mae = mean(abs(true_prop - est_prop), na.rm = TRUE)) %>% 
#   ungroup()
# all_df2 <- all_df %>% 
#   left_join(summary_df, by = "disease")

# ggplot(all_df2, aes(x = true_prop, y = est_prop, col = cell_type, shape = cell_type)) +
#   geom_point(size = 1.5) +
#   # geom_text(aes(label = paste0("Pearson's correlation = ", pearson)), x = 0.6, y = 0.2) +
#   geom_abline(slope = 1, intercept = 0, linetype = 2, col = "blue") + 
#   scale_color_manual(name = "Cell type", values = c("#D72638", "#3F88C5", "#F49D37", "#140F2D")) +
#   scale_shape_discrete(name = "Cell type") +
#   coord_fixed() +
#   xlim(0, 1) +
#   ylim(0, 1) +
#   facet_wrap(~disease) +
#   theme_linedraw() +
#   labs(x = "Single cell proportion", y = "Estimated proportion") +
#   theme(legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.box = "vertical") 

## can't know lymphopenia bc not whole blood
# P_long %>% 
#   filter(Cell_type == "T" | Cell_type == "B" | Cell_type == "NK") %>%
#   group_by(subject_id) %>% 
#   summarise(Lymph_prop = sum(Proportion)) %>% 
#   left_join(col_data, by = "subject_id") %>% 
#   ggplot(aes(x = disease_state, y = Lymph_prop, col = disease_state)) +
#   geom_boxplot() +
#   labs(title = "m = 50", x = "Cell type")

get_all_prop2 <- function(P) {
  T_ct_ind <- grepl("T", rownames(P))
  B_ct_ind <- grepl("B", rownames(P))
  NK_ct_ind <- grepl("NK", rownames(P))
  Mono_ct_ind <- grepl("Mono", rownames(P))
  P_T <- P[T_ct_ind, ]
  P_B <- P[B_ct_ind, ]
  P_NK <- P[NK_ct_ind, ]
  P_Mono <- P[Mono_ct_ind, ]
  all_prop <- rbind(colSums(P_T), colSums(P_B), P_NK, colSums(P_Mono))
  rownames(all_prop) <- c("T cells", "B cells", "NK cells", "Monocytes")
  return(all_prop)
}

## process the data
all_prop_ext <- get_all_prop2(P_ext_hat)
all_est <- all_prop_ext %>% 
  as_tibble(rownames = "cell_type") %>% 
  pivot_longer(-cell_type, 
               names_to = "subject_id", 
               values_to = "est_prop") %>% 
  mutate(subject_id = str_split_i(subject_id, "_", 1))
## combine the truths and the estimates
all_df_ext <- all_true %>% 
  left_join(all_est, by = c("subject_id", "cell_type"))
all_df_ext$cell_type <- factor(all_df_ext$cell_type, levels = c("T cells", "B cells", "NK cells", "Monocytes"))

ggplot(all_df_ext, aes(x = true_prop, y = est_prop, col = cell_type, shape = cell_type)) +
  geom_point(size = 1.5) +
  # geom_text(aes(label = paste0("Pearson's correlation = ", pearson)), x = 0.6, y = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, col = "blue") + 
  scale_color_manual(name = "Cell type", values = c("#D72638", "#3F88C5", "#F49D37", "#140F2D")) +
  scale_shape_discrete(name = "Cell type") +
  coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_wrap(~disease) +
  theme_base(base_size = 12, base_family = "Times") +
  # theme_linedraw() +
  labs(title = "Estimated vs. Single Cell Proportions",
       subtitle = "Incorporate Gamma2 from dSVA",
       x = "Single cell proportion", y = "Estimated proportion") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        plot.background = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5))
ggsave("plots/COVID_19/dSVA_Gamma2_compare_20240228.pdf")

## overall correlation
all_df_ext %>% 
  group_by(disease) %>% 
  summarise(CCC = DescTools::CCC(true_prop, est_prop, na.rm = TRUE)$rho.c[1],
            # MSE = mean((true_prop - est_prop)^2, na.rm = TRUE),
            MAE = mean(abs(true_prop - est_prop), na.rm = TRUE))

## correlation by cell type
Ext_by_CT <- all_df_ext %>% 
  group_by(cell_type, disease) %>% 
  summarise(CCC = DescTools::CCC(true_prop, est_prop, na.rm = TRUE)$rho.c[1],
            # MSE = mean((true_prop - est_prop)^2, na.rm = TRUE),
            MAE = mean(abs(true_prop - est_prop), na.rm = TRUE)) 

## correlation by samples
Ext_by_prop <- all_df_ext %>% 
  group_by(subject_id, disease) %>% 
  summarise(CCC = DescTools::CCC(true_prop, est_prop, na.rm = TRUE)$rho.c[1],
            # MSE = mean((true_prop - est_prop)^2, na.rm = TRUE),
            MAE = mean(abs(true_prop - est_prop), na.rm = TRUE)) %>% 
  dplyr::filter(subject_id != "S155")

## overall correlation - NNLS
all_df_nnls %>% 
  group_by(disease) %>% 
  summarise(CCC = DescTools::CCC(true_prop, est_prop, na.rm = TRUE)$rho.c[1],
            # MSE = mean((true_prop - est_prop)^2, na.rm = TRUE),
            MAE = mean(abs(true_prop - est_prop), na.rm = TRUE))

## correlation by cell type - NNLS
NNLS_by_CT <- all_df_nnls %>% 
  group_by(cell_type, disease) %>% 
  summarise(CCC = DescTools::CCC(true_prop, est_prop, na.rm = TRUE)$rho.c[1],
            # MSE = mean((true_prop - est_prop)^2, na.rm = TRUE),
            MAE = mean(abs(true_prop - est_prop), na.rm = TRUE))

## correlation by sbuject - NNLS
NNLS_by_prop <- all_df_nnls %>% 
  group_by(subject_id, disease) %>% 
  summarise(CCC = DescTools::CCC(true_prop, est_prop, na.rm = TRUE)$rho.c[1],
            # MSE = mean((true_prop - est_prop)^2, na.rm = TRUE),
            MAE = mean(abs(true_prop - est_prop), na.rm = TRUE)) %>% 
  dplyr::filter(subject_id != "S155")

## combine the tables
## by CT
all_by_CT <- NNLS_by_CT %>% 
  right_join(Ext_by_CT, by = c("disease", "cell_type")) %>% 
  arrange(disease)
colnames(all_by_CT) <- c("Cell Type", "Disease Status", "CCC (NNLS Only)", "MAE (NNLS Only)", "CCC (Adjusted)", "MAE (Adjusted)")
print(xtable::xtable(all_by_CT, type = "latex")) # print the latex code

## by sample
all_by_prop <- NNLS_by_prop %>% 
  right_join(Ext_by_prop, by = c("disease", "subject_id")) %>% 
  arrange(disease)
colnames(all_by_prop) <- c("Subject ID", "Disease Status", "CCC (NNLS Only)", "MAE (NNLS Only)", "CCC (Adjusted)", "MAE (Adjusted)")
print(xtable::xtable(all_by_prop, type = "latex")) # print the latex code

sev_col <- bulk_col2 %>% 
  dplyr::select(sample, severity) %>% 
  dplyr::rename(`Subject ID` = sample)
sev_tbl <- all_by_prop %>% 
  left_join(sev_col, by = "Subject ID") %>% 
  dplyr::select(c(1, 7, 3:6)) 
print(xtable::xtable(sev_tbl, type = "latex"))

## Do this again for NNLS only
## process the data
all_prop_nnls <- get_all_prop2(P_ls$P_hat)
all_est <- all_prop_nnls %>% 
  as_tibble(rownames = "cell_type") %>% 
  pivot_longer(-cell_type, 
               names_to = "subject_id", 
               values_to = "est_prop") %>% 
  mutate(subject_id = str_split_i(subject_id, "_", 1))
## combine the truths and the estimates
all_df_nnls <- all_true %>% 
  left_join(all_est, by = c("subject_id", "cell_type"))
all_df_nnls$cell_type <- factor(all_df_nnls$cell_type, levels = c("T cells", "B cells", "NK cells", "Monocytes"))

ggplot(all_df_nnls, aes(x = true_prop, y = est_prop, col = cell_type, shape = cell_type)) +
  geom_point(size = 1.5) +
  # geom_text(aes(label = paste0("Pearson's correlation = ", pearson)), x = 0.6, y = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, col = "blue") + 
  scale_color_manual(name = "Cell type", values = c("#D72638", "#3F88C5", "#F49D37", "#140F2D")) +
  scale_shape_discrete(name = "Cell type") +
  coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_wrap(~disease) +
  theme_base(base_size = 12, base_family = "Times") +
  # theme_linedraw() +
  labs(title = "Estimated vs. Single Cell Proportions",
       subtitle = "NNLS Estimates",
       x = "Single cell proportion", y = "Estimated proportion") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        plot.background = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5))
ggsave("plots/COVID_19/NNLS_compare_20240228.pdf")

P_long_all <- rbind(P_long_ext, P_long_nnls)
p_plasma2 <- P_long_all %>%
  filter(Cell_type == "Plasmablasts") %>% 
  ggplot(aes(x = disease_state_bin, y = Proportion, fill = disease_state_bin)) +
  geom_boxplot() +
  facet_grid(.~adjust) +
  labs(title = "Plasmablast Abundance", x = "Disease state"
       # ,caption = paste("Using ABIS-seq as reference; q_hat =", q_hat2)
  ) +
  theme_base(base_size = 12, base_family = "Times") +
  theme(axis.text.x = element_blank(),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5),
        legend.position = "bottom",
        plot.background=element_blank()) +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f")) 
p_plasma_ext <- p_plasma2 + ggpubr::stat_compare_means(method = "wilcox.test")
p_plasma_ext
ggsave("plots/COVID_19/Plasmablast_estimate_20240228.pdf")

p_CD82 <- P_long_all %>%
  filter(Cell_type == "T.CD8.Memory") %>%
  # group_by(subject_id) %>%
  # summarise(Proportion = sum(Proportion)) %>%
  # right_join(P_long %>% dplyr::select(-Proportion) %>% distinct()) %>%
  ggplot(aes(x = disease_state_bin, y = Proportion, fill = disease_state_bin)) +
  geom_boxplot() +
  facet_grid(.~adjust) +
  labs(title = "Memory CD8+ T Cell Abundances", x = "Disease state"
       # ,caption = paste("Using ABIS-seq as reference; q_hat =", q_hat2)
  ) +
  theme_base(base_size = 12, base_family = "Times") +
  theme(axis.text.x = element_blank(),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5),
        legend.position = "bottom",
        plot.background=element_blank()) +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f")) 
p_CD8_ext <- p_CD82 + ggpubr::stat_compare_means(method = "wilcox.test")
p_CD8_ext
ggsave("plots/COVID_19/Memory_CD8_estimate_20240228.pdf")

p_NK <- P_long_all %>%
  filter(Cell_type == "NK") %>%
  ggplot(aes(x = disease_state_bin, y = Proportion, fill = disease_state_bin)) +
  geom_boxplot() +
  facet_grid(.~adjust) +
  labs(title = "NK Cell Abundances", x = "Disease state") +
  theme_base(base_size = 12, base_family = "Times") +
  theme(axis.text.x = element_blank(),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5),
        legend.position = "bottom",
        plot.background=element_blank()) +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f")) 
p_NK_ext <- p_NK + ggpubr::stat_compare_means(method = "wilcox.test")
p_NK_ext

P_long_all %>%
  filter(grepl("CD8", Cell_type)) %>%
  filter(post_infection3 < 30) %>% 
  # group_by(subject_id) %>%
  # summarise(Proportion = sum(Proportion)) %>%
  # right_join(P_long %>% select(-Proportion) %>% distinct()) %>%
  # ungroup() %>% 
  ggplot(aes(x = post_infection3, y = Proportion)) +
  geom_point(aes(col = severity)) +
  xlim(0, 25) +
  geom_smooth() +
  facet_grid(adjust~Cell_type) +
  labs(title = "CD8+ T Cell Abundance",
       x = "Days post infection",
       caption = paste("Using ABIS-seq as reference; q_hat =", 1, "; convalescent patient filtered")) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_color_manual(name = "Severity", values = col_vals[5:1]) +
  theme(legend.position = "bottom",
        plot.background=element_blank())
ggsave("plots/COVID_19/Memory_CD8_time_20240228.pdf")

P_long_all %>%
  filter(grepl("CD4", Cell_type)) %>%
  filter(post_infection3 < 30) %>% 
  # group_by(subject_id) %>%
  # summarise(Proportion = sum(Proportion)) %>%
  # right_join(P_long %>% select(-Proportion) %>% distinct()) %>%
  # ungroup() %>% 
  ggplot(aes(x = post_infection3, y = Proportion)) +
  geom_point(aes(col = severity)) +
  xlim(0, 25) +
  geom_smooth() +
  facet_grid(adjust~Cell_type) +
  labs(title = "CD4+ T Cell Abundance",
       x = "Days post infection",
       caption = paste("Using ABIS-seq as reference; q_hat =", 1, "; convalescent patient filtered")) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_color_manual(name = "Severity", values = col_vals[5:1]) +
  theme(legend.position = "bottom",
        plot.background=element_blank())

P_long_all %>%
  filter(grepl("Mono", Cell_type)) %>%
  filter(post_infection3 < 30) %>% 
  # group_by(subject_id) %>%
  # summarise(Proportion = sum(Proportion)) %>%
  # right_join(P_long %>% select(-Proportion) %>% distinct()) %>%
  # ungroup() %>% 
  ggplot(aes(x = post_infection3, y = Proportion)) +
  geom_point(aes(col = severity)) +
  xlim(0, 25) +
  geom_smooth() +
  facet_grid(adjust~Cell_type) +
  labs(title = "Monocyte Abundance",
       x = "Days post infection",
       caption = paste("Using ABIS-seq as reference; q_hat =", 1, "; convalescent patient filtered")) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_color_manual(name = "Severity", values = col_vals[5:1]) +
  theme(legend.position = "bottom",
        plot.background=element_blank())

p_T_NNLS <- P_long_all %>%
  filter(grepl("T", Cell_type)) %>%
  filter(post_infection3 < 30) %>% 
  filter(adjust == "Unadjusted") %>% 
  # group_by(subject_id) %>%
  # summarise(Proportion = sum(Proportion)) %>%
  # right_join(P_long %>% select(-Proportion) %>% distinct()) %>%
  # ungroup() %>% 
  ggplot(aes(x = post_infection3, y = Proportion)) +
  geom_point(aes(col = severity)) +
  xlim(0, 25) +
  geom_smooth(method = "lm") +
  facet_wrap(Cell_type~., scales = "free") +
  labs(title = "T Cell Abundance",
       subtitle = "NNLS Estimates",
       x = "Days post infection") +
  theme_base(base_size = 12, base_family = "Times") +
  scale_color_manual(name = "Severity", values = col_vals[5:1]) +
  theme(legend.position = "bottom",
        plot.background=element_blank())

p_T_ext <- P_long_all %>%
  filter(grepl("T", Cell_type)) %>%
  filter(post_infection3 < 30) %>% 
  filter(adjust != "Unadjusted") %>% 
  # group_by(subject_id) %>%
  # summarise(Proportion = sum(Proportion)) %>%
  # right_join(P_long %>% select(-Proportion) %>% distinct()) %>%
  # ungroup() %>% 
  ggplot(aes(x = post_infection3, y = Proportion)) +
  geom_point(aes(col = severity)) +
  xlim(0, 25) +
  geom_smooth(method = "lm") +
  facet_wrap(Cell_type~., scales = "free") +
  labs(subtitle = "dSVA Estimates with Gamma2 Incorporated",
       x = "Days post infection",
       caption = paste("Using ABIS-seq as reference; q_hat =", 1, "; convalescent patient filtered")) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_color_manual(name = "Severity", values = col_vals[5:1]) +
  theme(legend.position = "bottom",
        plot.background=element_blank())

ggarrange(p_T_NNLS, p_T_ext, ncol = 1, common.legend = TRUE)
ggsave("plots/COVID_19/All_T_time_20240228.pdf", width = 7, height = 12)

# Analysis by Projecting The Residual Space -------------------------------

cor(sig_marker, P_ls2$Gamma_hat)
get_proj <- function(x) {
  if (is.vector(x)) x %*% t(x) / (sum(x^2)) else x %*% solve(t(x) %*% x) %*% t(x)  
}

M_cm <- get_proj(mono_c_sig)

Proj_mat <- diag(1, length(mono_c_sig), length(mono_c_sig)) - M_cm

Gamma_new <- Proj_mat %*% P_ls2$Gamma_hat 
cor(Gamma_new, sig_marker)


X <- cbind(Gamma_new, 1, sig_marker)
colnames(X)[1:3] <- c("Gh1", "Gh2", "intercept")

B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(sig_marker))$x})
P_ext2 <- B_hat[-seq(3), ]
P_ext_hat2 <- apply(P_ext2, 2, function(x) x/sum(x))
rownames(P_ext_hat2) <- colnames(sig_marker)

P_long_ext2 <- get_long_form(P_ext_hat2) %>% 
  mutate(adjust = "Adjusted (second comp)")

P_long_ext2 %>% 
  ggplot(aes(x = Cell_type, y = Proportion, fill = disease_state_bin)) +
  geom_boxplot() +
  labs(title = "Deconvoluted PBMC proportions", 
       subtitle = "Latent Factors Projected onto The Ortho. Space of Monocyte.C",
       x = "Cell type",
       caption = paste("Using ABIS-seq as reference; q =", P_ls2$q_hat)) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background = element_blank(),
        legend.position = "bottom",
        # legend.position = c(0.70, 0.98),
        # legend.direction = "horizontal",
        # legend.justification = c(0.01, 0.98),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5)) +
  ggsci::scale_color_tron(name = "Disease state")

all_prop_ext2 <- get_all_prop2(P_ext_hat2)
all_est2 <- all_prop_ext2 %>% 
  as_tibble(rownames = "cell_type") %>% 
  pivot_longer(-cell_type, 
               names_to = "subject_id", 
               values_to = "est_prop") %>% 
  mutate(subject_id = str_split_i(subject_id, "_", 1))
## combine the truths and the estimates
all_df_ext2 <- all_true %>% 
  left_join(all_est2, by = c("subject_id", "cell_type"))
all_df_ext2$cell_type <- factor(all_df_ext2$cell_type, levels = c("T cells", "B cells", "NK cells", "Monocytes"))

ggplot(all_df_ext2, aes(x = true_prop, y = est_prop, col = cell_type, shape = cell_type)) +
  geom_point(size = 1.5) +
  # geom_text(aes(label = paste0("Pearson's correlation = ", pearson)), x = 0.6, y = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, col = "blue") + 
  scale_color_manual(name = "Cell type", values = c("#D72638", "#3F88C5", "#F49D37", "#140F2D")) +
  scale_shape_discrete(name = "Cell type") +
  coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_wrap(~disease) +
  theme_base(base_size = 12, base_family = "Times") +
  # theme_linedraw() +
  labs(title = "Estimated vs. Single Cell Proportions",
       subtitle = "Latent Factors Projected onto The Ortho. Space of Monocyte.C",
       x = "Single cell proportion", y = "Estimated proportion") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        plot.background = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5))

all_df_ext2 %>% 
  group_by(cell_type, disease) %>% 
  summarise(CCC = DescTools::CCC(true_prop, est_prop, na.rm = TRUE)$rho.c[1],
            MSE = mean((true_prop - est_prop)^2, na.rm = TRUE),
            MAE = mean(abs(true_prop - est_prop), na.rm = TRUE))


# Projecting onto the Orthog. Space of Theta ------------------------------

M_theta <- get_proj(sig_marker)

Proj_mat2 <- diag(1, length(mono_c_sig), length(mono_c_sig)) - M_theta

Gamma_new2 <- Proj_mat2 %*% P_ls2$Gamma_hat 
cor(Gamma_new2, sig_marker)
X <- cbind(Gamma_new2, 1, sig_marker)
colnames(X)[1:3] <- c("Gh1", "Gh2", "intercept")

B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(sig_marker))$x})
P_ext3 <- B_hat[-seq(3), ]
P_ext_hat3 <- apply(P_ext3, 2, function(x) x/sum(x))
rownames(P_ext_hat3) <- colnames(sig_marker)

P_long_ext3 <- get_long_form(P_ext_hat3) %>% 
  mutate(adjust = "Adjusted (Theta Perp)")

P_long_ext3 %>% 
  ggplot(aes(x = Cell_type, y = Proportion, fill = disease_state_bin)) +
  geom_boxplot() +
  labs(title = "Deconvoluted PBMC proportions", 
       subtitle = "Latent Factors Projected onto Theta Perp",
       x = "Cell type",
       caption = paste("Using ABIS-seq as reference; q =", P_ls2$q_hat)) +
  theme_base(base_size = 12, base_family = "Times") +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background = element_blank(),
        legend.position = "bottom",
        # legend.position = c(0.70, 0.98),
        # legend.direction = "horizontal",
        # legend.justification = c(0.01, 0.98),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5)) +
  ggsci::scale_color_tron(name = "Disease state")

all_prop_ext3 <- get_all_prop2(P_ext_hat3)
all_est3 <- all_prop_ext3 %>% 
  as_tibble(rownames = "cell_type") %>% 
  pivot_longer(-cell_type, 
               names_to = "subject_id", 
               values_to = "est_prop") %>% 
  mutate(subject_id = str_split_i(subject_id, "_", 1))
## combine the truths and the estimates
all_df_ext3 <- all_true %>% 
  left_join(all_est3, by = c("subject_id", "cell_type"))
all_df_ext3$cell_type <- factor(all_df_ext3$cell_type, levels = c("T cells", "B cells", "NK cells", "Monocytes"))

ggplot(all_df_ext3, aes(x = true_prop, y = est_prop, col = cell_type, shape = cell_type)) +
  geom_point(size = 1.5) +
  # geom_text(aes(label = paste0("Pearson's correlation = ", pearson)), x = 0.6, y = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, col = "blue") + 
  scale_color_manual(name = "Cell type", values = c("#D72638", "#3F88C5", "#F49D37", "#140F2D")) +
  scale_shape_discrete(name = "Cell type") +
  coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_wrap(~disease) +
  theme_base(base_size = 12, base_family = "Times") +
  # theme_linedraw() +
  labs(title = "Estimated vs. Single Cell Proportions",
       subtitle = "Latent Factors Projected onto Theta Perp",
       x = "Single cell proportion", y = "Estimated proportion") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        plot.background = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5))

all_df_ext3 %>% 
  group_by(cell_type, disease) %>% 
  summarise(CCC = DescTools::CCC(true_prop, est_prop, na.rm = TRUE)$rho.c[1],
            MSE = mean((true_prop - est_prop)^2, na.rm = TRUE),
            MAE = mean(abs(true_prop - est_prop), na.rm = TRUE))

# Chuwen's matrix ---------------------------------------------------------

## use the controls' scRNA data from Chuwen to do this again
# # chuwen_full <- readRDS("../dSVA_datasets/Chuwen_sc_sig_full.rds") # signature using all cells
# # chuwen_full <- readRDS("../dSVA_datasets/Chuwen_sc_sig_mid.rds") # signature using cells whose total counts are in the middle
# chuwen_full <- readRDS("../dSVA_datasets/Chuwen_sc_sig_sum_mid_500_cells.rds") # signature using sum of expression from cells whose total counts are in the middle
# # chuwen_full <- readRDS("../dSVA_datasets/Chuwen_sc_sig_sum_mid_500_cells_control.rds") # CONTROLS: healthy sample signature using sum of expression from cells whose total counts are in the middle
# 
# chuwen_full <- chuwen_full[, 1:4] # remove DCs
# 
# chuwen_full2 <- chuwen_full[rowSums(chuwen_full) > 5, ]
# overlap_genes2 <- intersect(rownames(bulk_full), rownames(chuwen_full2)) # 1,161 genes
# bulk_sub2 <- bulk_full[overlap_genes2, ]
# sig_sub2 <- chuwen_full2[overlap_genes2, ]
# 
# ## check the signature matrix further
# dim(sig_sub2)
# dim(sig_full2)
# 
# cor(sig_sub2[grepl("CD4", rownames(sig_sub2)), "T"], sig_sub2[grepl("CD4", rownames(sig_sub2)), "NK"])
# 
# sig_sub2[grepl("CD3", rownames(sig_sub2)), ] # high on T cell
# sig_sub2[grepl("CD4", rownames(sig_sub2)), ] # high on monocyte (not normal)
# sig_sub2[grepl("CD8", rownames(sig_sub2)), ] # high on T cell
# sig_sub2[grepl("CD14", rownames(sig_sub2)), ] # high on monocyte
# 
# sig_sub2[grepl("CD38", rownames(sig_sub2)), ] # high on NK in the healthy reference (should be on B?) - High on B in the COVID-19 reference!
# sig_sub2[grepl("CD19", rownames(sig_sub2)), ] # high on B 
# 
# 
# ## find the marker genes (100 per cell type)
# rat_df2 <- tibble(gene = rownames(sig_sub2),
#                  cell_type = apply(sig_sub2, 1, function(x) colnames(sig_sub2)[which(x == max(x))[1]]))
# max_val2 <- apply(sig_sub2, 1, max)
# sec_val2 <- apply(sig_sub2, 1, function(x) sort(x, decreasing = TRUE)[2])
# rat_df2 <- rat_df2 %>% 
#   mutate(ratio = max_val2/sec_val2)
# marker_genes2 <- rat_df2 %>%
#   group_by(cell_type) %>%
#   arrange(desc(ratio)) %>%
#   filter(row_number() <= 200) %>%
#   # filter(ratio >= 7) %>% 
#   ungroup %>%
#   dplyr::select(gene) %>%
#   unlist
# 
# # rat_df2 %>%
# #   group_by(cell_type) %>%
# #   arrange(desc(ratio)) %>%
# #   filter(ratio >= 7) %>%
# #   group_by(cell_type) %>%
# #   summarise(n = n())
# 
# # marker_genes3 <- rat_df2 %>% 
# #   group_by(cell_type) %>% 
# #   arrange(desc(ratio)) %>% 
# #   filter(ratio >= 4) %>% 
# #   ungroup %>% 
# #   dplyr::select(gene) %>% 
# #   unlist
# 
# bulk_marker2 <- bulk_sub2[marker_genes2, ]
# sig_marker2 <- sig_sub2[marker_genes2, ]
# cor(sig_marker2)
# # sig_marker2 <- sig_marker[, !grepl("DC", colnames(sig_marker))]
# 
# ## run dSVA deconv
# q_hat3 <- estimate_n_comp(Y = bulk_marker2, Theta = sig_marker2, intercept = TRUE, method = "be", B = 49, seed = 100) 
# # q_hat2 <- estimate_n_comp(Y = bulk_marker2, Theta = sig_marker2, intercept = TRUE, method = "cutoff")
# 
# # q_hat2 <- estimate_n_comp(Y = bulk_marker, Theta = sig_marker, intercept = FALSE, method = "be", B = 49, seed = 100)
# 
# P_nnls3 <- dsva_for_sim(Y = bulk_marker2, 
#                         Theta = sig_marker2, 
#                         n_comp = q_hat3, 
#                         intercept = TRUE, 
#                         alg = "nnls", 
#                         solver = "lsei") # using pnnls the results get interesting
# 
# rownames(P_nnls3) <- colnames(sig_marker2)
# colnames(P_nnls3) <- colnames(bulk_marker2)
# 
# P_df2 <- as_tibble(t(P_nnls3), rownames = "subject_id")
# P_long2 <- P_df2 %>% 
#   pivot_longer(-subject_id, 
#                names_to = "Cell_type",
#                values_to = "Proportion") %>% 
#   left_join(bulk_cols, by = "subject_id")
# P_long2 %>% 
#   ggplot(aes(x = Cell_type, y = Proportion, col = disease_state)) +
#   geom_boxplot() +
#   ggpubr::theme_pubr() +
#   labs(title = "m = 1600 (COVID-19 signature from Chuwen), Remove 2 Samples \n Using Middle Cells as References", x = "Cell type")
# 
# 
# P_long2$severity <- factor(P_long2$severity, levels = c("Healthy",
#                                                         "Convalescent",
#                                                         "Moderate",
#                                                         "Severe",
#                                                         "ICU"))
# P_long2 %>% 
#   filter(Cell_type == "T") %>% 
#   ggplot(aes(x = severity, y = Proportion, col = severity)) +
#   geom_boxplot(alpha = 0.4) +
#   geom_jitter() +
#   labs(title = "T-cell Proportion, m = 1600 (COVID-19 signature from Chuwen),  \n Remove 2 Samples \n Using Middle Cells as References", x = "Disease Severity") +
#   scale_color_viridis_d() +
#   ggpubr::theme_pubr() 
# 
# P_long2 %>% 
#   filter(Cell_type == "Mono") %>% 
#   ggplot(aes(x = severity, y = Proportion, col = severity)) +
#   geom_boxplot(alpha = 0.4) +
#   geom_jitter() +
#   labs(title = "Monocyte Proportion, m = 1600 (COVID-19 signature from Chuwen),  \n Remove 2 Samples \n Using Middle Cells as References", x = "Disease Severity") +
#   scale_color_viridis_d() +
#   ggpubr::theme_pubr() 
# 
# # P_long2 %>% 
# #   filter(Cell_type == "T") %>% 
# #   ggplot(aes(x = post_infection, y = Proportion, col = severity)) +
# #   geom_jitter() +
# #   labs(title = "Selected ratio >= 4", x = "Disease Severity") +
# #   scale_color_viridis_d() +
# #   ggpubr::theme_pubr() 
# 
# ## draw the comparison plots
# P_null <- dsva_for_sim(Y = bulk_marker2, 
#                        Theta = sig_marker2, 
#                        n_comp = 0, 
#                        intercept = TRUE, 
#                        alg = "nnls", 
#                        solver = "lsei") # using pnnls the results get interesting
# 
# P_adjust <- dsva_for_sim(Y = bulk_marker2, 
#                          Theta = sig_marker2, 
#                          n_comp = q_hat3, 
#                          intercept = TRUE, 
#                          alg = "nnls", 
#                          solver = "lsei")
# 
# rownames(P_null) <- rownames(P_adjust) <- colnames(sig_marker2)
# colnames(P_null) <- colnames(P_adjust) <- colnames(bulk_marker2)
# 
# P_null_long <- get_long_form(P_null) %>% 
#   rename(Proportion_unadjusted = Proportion)
# P_adjust_long <- get_long_form(P_adjust) %>% 
#   rename(Proportion_adjusted = Proportion)
# 
# P_all <- tibble(subject_id = P_null_long$subject_id,
#                 prop_unadj = P_null_long$Proportion_unadjusted,
#                 prop_adj = P_adjust_long$Proportion_adjusted,
#                 severity = P_null_long$severity,
#                 disease_state = P_null_long$disease_state,
#                 cell_type = P_null_long$Cell_type)
# 
# ggplot(P_all, aes(x = prop_unadj, y = prop_adj, col = severity)) +
#   geom_point(alpha = 0.6) +
#   geom_abline(slope = 1, intercept = 0, linetype = 2) +
#   coord_fixed(ratio = 1) +
#   labs(x = "Unadjusted proportion (q = 0)", y = paste0("Adjusted proportion (q = ", q_hat3, ")")) +
#   facet_wrap(~cell_type) +
#   ggpubr::theme_pubr()


# P_long %>% 
#   filter(Cell_type == "B") %>% 
#   ggplot(aes(x = severity, y = Proportion, col = severity)) +
#   geom_boxplot(alpha = 0.4) +
#   geom_jitter() +
#   labs(title = "B-Cell Proportion (m = 500)", x = "Disease Severity") +
#   scale_color_viridis_d() +
#   ggpubr::theme_pubr()