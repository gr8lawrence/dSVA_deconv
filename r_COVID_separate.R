library(tidyverse)
library(ggpubr)
library(ggfortify)
library(ggthemes)
source("a_dSVA_functions.R")
source("r_real_data_funs.R")
source("r_real_data_support_funs.R")

## NOTE: number of marker per cell type = 300, 350, 400 works! 200, 450, 500 do not
## read the bulk data
bulk_full <- readRDS("../dSVA_datasets/covid_bulk_count_processed.rds")

## gene signature from controls
sig_controls <- readRDS("../dSVA_datasets/Chuwen_sc_sig_sum_mid_500_cells_control.rds")
con_genes <- intersect(rownames(sig_controls), rownames(bulk_full))
sig_con2 <- sig_controls[con_genes, ]

## cases
sig_cases <- readRDS("../dSVA_datasets/Chuwen_sc_sig_sum_mid_500_cells.rds")
case_genes <- intersect(rownames(sig_cases), rownames(bulk_full))
sig_case2 <- sig_cases[case_genes, ]

## bulk column data
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
bulk_cols$sample <- str_split_i(bulk_cols$subject_id, "_", 1)

## getting the single-cell proportions (ground truths)
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
all_true <- all_true %>% 
  mutate(cell_type = case_when(
    cell_type == "T cells" ~ "T",
    cell_type == "B cells" ~ "B",
    cell_type == "NK cells" ~ "NK",
    cell_type == "Monocytes" ~ "Mono",
  ))


# Control Signatures ------------------------------------------------------

# sig_con_300 <- get_marker_genes(sig_con2, 300) # include DC (1,500 genes)
sig_con_300 <- get_marker_genes(sig_con2[, 1:4], 350) # NOT include DC (800 genes)

P_ls_con_300_nnls <- get_P_ls(bulk_full[sig_con_300, ], sig_con2[sig_con_300, 1:4], "none") # remove 1:4 if want to include DC
P_ls_con_300 <- get_P_ls(bulk_full[sig_con_300, ], sig_con2[sig_con_300, 1:4], "cutoff")

## check the correlation
P_nnls <- P_ls_con_300_nnls$P_hat
D_nnls <- cbind(matrix(0, nrow(P_nnls), ncol(P_nnls)/2), P_nnls[, seq(ncol(P_nnls)/2 + 1, ncol(P_nnls))])
cor(t(P_nnls), t(D_nnls))

# P_ls_con_300 <- run_dSVA(bulk_full[sig_con_300, ], 
#                          sig_con2[sig_con_300, 1:3], 
#                          intercept = TRUE, 
#                          exclude = NULL, 
#                          alg = "nnls", 
#                          solver = "lsei", 
#                          q = 2)

cor(sig_con2[sig_con_300, 1:4])
cor(sig_con2[sig_con_300, 1:4] , P_ls_con_300$Gamma_hat)

## get the R hat PCA plots
R_hat_con <- get_R_hat(bulk_full[sig_con_300, ],  sig_con2[sig_con_300, 1:4])
pca_con <- prcomp(t(R_hat_con), scale. = TRUE, center = TRUE)
p_pca_con <- autoplot(pca_con, data = bulk_cols, colour = "disease_state") +
  labs(title = "PCA Plot on COVID-19 Bulk Residuals", 
       subtitle = "Single Cell Data from Healthy Controls") +
  ggrepel::geom_label_repel(aes(label = sample)) +
  theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Disease state", values = c("#8a508f", "#ffd380", "#2c4875")) 
p_pca_con

## proportion boxplots
P_con_300_long <- get_long_form(P_ls_con_300$P_hat)
box_con_300 <- get_prop_boxplots(P_con_300_long, "Single Cell Data from Healthy Controls + dSVA (q = 1)")
box_con_300

P_con_300_NNLS_long <- get_long_form(P_ls_con_300_nnls$P_hat)
box_con_300_NNLS <- get_prop_boxplots(P_con_300_NNLS_long, "Single Cell Data from Healthy Controls + NNLS")
box_con_300_NNLS

## relative proportions of the first four cell types
P_con_300_rel_df <- apply(P_ls_con_300$P_hat[1:4, ], 2, function(x) x/sum(x)) %>% 
  get_prop_df() # dSVA/cutoff + NNLS
P_con_300_nnls_rel_df <- apply(P_ls_con_300_nnls$P_hat[1:4, ], 2, function(x) x/sum(x)) %>% 
  get_prop_df() # NNLS

## combine the truths and the estimates
p_con_300_dSVA <- get_compare_plot(P_con_300_rel_df, all_true, "Single Cell Data from Healthy Controls + dSVA (q = 1)")
p_con_300_dSVA

p_con_300_NNLS <- get_compare_plot(P_con_300_nnls_rel_df, all_true, "Single Cell Data from Healthy Controls + NNLS")
p_con_300_NNLS

# get the stats (promising results only)
stat_con_300_dSVA_ls <- get_compare_stat(P_con_300_rel_df, all_true)
stat_con_300_dSVA_ls$subject
stat_con_300_dSVA_ls$CT


stat_con_300_NNLS_ls <- get_compare_stat(P_con_300_nnls_rel_df, all_true)
stat_con_300_NNLS_ls$subject
stat_con_300_NNLS_ls$CT


# Case signatures ---------------------------------------------------------

# sig_case_300 <- get_marker_genes(sig_case2, 300) #1,500 genes
sig_case_300 <- get_marker_genes(sig_case2[, 1:4], 350) #1,200 genes

## get the R hat PCA plot
R_hat_cases <- get_R_hat(bulk_full[sig_case_300, ],  sig_case2[sig_case_300, 1:4])
pca_cases <- prcomp(t(R_hat_cases), scale. = TRUE, center = TRUE)
p_pca_cases <- autoplot(pca_cases, data = bulk_cols, colour = "disease_state") +
  labs(title = "PCA Plot on COVID-19 Bulk Residuals", 
       subtitle = "Single Cell Data from COVID-19 Subjects") +
  theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Disease state", values = c("#8a508f", "#ffd380", "#2c4875")) 
p_pca_cases


P_ls_case_300_nnls <- get_P_ls(bulk_full[sig_case_300, ], sig_case2[sig_case_300, 1:4], "none")
P_ls_case_300 <- get_P_ls(bulk_full[sig_case_300, ], sig_case2[sig_case_300, 1:4], "cutoff")

cor(sig_case2[sig_case_300, 1:4])
cor(sig_case2[sig_case_300, 1:4] , P_ls_case_300$Gamma_hat)

P_case_300_long <- get_long_form(P_ls_case_300$P_hat)
box_case_300 <- get_prop_boxplots(P_case_300_long, "Single Cell Data from COVID-19 Subjects + dSVA (q = 1)")
box_case_300

P_case_300_NNLS_long <- get_long_form(P_ls_case_300_nnls$P_hat)
box_case_300_NNLS <- get_prop_boxplots(P_case_300_NNLS_long, "Single Cell Data from COVID-19 Subjects + NNLS")
box_case_300_NNLS

## get the relative proportions of the first four cell types
P_case_300_rel_df <- apply(P_ls_case_300$P_hat[1:4, ], 2, function(x) x/sum(x)) %>% 
  get_prop_df() # dSVA/cutoff + NNLS
P_case_300_nnls_rel_df <- apply(P_ls_case_300_nnls$P_hat[1:4, ], 2, function(x) x/sum(x)) %>% 
  get_prop_df() # NNLS

## combine the truths and the estimates
p_case_300_dSVA <- get_compare_plot(P_case_300_rel_df, all_true, "Single Cell Data from COVID-19 Subjects + dSVA (q = 1)")
p_case_300_dSVA

p_case_300_NNLS <- get_compare_plot(P_case_300_nnls_rel_df, all_true, "Single Cell Data from COVID-19 Subjects + NNLS")
p_case_300_NNLS


stat_case_300_NNLS_ls <- get_compare_stat(P_case_300_nnls_rel_df, all_true)
stat_case_300_NNLS_ls$subject
stat_case_300_NNLS_ls$CT

# Scratch code ------------------------------------------------------------

## Show the results
## COVID-19 Subjects
P1 <- P_con_300_rel_df %>% 
  left_join(meta_data, by = "subject_id") %>% 
  filter(disease == "COVID-19") %>% 
  mutate(method = "dSVA + NNLS", ref = "Ref: Healthy") %>% 
  dplyr::select(-sample, -disease) %>% 
  left_join(all_true, by = c("subject_id", "cell_type"))
P2 <- P_con_300_nnls_rel_df %>% 
  left_join(meta_data, by = "subject_id") %>% 
  filter(disease == "COVID-19") %>% 
  mutate(method = "NNLS", ref = "Ref: Healthy") %>% 
  dplyr::select(-sample, -disease) %>% 
  left_join(all_true, by = c("subject_id", "cell_type"))
P3 <- P_case_300_nnls_rel_df %>% 
  left_join(meta_data, by = "subject_id") %>% 
  filter(disease == "COVID-19") %>% 
  mutate(method = "NNLS", ref = "Ref: COVID-19") %>% 
  dplyr::select(-sample, -disease) %>% 
  left_join(all_true, by = c("subject_id", "cell_type"))
P_summary <- rbind(P1, P2, P3)
P_summary$cell_type <- factor(P_summary$cell_type, levels = c("T", "B", "NK", "Mono"))
P_summary$ref <- factor(P_summary$ref, levels = c("Ref: Healthy", "Ref: COVID-19"))
P_summary$method <- factor(P_summary$method, levels = c("NNLS", "dSVA + NNLS"))


ggplot(P_summary, aes(x = true_prop, y = est_prop, col = cell_type, shape = cell_type)) +
  geom_point(size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, col = "blue") + 
  scale_color_manual(name = "Cell type", values = c("#D72638", "#3F88C5", "#F49D37", "#140F2D")) +
  scale_shape_discrete(name = "Cell type") +
  coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_wrap(ref~method) +
  theme_base(base_size = 12, base_family = "Times") +
  labs(title = "Estimated vs. Single Cell Proportions",
       subtitle = "COVID-19 Diagnosed Patients",
       x = "Single cell proportion", y = "Estimated proportion") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        plot.background = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5))
ggsave("plots/COVID_19/Compare_COVID_patients_20240313.pdf", width = 10.5, height = 6.5)

## Healthy subjects
## Show the results
P4 <- P_con_300_rel_df %>% 
  left_join(meta_data, by = "subject_id") %>% 
  filter(disease == "Control") %>% 
  mutate(method = "dSVA + NNLS", ref = "Ref: Healthy") %>% 
  dplyr::select(-sample, -disease) %>% 
  left_join(all_true, by = c("subject_id", "cell_type"))
P5 <- P_con_300_nnls_rel_df %>% 
  left_join(meta_data, by = "subject_id") %>% 
  filter(disease == "Control") %>% 
  mutate(method = "NNLS", ref = "Ref: Healthy") %>% 
  dplyr::select(-sample, -disease) %>% 
  left_join(all_true, by = c("subject_id", "cell_type"))

P_con_summary <- rbind(P4, P5)
P_con_summary$cell_type <- factor(P_con_summary$cell_type, levels = c("T", "B", "NK", "Mono"))
P_con_summary$ref <- factor(P_con_summary$ref, levels = c("Ref: Healthy", "Ref: COVID-19"))
P_con_summary$method <- factor(P_con_summary$method, levels = c("NNLS", "dSVA + NNLS"))


ggplot(P_con_summary, aes(x = true_prop, y = est_prop, col = cell_type, shape = cell_type)) +
  geom_point(size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, col = "blue") + 
  scale_color_manual(name = "Cell type", values = c("#D72638", "#3F88C5", "#F49D37", "#140F2D")) +
  scale_shape_discrete(name = "Cell type") +
  coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_wrap(ref~method) +
  theme_base(base_size = 12, base_family = "Times") +
  labs(title = "Estimated vs. Single Cell Proportions",
       subtitle = "Healthy Controls",
       x = "Single cell proportion", y = "Estimated proportion") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        plot.background = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5))
ggsave("plots/COVID_19/Compare_healthy_patients_20240313.pdf")


# theme_minimal()# all_con_300_rel_df <- all_true %>% 
#   left_join(P_con_300_rel_df , by = c("subject_id", "cell_type"))
# all_con_300_rel_df$cell_type <- factor(all_con_300_rel_df$cell_type, levels = c("T", "B", "NK", "Mono"))
# 
# ggplot(all_con_300_rel_df, aes(x = true_prop, y = est_prop, col = cell_type, shape = cell_type)) +
#   geom_point(size = 1.5) +
#   # geom_text(aes(label = paste0("Pearson's correlation = ", pearson)), x = 0.6, y = 0.2) +
#   geom_abline(slope = 1, intercept = 0, linetype = 2, col = "blue") + 
#   scale_color_manual(name = "Cell type", values = c("#D72638", "#3F88C5", "#F49D37", "#140F2D")) +
#   scale_shape_discrete(name = "Cell type") +
#   coord_fixed() +
#   xlim(0, 1) +
#   ylim(0, 1) +
#   facet_wrap(~disease) +
#   theme_base(base_size = 12, base_family = "Times") +
#   # theme_linedraw() +
#   labs(title = "Estimated vs. Single Cell Proportions",
#        subtitle = "Single Cell Data from Healthy Controls + dSVA (q = 1)",
#        x = "Single cell proportion", y = "Estimated proportion") +
#   theme(legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.box = "vertical",
#         plot.background = element_blank(),
#         legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
#         panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5))
# 
