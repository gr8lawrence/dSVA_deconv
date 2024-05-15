## This file is for analyzing the longitudinal data in GSE52428 separately for H1N1 and H3N2 samples
library(tidyverse)
library(ggfortify)
library(ggrepel)
library(gprofiler2)
library(ggthemes)
source("a_dSVA_functions.R")
source("r_real_data_funs.R")
source("r_real_data_support_funs.R")


# Preprocessing -----------------------------------------------------------

## read the full data
all_data_ls <- readRDS("../dSVA_datasets/Longitudinal_viral_challenge_data.rds")
all_bulk <- all_data_ls$bulk
all_sig <- all_data_ls$sig
cell_types <- c("T cell", "B cell", "NK cell", "Monocyte", "Neutrophil")
all_sig <- all_sig[, cell_types]
scale_factor <- c(1, 0.946, 0.742, 0.756, 0.115) 
all_sig1 <- all_sig %*% diag(1 / scale_factor) # add the scaling factors
colnames(all_sig1) <- colnames(all_sig)
all_meta <- all_data_ls$meta

### separate pathogens, and choose the last time point of Day 0 1 3 5 for each subject
## H1N1 data (83 samples from 21 subjects; H1N1_021 has one data point missing)
h1n1_meta <- all_meta %>% 
  filter(pathogen == "H1N1") %>% 
  filter(time_hr == 0 | time_hr == 21.5 | time_hr == 69.5 | time_hr == 108) 
h1n1_id <- h1n1_meta$subject_id
h1n1_bulk <- all_bulk[, h1n1_id]

## H3N2 data (67 samples from 17 subjects; H3N2 008 has one data point missing)
h3n2_meta <- all_meta %>% 
  filter(pathogen == "H3N2") %>% 
  filter(time_hr == 0 | time_hr == 21.5 | time_hr == 69.5 | time_hr == 108) 
h3n2_id <- h3n2_meta$subject_id
h3n2_bulk <- all_bulk[, h3n2_id]

## gene signatures
marker_id <- get_marker_genes(sig_mat, 160) # NOT include DC (800 genes)
sig_marker <- all_sig1[marker_id, ]
cor(sig_marker)

h1n1_bulk_marker <- h1n1_bulk[marker_id, ]
h3n2_bulk_marker <- h3n2_bulk[marker_id, ]

alg <-  "nnls"
intercept <- TRUE

# H1N1 --------------------------------------------------------------------

if (FALSE) {
  h1n1_R_hat <- get_R_hat(Y = h1n1_bulk_marker, Theta = sig_marker)
  h1n1_pca_R_hat <- prcomp(t(h1n1_R_hat), scale. = TRUE, center = TRUE)
  
  autoplot(h1n1_pca_R_hat, data = h1n1_meta, colour = "infected_ID") +
    labs(title = "PCA Plot on H1N1-Challenged Samples \nin Longitudinal Inoculation Study", subtitle = "LS Residual") +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_manual(name = "Infection Status", values = c("#8a508f", "#ffd380", "#2c4875")) 
  
  autoplot(h1n1_pca_R_hat, data = h1n1_meta, colour = "symptom") +
    labs(title = "PCA Plot on H1N1-Challenge Samples \nin Longitudinal Inoculation Study", subtitle = "LS Residual") +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_manual(name = "Symptom", values = c("#8a508f", "#ffd380", "#2c4875")) 
  
  autoplot(h1n1_pca_R_hat, data = h1n1_meta, colour = "day_fct") +
    labs(title = "PCA Plot on H1N1-Challenge Samples \nin Longitudinal Inoculation Study", subtitle = "LS Residual") +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_brewer(name = "Days Post Infection", palette = 3) 
  
  autoplot(h1n1_pca_R_hat, data = h1n1_meta, colour = "race") +
    labs(title = "PCA Plot on H1N1-Challenge Samples \nin Longitudinal Inoculation Study", subtitle = "LS Residual") +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_viridis_d(name = "Race")
  
  autoplot(h1n1_pca_R_hat, data = h1n1_meta, colour = "age_group") +
    labs(title = "PCA Plot on H1N1-Challenge Samples \nin Longitudinal Inoculation Study", subtitle = "LS Residual") +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_viridis_d(name = "Age Group")
}


h1n1_P_ls_nnls <- get_P_ls(h1n1_bulk_marker, sig_marker, "none", alg = alg, intercept = intercept) # q = 0
h1n1_P_ls_cutoff <- get_P_ls(h1n1_bulk_marker, sig_marker, "cutoff", alg = alg, intercept = intercept) # q = 1
h1n1_P_ls_be <- get_P_ls(h1n1_bulk_marker, sig_marker, "be", alg = alg, intercept = intercept) # q = 7 based on B = 49
h1n1_P_ls_tw <- get_P_ls(h1n1_bulk_marker, sig_marker, "tw", alg = alg, intercept = intercept) # q = 14 (obviously blew up)

cor(sig_marker, h1n1_P_ls_cutoff$Gamma_hat)

h1n1_P_nnls_long <- get_long_form(h1n1_P_ls_nnls$P_hat, h1n1_meta)
h1n1_P_nnls_long$Cell_type <- factor(h1n1_P_nnls_long$Cell_type, levels = cell_types)
h1n1_P_nnls_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr, y = Proportion, col = pathogen_symptom)) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), alpha = .6) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 12, base_family = "Times") +
  labs(x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_viridis_d(name = "Pathogen/Symptom") +
  guides(color = guide_legend(nrow = 2))

h1n1_P_nnls_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = day_fct, y = Proportion, fill = pathogen_symptom)) +
  geom_boxplot() +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 12, base_family = "Times") +
  labs(x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_viridis_d(name = "Pathogen/Symptom") +
  guides(fill = guide_legend(nrow = 2))


h1n1_P_cutoff_long <- get_long_form(h1n1_P_ls_cutoff$P_hat, h1n1_meta)
h1n1_P_cutoff_long$Cell_type <- factor(h1n1_P_cutoff_long$Cell_type, levels = cell_types)

h1n1_P_cutoff_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr, y = Proportion, col = pathogen_symptom)) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), alpha = .6) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 12, base_family = "Times") +
  labs(x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_viridis_d(name = "Pathogen/Symptom") +
  guides(color = guide_legend(nrow = 2))

h1n1_P_cutoff_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = day_fct, y = Proportion, fill = pathogen_symptom)) +
  geom_boxplot() +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 12, base_family = "Times") +
  labs(x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_viridis_d(name = "Pathogen/Symptom") +
  guides(fill = guide_legend(nrow = 2))

### Extra work
h1n1_P_ls_q2 <- get_P_ls(h1n1_bulk_marker, sig_marker, q = 2, alg = alg) # q = 2
cor(sig_marker, h1n1_P_ls_q2$Gamma_hat)

Gamma2 <- h1n1_P_ls_q2$Gamma_hat[, 2]
X <- cbind(Gamma2, 1, sig_marker)
rownames(X) <- rownames(sig_marker)
colnames(X)[2] <- "intercept"

B_hat <- apply(h1n1_bulk_marker, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(sig_marker))$x})
P_ext <- B_hat[-seq(2), ]
P_ext_hat <- apply(P_ext, 2, function(x) x/sum(x))
rownames(P_ext_hat) <- colnames(sig_marker)

h1n1_P_nnls_long2 <- get_long_form(P_ext_hat, h1n1_meta)
h1n1_P_nnls_long2$Cell_type <- factor(h1n1_P_nnls_long2$Cell_type, levels = cell_types)
h1n1_P_nnls_long2 %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr, y = Proportion, col = pathogen_symptom)) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), alpha = .6) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 12, base_family = "Times") +
  labs(x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_viridis_d(name = "Pathogen/Symptom") +
  guides(color = guide_legend(nrow = 2))

h1n1_P_nnls_long2 %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = day_fct, y = Proportion, fill = pathogen_symptom)) +
  geom_boxplot() +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 12, base_family = "Times") +
  labs(x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_viridis_d(name = "Pathogen/Symptom") +
  guides(fill = guide_legend(nrow = 2))


# H3N2 --------------------------------------------------------------------


