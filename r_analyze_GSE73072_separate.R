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
all_data_ls <- readRDS("../dSVA_datasets/GSE73072_processed.rds")
all_bulk_ls <- all_data_ls$bulk_ls
virus_names <- names(all_data_ls$bulk_ls)
all_sig1 <- all_data_ls$sig[, -5]
# cell_types <- c("T cell", "B cell", "NK cell", "Monocyte", "Neutrophil")
cell_types <- colnames(all_sig1)
# all_sig1 <- all_sig <- all_sig[, cell_types]
# colnames(all_sig1) <- colnames(all_sig)
all_meta <- all_data_ls$col_data

### separate pathogens, and choose the last time point of Day 0 1 3 5 for each subject
## H1N1 data 
h1n1_meta_full <- all_meta$H1N1 %>% 
  as_tibble()
h1n1_meta <- h1n1_meta_full %>% 
  as_tibble() %>% 
  filter(time_hr == 0 | time_hr == 21 | time_hr == 22 | time_hr == 45 | time_hr == 46 | time_hr == 69 | time_hr == 70 | time_hr == 108) %>%
  filter(infected_ID != "Undecided") %>% # Now we want to clean the samples whose infected_ID is not undecided %>% 
  mutate(time_hr_char = as.character(time_hr),
         time_hr_char = fct_reorder(time_hr_char, time_hr)) %>% 
  mutate(time_day = case_when(
    time_hr == 0 ~ "Day 0",
    time_hr == 21 | time_hr == 22 ~ "Day 1",
    time_hr == 45 | time_hr == 46 ~ "Day 2",
    time_hr == 69 | time_hr == 70 ~ "Day 3",
    time_hr == 108 ~ "Day 5"
  )) %>% 
  arrange(time_hr)
h1n1_bulk <- all_bulk_ls$H1N1

# table(h1n1_meta_full$ID, h1n1_meta_full$time_hr)

## H3N2 data
h3n2_meta <- all_meta$H3N2 %>% 
  as_tibble() %>% 
  filter(time_hr == 0 | time_hr == 21 | time_hr == 22 | time_hr == 45 | time_hr == 46 | time_hr == 69 | time_hr == 70 | time_hr == 108) %>%
  filter(infected_ID != "Undecided") %>% # Now we want to clean the samples whose infected_ID is not undecided
  mutate(time_hr_char = as.character(time_hr),
         time_hr_char = fct_reorder(time_hr_char, time_hr)) %>% 
  mutate(time_day = case_when(
    time_hr == 0 ~ "Day 0",
    time_hr == 21 | time_hr == 22 ~ "Day 1",
    time_hr == 45 | time_hr == 46 ~ "Day 2",
    time_hr == 69 | time_hr == 70 ~ "Day 3",
    time_hr == 108 ~ "Day 5"
  )) %>% 
  arrange(time_hr) 
h3n2_bulk <- all_bulk_ls$H3N2

## HRV
hrv_meta <- all_meta$HRV %>% 
  as_tibble() %>% 
  filter(time_hr == 0 | time_hr == 24 | time_hr == 48 | time_hr == 72 | time_hr == 108) %>%
  filter(infected_ID != "Undecided") %>% # Now we want to clean the samples whose infected_ID is not undecided
  mutate(time_hr_char = as.character(time_hr),
         time_hr_char = fct_reorder(time_hr_char, time_hr)) %>% 
  mutate(time_day = case_when(
    time_hr == 0 ~ "Day 0",
    time_hr == 24 ~ "Day 1",
    time_hr == 48 ~ "Day 2",
    time_hr == 72 ~ "Day 3",
    time_hr == 108 ~ "Day 5"
  )) %>% 
  arrange(time_hr) 
hrv_bulk <- all_bulk_ls$HRV

## RSV
rsv_meta <- all_meta$RSV %>% 
  as_tibble() %>% 
  filter(time_hr == 0 | time_hr == 22 | time_hr == 46 | time_hr == 70 | time_hr == 108) %>%
  filter(infected_ID != "Undecided") %>% # Now we want to clean the samples whose infected_ID is not undecided
  mutate(time_hr_char = as.character(time_hr),
         time_hr_char = fct_reorder(time_hr_char, time_hr)) %>% 
  mutate(time_day = case_when(
    time_hr == 0 ~ "Day 0",
    time_hr == 22 ~ "Day 1",
    time_hr == 46 ~ "Day 2",
    time_hr == 70 ~ "Day 3",
    time_hr == 108 ~ "Day 5"
  )) %>% 
  arrange(time_hr) 
rsv_bulk <- all_bulk_ls$RSV

## gene signatures
marker_id <- get_marker_genes(all_sig1, 150) # NOT include DC (800 genes)
sig_marker <- all_sig1[marker_id, ]
cor(sig_marker)

## subset the data
h1n1_bulk_marker <- h1n1_bulk[marker_id, ]
h3n2_bulk_marker <- h3n2_bulk[marker_id, ]
hrv_bulk_marker <- hrv_bulk[marker_id, ]
rsv_bulk_marker <- rsv_bulk[marker_id, ]

## set method parameters
alg <-  "nnls"
intercept <- TRUE

## unique time points
timepoints <- c("Day 0", "Day 1", "Day 2", "Day 3", "Day 5")

# H1N1 --------------------------------------------------------------------

print("H1N1")
## Use four lists to save the results
h1n1_P_ls_nnls <- list()
h1n1_P_ls_cutoff <- list()
h1n1_P_ls_be <- list()
h1n1_P_ls_tw <- list()
plt_ls <- list()

## Analysis separated by time point
for (i in 1:length(timepoints)) {
  # i = 1
  timepoint <- timepoints[i]
  subject_id <- h1n1_meta$geo_accession[h1n1_meta$time_day == timepoint]
  Y_t <- h1n1_bulk_marker[, subject_id]
  Meta_t <- h1n1_meta[h1n1_meta$geo_accession %in% subject_id, ]
  # all.equal(Meta_t$geo_accession, colnames(Y_t)) # TRUE
  
  h1n1_R_hat <- get_R_hat(Y = Y_t, Theta = sig_marker)
  h1n1_pca_R_hat <- prcomp(t(h1n1_R_hat), scale. = TRUE, center = TRUE)
  plt_ls[[i]] <-  autoplot(h1n1_pca_R_hat, data = Meta_t, colour = "infected_ID") +
    labs(subtitle = paste("Timepoint =", timepoint)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_manual(name = "Infection Status", values = c("#ffd380", "#2c4875")) 
    # scale_color_manual(name = "Infection Status", values = c("#ffd380", "#8a508f", "#2c4875")) 
  
  h1n1_P_ls_nnls[[i]] <- get_P_ls(Y_t, sig_marker, "none", alg = alg, intercept = intercept) # q = 0
  h1n1_P_ls_cutoff[[i]] <- get_P_ls(Y_t, sig_marker, "cutoff", alg = alg, intercept = intercept) # q = 1
  h1n1_P_ls_be[[i]] <- get_P_ls(Y_t, sig_marker, "be", alg = alg, intercept = intercept) # q = 7 based on B = 49
  h1n1_P_ls_tw[[i]] <- get_P_ls(Y_t, sig_marker, "tw", alg = alg, intercept = intercept) # q = 14 (obviously blew up)
}

ggpubr::ggarrange(plotlist = plt_ls, nrow = 3, ncol = 2, common.legend = FALSE)
ggsave("GSE73072_plots/H1N1_Residual_PC1.pdf", width = 7.5, height = 9.5)

## check correlation
lapply(h1n1_P_ls_cutoff, function(x) cor(x$Gamma_hat, sig_marker))

# other two methods
# lapply(h1n1_P_ls_tw, function(x) cor(x$Gamma_hat, sig_marker))
lapply(h1n1_P_ls_be, function(x) cor(x$Gamma_hat, sig_marker))

## combine the results
P_nnls <- Reduce(cbind, lapply(h1n1_P_ls_nnls, function(x) x$P_hat))
P_cutoff <- Reduce(cbind, lapply(h1n1_P_ls_cutoff, function(x) x$P_hat))

## re-arrange columns
P_nnls <- P_nnls[, h1n1_meta$geo_accession]
P_cutoff <- P_cutoff[, h1n1_meta$geo_accession]

## plot the results
h1n1_P_nnls_long <- get_long_form2(P_nnls, h1n1_meta)
h1n1_P_nnls_long$Cell_type <- factor(h1n1_P_nnls_long$Cell_type, levels = cell_types)
line_plt_nnls <- h1n1_P_nnls_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr, y = Proportion, col = infected_ID)) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "H1N1: NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(color = guide_legend(nrow = 1))
ggsave("GSE73072_plots/H1N1_Props_Line_Plots_NNLS.pdf", line_plt_nnls, width = 9.5, height = 7.5)
ggsave("GSE73072_plots/H1N1_Props_Line_Plots_NNLS.tiff", line_plt_nnls, width = 19, height = 15)


boxplot_nnls <- h1n1_P_nnls_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr_char, y = Proportion, fill = infected_ID)) +
  geom_boxplot(alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "H1N1: NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("GSE73072_plots/H1N1_Props_Boxplots_NNLS.pdf", boxplot_nnls, width = 9.5, height = 7.5)
ggsave("GSE73072_plots/H1N1_Props_Boxplots_NNLS.tiff", boxplot_nnls, width = 19, height = 15)

h1n1_P_cutoff_long <- get_long_form2(P_cutoff, h1n1_meta)
h1n1_P_cutoff_long$Cell_type <- factor(h1n1_P_cutoff_long$Cell_type, levels = cell_types)
line_plt_cutoff <- h1n1_P_cutoff_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr, y = Proportion, col = infected_ID)) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "H1N1: dSVA + NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(color = guide_legend(nrow = 1))
ggsave("GSE73072_plots/H1N1_Props_Line_Plots_Cutoff.pdf", line_plt_cutoff, width = 9.5, height = 7.5)

boxplot_cutoff <- h1n1_P_cutoff_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr_char, y = Proportion, fill = infected_ID)) +
  geom_boxplot(alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "H1N1: dSVA + NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("GSE73072_plots/H1N1_Props_Boxplots_Cutoff.pdf", boxplot_cutoff, width = 9.5, height = 7.5)


# H3N2 --------------------------------------------------------------------

print("H3N2")
## Use four lists to save the results
h3n2_P_ls_nnls <- list()
h3n2_P_ls_cutoff <- list()
h3n2_P_ls_be <- list()
h3n2_P_ls_tw <- list()
plt_ls <- list()

## Analysis separated by time point
for (i in 1:length(timepoints)) {
  # i = 1
  timepoint <- timepoints[i]
  subject_id <- h3n2_meta$geo_accession[h3n2_meta$time_day == timepoint]
  Y_t <- h3n2_bulk_marker[, subject_id]
  Meta_t <- h3n2_meta[h3n2_meta$geo_accession %in% subject_id, ]
  # all.equal(Meta_t$geo_accession, colnames(Y_t)) # TRUE
  
  h3n2_R_hat <- get_R_hat(Y = Y_t, Theta = sig_marker)
  h3n2_pca_R_hat <- prcomp(t(h3n2_R_hat), scale. = TRUE, center = TRUE)
  plt_ls[[i]] <-  autoplot(h3n2_pca_R_hat, data = Meta_t, colour = "infected_ID") +
    labs(subtitle = paste("Timepoint =", timepoint)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_manual(name = "Infected", values = c("#ffd380", "#2c4875")) 
  # scale_color_manual(name = "Infection Status", values = c("#ffd380", "#8a508f", "#2c4875")) 
  
  h3n2_P_ls_nnls[[i]] <- get_P_ls(Y_t, sig_marker, "none", alg = alg, intercept = intercept) # q = 0
  h3n2_P_ls_cutoff[[i]] <- get_P_ls(Y_t, sig_marker, "cutoff", alg = alg, intercept = intercept) # q = 1
  h3n2_P_ls_be[[i]] <- get_P_ls(Y_t, sig_marker, "be", alg = alg, intercept = intercept) # q = 7 based on B = 49
  h3n2_P_ls_tw[[i]] <- get_P_ls(Y_t, sig_marker, "tw", alg = alg, intercept = intercept) # q = 14 (obviously blew up)
}

ggpubr::ggarrange(plotlist = plt_ls, nrow = 3, ncol = 2, common.legend = FALSE)
ggsave("GSE73072_plots/H3N2_Residual_PC1.pdf", width = 7.5, height = 9.5)

## check correlation
lapply(h3n2_P_ls_cutoff, function(x) cor(x$Gamma_hat, sig_marker))

# other two methods
# lapply(h3n2_P_ls_tw, function(x) cor(x$Gamma_hat, sig_marker)) # q_hat = 0 for certain timepoints
lapply(h3n2_P_ls_be, function(x) cor(x$Gamma_hat, sig_marker))

## combine the results
P_nnls <- Reduce(cbind, lapply(h3n2_P_ls_nnls, function(x) x$P_hat))
P_cutoff <- Reduce(cbind, lapply(h3n2_P_ls_cutoff, function(x) x$P_hat))

## re-arrange columns
P_nnls <- P_nnls[, h3n2_meta$geo_accession]
P_cutoff <- P_cutoff[, h3n2_meta$geo_accession]

## plot the results
h3n2_P_nnls_long <- get_long_form2(P_nnls, h3n2_meta)
h3n2_P_nnls_long$Cell_type <- factor(h3n2_P_nnls_long$Cell_type, levels = cell_types)
line_plt_nnls <- h3n2_P_nnls_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr, y = Proportion, col = infected_ID)) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "H3N2: NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(color = guide_legend(nrow = 1))
ggsave("GSE73072_plots/H3N2_Props_Line_Plots_NNLS.pdf", line_plt_nnls, width = 9.5, height = 7.5)
ggsave("GSE73072_plots/H3N2_Props_Line_Plots_NNLS.tiff", line_plt_nnls, width = 19, height = 15)


boxplot_nnls <- h3n2_P_nnls_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr_char, y = Proportion, fill = infected_ID)) +
  geom_boxplot(alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "H3N2: NNLS", x = "Time (day)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("GSE73072_plots/H3N2_Props_Boxplots_NNLS.pdf", boxplot_nnls, width = 9.5, height = 7.5)
ggsave("GSE73072_plots/H3N2_Props_Boxplots_NNLS.tiff", boxplot_nnls, width = 19, height = 15)

h3n2_P_cutoff_long <- get_long_form2(P_cutoff, h3n2_meta)
h3n2_P_cutoff_long$Cell_type <- factor(h3n2_P_cutoff_long$Cell_type, levels = cell_types)
line_plt_cutoff <- h3n2_P_cutoff_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr, y = Proportion, col = infected_ID)) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "H3N2: dSVA + NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(color = guide_legend(nrow = 1))
ggsave("GSE73072_plots/H3N2_Props_Line_Plots_Cutoff.pdf", line_plt_cutoff, width = 9.5, height = 7.5)

boxplot_cutoff <- h3n2_P_cutoff_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr_char, y = Proportion, fill = infected_ID)) +
  geom_boxplot(alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "H3N2: dSVA + NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("GSE73072_plots/H3N2_Props_Boxplots_Cutoff.pdf", boxplot_cutoff, width = 9.5, height = 7.5)

# HRV ---------------------------------------------------------------------

print("HRV")
## Use four lists to save the results
hrv_P_ls_nnls <- list()
hrv_P_ls_cutoff <- list()
hrv_P_ls_be <- list()
hrv_P_ls_tw <- list()
plt_ls <- list()

## Analysis separated by time point
for (i in 1:length(timepoints)) {
  # i = 2
  timepoint <- timepoints[i]
  subject_id <- hrv_meta$geo_accession[hrv_meta$time_day == timepoint]
  Y_t <- hrv_bulk_marker[, subject_id]
  Meta_t <- hrv_meta[hrv_meta$geo_accession %in% subject_id, ]
  # all.equal(Meta_t$geo_accession, colnames(Y_t)) # TRUE
  
  hrv_R_hat <- get_R_hat(Y = Y_t, Theta = sig_marker)
  hrv_pca_R_hat <- prcomp(t(hrv_R_hat), scale. = TRUE, center = TRUE)
  plt_ls[[i]] <-  autoplot(hrv_pca_R_hat, data = Meta_t, colour = "infected_ID") +
    labs(subtitle = paste("Timepoint =", timepoint)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    # scale_color_manual(name = "Infection Status", values = c("#ffd380", "#2c4875")) 
  scale_color_manual(name = "Infected", values = c("#ffd380", "#8a508f", "#2c4875"))
  
  hrv_P_ls_nnls[[i]] <- get_P_ls(Y_t, sig_marker, "none", alg = alg, intercept = intercept) # q = 0
  hrv_P_ls_cutoff[[i]] <- get_P_ls(Y_t, sig_marker, "cutoff", alg = alg, intercept = intercept) # q = 1
  hrv_P_ls_be[[i]] <- get_P_ls(Y_t, sig_marker, "be", alg = alg, intercept = intercept) # q = 7 based on B = 49
  hrv_P_ls_tw[[i]] <- get_P_ls(Y_t, sig_marker, "tw", alg = alg, intercept = intercept) # q = 14 (obviously blew up)
}

ggpubr::ggarrange(plotlist = plt_ls, nrow = 3, ncol = 2, common.legend = FALSE)
ggsave("GSE73072_plots/HRV_Residual_PC1.pdf", width = 7.5, height = 9.5)

## check correlation
lapply(hrv_P_ls_cutoff, function(x) cor(x$Gamma_hat, sig_marker))

# other two methods
# lapply(hrv_P_ls_tw, function(x) cor(x$Gamma_hat, sig_marker)) # q_hat = 0
lapply(hrv_P_ls_be, function(x) cor(x$Gamma_hat, sig_marker))

## combine the results
P_nnls <- Reduce(cbind, lapply(hrv_P_ls_nnls, function(x) x$P_hat))
P_cutoff <- Reduce(cbind, lapply(hrv_P_ls_cutoff, function(x) x$P_hat))

## re-arrange columns
P_nnls <- P_nnls[, hrv_meta$geo_accession]
P_cutoff <- P_cutoff[, hrv_meta$geo_accession]

## plot the results
hrv_P_nnls_long <- get_long_form2(P_nnls, hrv_meta)
hrv_P_nnls_long$Cell_type <- factor(hrv_P_nnls_long$Cell_type, levels = cell_types)
line_plt_nnls <- hrv_P_nnls_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr, y = Proportion, col = infected_ID)) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "HRV: NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Infected", values = c("#ffd380", "#8a508f", "#2c4875")) +
  guides(color = guide_legend(nrow = 1))
ggsave("GSE73072_plots/HRV_Props_Line_Plots_NNLS.pdf", line_plt_nnls, width = 9.5, height = 7.5)
ggsave("GSE73072_plots/HRV_Props_Line_Plots_NNLS.tiff", line_plt_nnls, width = 19, height = 15)


boxplot_nnls <- hrv_P_nnls_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr_char, y = Proportion, fill = infected_ID)) +
  geom_boxplot(alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "HRV: NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_manual(name = "Infected", values = c("#ffd380", "#8a508f", "#2c4875")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("GSE73072_plots/HRV_Props_Boxplots_NNLS.pdf", boxplot_nnls, width = 9.5, height = 7.5)
ggsave("GSE73072_plots/HRV_Props_Boxplots_NNLS.tiff", boxplot_nnls, width = 19, height = 15)


hrv_P_cutoff_long <- get_long_form2(P_cutoff, hrv_meta)
hrv_P_cutoff_long$Cell_type <- factor(hrv_P_cutoff_long$Cell_type, levels = cell_types)
line_plt_cutoff <- hrv_P_cutoff_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr, y = Proportion, col = infected_ID)) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "HRV: dSVA + NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Infected", values = c("#ffd380", "#8a508f", "#2c4875")) +
  guides(color = guide_legend(nrow = 1))
ggsave("GSE73072_plots/HRV_Props_Line_Plots_Cutoff.pdf", line_plt_cutoff, width = 9.5, height = 7.5)

boxplot_cutoff <- hrv_P_cutoff_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr_char, y = Proportion, fill = infected_ID)) +
  geom_boxplot(alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "HRV: dSVA + NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_manual(name = "Infected", values = c("#ffd380", "#8a508f", "#2c4875")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("GSE73072_plots/HRV_Props_Boxplots_Cutoff.pdf", boxplot_cutoff, width = 9.5, height = 7.5)

# RSV ---------------------------------------------------------------------

print("RSV")
## Use four lists to save the results
rsv_P_ls_nnls <- list()
rsv_P_ls_cutoff <- list()
rsv_P_ls_be <- list()
rsv_P_ls_tw <- list()
plt_ls <- list()

## Analysis separated by time point
for (i in 1:length(timepoints)) {
  # i = 1
  timepoint <- timepoints[i]
  subject_id <- rsv_meta$geo_accession[rsv_meta$time_day == timepoint]
  Y_t <- rsv_bulk_marker[, subject_id]
  Meta_t <- rsv_meta[rsv_meta$geo_accession %in% subject_id, ]
  # all.equal(Meta_t$geo_accession, colnames(Y_t)) # TRUE
  
  rsv_R_hat <- get_R_hat(Y = Y_t, Theta = sig_marker)
  rsv_pca_R_hat <- prcomp(t(rsv_R_hat), scale. = TRUE, center = TRUE)
  plt_ls[[i]] <-  autoplot(rsv_pca_R_hat, data = Meta_t, colour = "infected_ID") +
    labs(subtitle = paste("Timepoint =", timepoint)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_manual(name = "Infected", values = c("#ffd380", "#2c4875")) 
  # scale_color_manual(name = "Infection Status", values = c("#ffd380", "#8a508f", "#2c4875")) 
  
  rsv_P_ls_nnls[[i]] <- get_P_ls(Y_t, sig_marker, "none", alg = alg, intercept = intercept) # q = 0
  rsv_P_ls_cutoff[[i]] <- get_P_ls(Y_t, sig_marker, "cutoff", alg = alg, intercept = intercept) # q = 1
  rsv_P_ls_be[[i]] <- get_P_ls(Y_t, sig_marker, "be", alg = alg, intercept = intercept) # q = 7 based on B = 49
  rsv_P_ls_tw[[i]] <- get_P_ls(Y_t, sig_marker, "tw", alg = alg, intercept = intercept) # q = 14 (obviously blew up)
}

ggpubr::ggarrange(plotlist = plt_ls, nrow = 3, ncol = 2, common.legend = FALSE)
ggsave("GSE73072_plots/RSV_Residual_PC1.pdf", width = 7.5, height = 9.5)

## check correlation
lapply(rsv_P_ls_cutoff, function(x) cor(x$Gamma_hat, sig_marker))

# other two methods
# lapply(rsv_P_ls_tw, function(x) cor(x$Gamma_hat, sig_marker)) # q_hat = 0
lapply(rsv_P_ls_be, function(x) cor(x$Gamma_hat, sig_marker))

## combine the results
P_nnls <- Reduce(cbind, lapply(rsv_P_ls_nnls, function(x) x$P_hat))
P_cutoff <- Reduce(cbind, lapply(rsv_P_ls_cutoff, function(x) x$P_hat))

## re-arrange columns
P_nnls <- P_nnls[, rsv_meta$geo_accession]
P_cutoff <- P_cutoff[, rsv_meta$geo_accession]

## plot the results
rsv_P_nnls_long <- get_long_form2(P_nnls, rsv_meta)
rsv_P_nnls_long$Cell_type <- factor(rsv_P_nnls_long$Cell_type, levels = cell_types)
line_plt_nnls <- rsv_P_nnls_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr, y = Proportion, col = infected_ID)) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "RSV: NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(color = guide_legend(nrow = 1))
ggsave("GSE73072_plots/RSV_Props_Line_Plots_NNLS.pdf", line_plt_nnls, width = 9.5, height = 7.5)
ggsave("GSE73072_plots/RSV_Props_Line_Plots_NNLS.tiff", line_plt_nnls, width = 19, height = 15)


boxplot_nnls <- rsv_P_nnls_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr_char, y = Proportion, fill = infected_ID)) +
  geom_boxplot(alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "RSV: NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("GSE73072_plots/RSV_Props_Boxplots_NNLS.pdf", boxplot_nnls, width = 9.5, height = 7.5)
ggsave("GSE73072_plots/RSV_Props_Boxplots_NNLS.tiff", boxplot_nnls, width = 19, height = 15)

rsv_P_cutoff_long <- get_long_form2(P_cutoff, rsv_meta)
rsv_P_cutoff_long$Cell_type <- factor(rsv_P_cutoff_long$Cell_type, levels = cell_types)
line_plt_cutoff <- rsv_P_cutoff_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr, y = Proportion, col = infected_ID)) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "RSV: dSVA + NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(color = guide_legend(nrow = 1))
ggsave("GSE73072_plots/RSV_Props_Line_Plots_Cutoff.pdf", line_plt_cutoff, width = 9.5, height = 7.5)

boxplot_cutoff <- rsv_P_cutoff_long %>% 
  # group_by(subject_id, Cell_type) %>% 
  ggplot(aes(x = time_hr_char, y = Proportion, fill = infected_ID)) +
  geom_boxplot(alpha = .8) +
  facet_wrap(~Cell_type) +
  theme_base(base_size = 18, base_family = "Times") +
  labs(title = "RSV: dSVA + NNLS", x = "Time (hour)") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_fill_manual(name = "Infected", values = c("#ffd380", "#2c4875")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("GSE73072_plots/RSV_Props_Boxplots_Cutoff.pdf", boxplot_cutoff, width = 9.5, height = 7.5)

