## This file is for analyzing the longitudinal data in GSE52428
library(tidyverse)
library(ggfortify)
library(ggrepel)
library(gprofiler2)
library(ggthemes)
source("a_dSVA_functions.R")
source("r_real_data_funs.R")
source("r_real_data_support_funs.R")


# Preprocessing -----------------------------------------------------------


## read bulk data, sig data, and metadata
bulk_df <- read_tsv("../dSVA_datasets/Longitudinal_bulk_final.txt") %>% 
  arrange(Gene)
sig_df <- read_tsv("../dSVA_datasets/Longitudinal_sig_final.txt") %>% 
  arrange(Gene)
meta_df <- read_csv("../dSVA_datasets/GSE52428_Metadata.csv") 

bulk_mat <- data.matrix(bulk_df[, -1]); rownames(bulk_mat) <- bulk_df$Gene
sig_mat <- data.matrix(sig_df[, -1]); rownames(sig_mat) <- sig_df$Gene
# all.equal(rownames(bulk_mat), rownames(sig_mat)) #TRUE

## clean meta_df
meta_clean <- meta_df %>% 
  dplyr::select(geo_accession, ID, pathogen, time_hr, symptom, symptom_score, 
                gender, age, race, infected_ID, infected_sample)
meta_clean <- meta_clean %>% 
  filter(infected_sample == "Y" | infected_sample == "N" | time_hr == 0) %>% 
  filter(ID != "H1N1_002" & ID != "H1N1_003" & ID != "H1N1_007")
include_id <- meta_clean$geo_accession # include only samples with clear inoculation designation
bulk_mat_inc <- bulk_mat[, include_id]
# all.equal(colnames(bulk_mat_inc), meta_clean$geo_accession) #TRUE

hr_levels <- meta_clean$time_hr %>% unique() %>% sort %>% as.character()
meta_clean$time_fct <- factor(meta_clean$time_hr, levels = hr_levels)
meta_clean$day_fct <- case_when(
  meta_clean$time_hr <= 0 ~ "Day 0",
  meta_clean$time_hr > 0 & meta_clean$time_hr <= 24 ~ "Day 1",
  meta_clean$time_hr > 24 & meta_clean$time_hr <= 48 ~ "Day 2",
  meta_clean$time_hr > 48 & meta_clean$time_hr <= 72 ~ "Day 3",
  meta_clean$time_hr > 72 & meta_clean$time_hr <= 96 ~ "Day 4",
  meta_clean$time_hr > 96 & meta_clean$time_hr <= 120 ~ "Day 5",
  meta_clean$time_hr > 120 & meta_clean$time_hr <= 144 ~ "Day 6",
  meta_clean$time_hr > 144 ~ "Day 7+",
)
meta_clean$day_fct <- as.factor(meta_clean$day_fct)

meta_clean$age_group <- case_when(
  meta_clean$age <= 22 ~ "19-22",
  meta_clean$age > 22 & meta_clean$age <= 25 ~ "23-25",
  meta_clean$age > 25 & meta_clean$age <= 28 ~ "26-28",
  meta_clean$age > 28 ~ "28+"
)
meta_clean$age_group <- as.factor(meta_clean$age_group)
meta_clean$pathogen_symptom <- paste(meta_clean$pathogen, meta_clean$symptom, sep = "/")
# table(meta_clean$infected_ID, meta_clean$symptom)

meta_clean$symptom_score_group <- case_when(
  meta_clean$symptom_score <= 5 ~ "0-5",
  meta_clean$symptom_score > 5 & meta_clean$symptom_score <= 15 ~ "6-15",
  meta_clean$symptom_score > 15 & meta_clean$symptom_score <= 25 ~ "16-25",
  meta_clean$symptom_score > 25 ~ "25+"
)
meta_clean$symptom_score_group <- factor(meta_clean$symptom_score_group, levels = c("0-5", "6-15", "16-25", "25+"))
meta_clean <- meta_clean %>% 
  dplyr::rename(subject_id = geo_accession)

table(meta_clean$ID, meta_clean$time_hr)

## save the data
if (FALSE) {
  all_data_ls <- list(
    bulk = bulk_mat_inc,
    sig = sig_mat,
    meta = meta_clean
  )
  saveRDS(all_data_ls, "../dSVA_datasets/Longitudinal_viral_challenge_data.rds")
}

# LS Residuals ------------------------------------------------------------

### Pathogens combined
## explore the residuals (if (FALSE) protects the code from being accidentally executed)
if (FALSE) {
  R_hat <- get_R_hat(Y = bulk_mat_inc, Theta = sig_mat)
  pca_R_hat <- prcomp(t(R_hat), scale. = TRUE, center = TRUE)
  p_infected_ID <- autoplot(pca_R_hat, data = meta_clean, colour = "infected_ID") +
    labs(title = "PCA Plot on Longitudinal Inoculation Study", subtitle = "LS Residual") +
    # ggrepel::geom_label_repel(aes(label = sample)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_manual(name = "Infection Status", values = c("#8a508f", "#ffd380", "#2c4875")) 
  p_infected_ID
  
  p_symptom <- autoplot(pca_R_hat, data = meta_clean, colour = "symptom") +
    labs(title = "PCA Plot on Longitudinal Inoculation Study", subtitle = "LS Residual") +
    # ggrepel::geom_label_repel(aes(label = sample)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_manual(name = "Symptom", values = c("#8a508f", "#ffd380", "#2c4875")) 
  p_symptom
  
  p_symptom_score_group <- autoplot(pca_R_hat, data = meta_clean, colour = "symptom_score_group") +
    labs(title = "PCA Plot on Longitudinal Inoculation Study", subtitle = "LS Residual") +
    # ggrepel::geom_label_repel(aes(label = sample)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_viridis_d(name = "Symptom Score (Group)") 
  p_symptom_score_group
  
  p_pathogen <- autoplot(pca_R_hat, data = meta_clean, colour = "pathogen") +
    labs(title = "PCA Plot on Longitudinal Inoculation Study", subtitle = "LS Residual") +
    # ggrepel::geom_label_repel(aes(label = sample)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_manual(name = "Pathogen", values = c("#8a508f", "#ffd380", "#2c4875")) 
  p_pathogen
  
  p_time_day <- autoplot(pca_R_hat, data = meta_clean, colour = "day_fct") +
    labs(title = "PCA Plot on Longitudinal Inoculation Study", subtitle = "LS Residual") +
    # ggrepel::geom_label_repel(aes(label = sample)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_viridis_d(name = "Time (Day)")
  p_time_day
  
  p_infected_sample <- autoplot(pca_R_hat, data = meta_clean, colour = "infected_sample") +
    labs(title = "PCA Plot on Longitudinal Inoculation Study", subtitle = "LS Residual") +
    # ggrepel::geom_label_repel(aes(label = sample)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_viridis_d(name = "Infected Sample")
  p_infected_sample
  
  p_race <- autoplot(pca_R_hat, data = meta_clean, colour = "race") +
    labs(title = "PCA Plot on Longitudinal Inoculation Study", subtitle = "LS Residual") +
    # ggrepel::geom_label_repel(aes(label = sample)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_viridis_d(name = "Race")
  p_race
  
  p_age_group <- autoplot(pca_R_hat, data = meta_clean, colour = "age_group") +
    labs(title = "PCA Plot on Longitudinal Inoculation Study", subtitle = "LS Residual") +
    # ggrepel::geom_label_repel(aes(label = sample)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_viridis_d(name = "Age Group")
  p_age_group
  
  p_path_symp <- autoplot(pca_R_hat, data = meta_clean, colour = "pathogen_symptom") +
    labs(title = "PCA Plot on Longitudinal Inoculation Study", subtitle = "LS Residual") +
    # ggrepel::geom_label_repel(aes(label = sample)) +
    theme_base(base_size = 12, base_family = "Times") +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background=element_blank()) +
    scale_color_viridis_d(name = "Pathogen/Symptom") +
    guides(color = guide_legend(nrow = 2))
  p_path_symp
}

# Deconvolution -----------------------------------------------------------

# sig_400 <- get_marker_genes(sig_mat, 100) # NOT include DC (800 genes)


P_ls_nnls <- get_P_ls(bulk_mat_inc, sig_mat, "none") # q = 0
P_ls <- get_P_ls(bulk_mat_inc, sig_mat, "cutoff") # q = 1
P_ls_be <- get_P_ls(bulk_mat_inc, sig_mat, "be") # q = 17 based on B = 49
P_ls_tw <- get_P_ls(bulk_mat_inc, sig_mat, "tw") # q = 115 (obviously blew up)

## save these results
if (FALSE) {
  saveRDS(P_ls_nnls, "../dSVA_datasets/GSE52428_results/NNLS_results.rds")
  saveRDS(P_ls, "../dSVA_datasets/GSE52428_results/dSVA_cutoff_results.rds")
  saveRDS(P_ls_be, "../dSVA_datasets/GSE52428_results/dSVA_be_results.rds")
  saveRDS(P_ls_tw, "../dSVA_datasets/GSE52428_results/dSVA_tw_results.rds")
}

P_nnls_long <- get_long_form(P_ls_nnls$P_hat, meta_clean)
p_nnls <- P_nnls_long %>% 
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
p_nnls

p_nnls_box <- P_nnls_long %>% 
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
p_nnls_box

P_cutoff_long <- get_long_form(P_ls$P_hat, meta_clean)
p_cutoff <- P_cutoff_long %>% 
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
p_cutoff

p_cutoff_box <- P_cutoff_long %>% 
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
p_cutoff_box
