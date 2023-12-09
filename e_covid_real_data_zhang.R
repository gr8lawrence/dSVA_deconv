library(tidyverse)
library(ggpubr)
library(GEOquery)
library(ggfortify)
library(ggrepel)
library(gprofiler2)
library(sva)
source("a_dSVA_functions.R")
# source("f_customized_package_functions.R")

## pre-process the data and organize the column data
bulk_zhang <- read_delim("../dSVA_datasets/GSE179627_gene_reads_count.txt")
covid_series_zhang <- getGEO(GEO = "GSE179627", GSEMatrix = TRUE)
covid_sm <- covid_series_zhang$GSE179627_series_matrix.txt.gz
covid_bulk_pdata <- as_tibble(Biobase::pData(covid_sm))
covid_bulk_cols <- covid_bulk_pdata %>% 
  dplyr::select(geo_accession, title, matches(":ch1"))
colnames(covid_bulk_cols) <- c("acc_id" ,"subject_id", "age", "cell_type", "disease_state", 
                               "immune_response", "sex")
covid_bulk_cols$title <- str_split_i(covid_bulk_cols$title, ", ", 2)
covid_bulk_cols$subject_id <- str_replace(covid_bulk_cols$subject_id, "-", "_")
covid_bulk_cols$immune_response <- factor(covid_bulk_cols$immune_response,
                                          levels = c("uninfected", "asymptomatic", "symptomatic",
                                                     "recovering",
                                                     "re-detectable positive patients"))
covid_bulk_cols %>% 
  group_by(immune_response) %>% 
  summarise(n = n())
# covid_bulk_cols %>% 
#   group_by(disease_state) %>% 
#   summarise(n = n()) # they only coded "COVID-19 patient"

covid_raw_mat <- data.matrix(bulk_zhang[-(1:2), -(1:4)])
rownames(covid_raw_mat) <- bulk_zhang$Gene_symbol[-(1:2)]
covid_raw_mat <- covid_raw_mat[, covid_bulk_cols$subject_id] # reorder the columns of the bulk matrix
covid_mat_filtered <- covid_raw_mat[rowSums(covid_raw_mat) > 0, ]

saveRDS(covid_mat_filtered, "../dSVA_datasets/covid19_bulk_zhang.rds") # save the data
saveRDS(covid_bulk_cols, "../dSVA_datasets/covid19_col_zhang.rds")

## draw QC and PCA plots
## pair-wise heatmap
seqUtils::pw.cor.heatmap(covid_mat_filtered) +
  ggsci::scale_fill_gsea()

## PCA plots
pca <- prcomp(t(covid_mat_filtered), scale. = TRUE, center = TRUE)
samp_id <- str_split_i(colnames(covid_mat_filtered), "_", 1)
plot(pca$sdev[1:20]^2, main = "Elbow Plot \nEigenvalues of Raw Bulk Matrix (COVID-19 Study, Zhang et al.)", ylab = "Eigenvalue",
     log = "y")
autoplot(pca, data = covid_bulk_cols, colour = "immune_response") +
  labs(title = "PC of COVID-19 Bulk Samples", subtitle = "By Immune Response") +
  geom_label_repel(aes(label = subject_id)) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Immune response")

PCAloadings <- data.frame(variables = rownames(pca$rotation), pca$rotation) %>% 
  as_tibble()
PCAloadings_sub <- PCAloadings %>% 
  arrange(desc(abs(PC1))) %>% 
  filter(row_number() <= 20)
PCAloadings_sub

autoplot(pca, data = covid_bulk_cols, colour = "sex") +
  labs(title = "PC of COVID-19 Bulk Samples", subtitle = "By Sex") +
  geom_label_repel(aes(label = subject_id)) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Sex")


covid_mat_filtered[, c("Con_1", "Con_2", "RP_1", "RP_2", "Mo_1", "Re_1")]

## ComBat-seq corrected bulk matrix
batch <- ifelse(grepl("Con_", colnames(covid_raw_mat)) | grepl("RP_", colnames(covid_raw_mat)), 1, 2)
cov_ir <- as.numeric(covid_bulk_cols$immune_response) - 1
cov_sex <- as.numeric(as.factor(covid_bulk_cols$sex)) - 1
cov_mat <- cbind(cov_ir, cov_sex)

covid_adjusted_mat <- ComBat_seq(covid_raw_mat, batch = batch, group = NULL)
# covid_adjusted_mat <- ComBat_seq(covid_raw_mat, batch = batch, covar_mod = cov_mat)

covid_adjusted_filtered <- covid_adjusted_mat[rowSums(covid_adjusted_mat) > 0, ]

saveRDS(covid_adjusted_filtered, "../dSVA_datasets/covid19_bulk_zhang_ComBat.rds") # save the data

## visualize the correlations now
seqUtils::pw.cor.heatmap(covid_adjusted_filtered) +
  ggsci::scale_fill_gsea()

## PCA plots
pca <- prcomp(t(covid_adjusted_filtered), scale. = TRUE, center = TRUE)
samp_id <- str_split_i(colnames(covid_adjusted_filtered), "_", 1)
plot(pca$sdev[1:20]^2, main = "Elbow Plot \nEigenvalues of Raw Bulk Matrix (COVID-19 Study, Zhang et al.)", ylab = "Eigenvalue",
     log = "y")
autoplot(pca, data = covid_bulk_cols, colour = "immune_response") +
  labs(title = "PC of COVID-19 Bulk Samples (ComBat-seq Corrected)", subtitle = "By Immune Response") +
  geom_label_repel(aes(label = subject_id)) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Immune response")

autoplot(pca, data = covid_bulk_cols, colour = "sex") +
  labs(title = "PC of COVID-19 Bulk Samples (ComBat-seq Corrected)", subtitle = "By Sex") +
  geom_label_repel(aes(label = subject_id)) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Sex")

PCAloadings <- data.frame(variables = rownames(pca$rotation), pca$rotation) %>% 
  as_tibble()
PCAloadings_sub <- PCAloadings %>% 
  arrange(desc(abs(PC1))) %>% 
  filter(row_number() <= 20)
PCAloadings_sub
