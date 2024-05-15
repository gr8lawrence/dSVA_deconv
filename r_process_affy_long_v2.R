## This file converts the probe IDs of the affymetrix file to gene symbols
## and subset bulk data and gene signatures to common genes
library(tidyverse)

## helper function
get_unique_gene <- function(df) {
  total_exp <- df[, -1] %>% 
    data.matrix() %>% 
    rowSums()
  df$total_exp <- total_exp
  df_unique_gene <- df %>% 
    group_by(Gene) %>% 
    arrange(desc(total_exp)) %>% 
    filter(row_number() == 1) %>% 
    ungroup() %>% 
    dplyr::select(-total_exp)
  return(df_unique_gene)
}

## read the bulk data and gene signatures (133A 2.0 array)
bulk_df_probes <- read_tsv("../dSVA_datasets/GSE52428_bulk_raw.txt") #22,277 probes each

## annotation file
probe_name_df <- read_csv("../dSVA_datasets/Affy_probe_name.csv") %>% 
  dplyr::rename(ProbeID = probe) %>% 
  rename(Gene = name) %>% 
  dplyr::select(Gene, ProbeID) %>% 
  filter(!is.na(Gene))
# probe_name_df

bulk_df_gene_symbol <- probe_name_df %>% 
  inner_join(bulk_df_probes, by = "ProbeID") %>% 
  dplyr::select(-ProbeID) 

bulk_df_unique <- get_unique_gene(bulk_df_gene_symbol) #14,387 genes


## signature matrix (133A array)
sig_df_probes <- read_csv("../dSVA_datasets/GSE22886_GPL96.csv")  #22,283 probes
proba_name_133A_df <- read_tsv("../dSVA_datasets/GPL96_57554.txt", skip = 16) %>% 
  dplyr::select(`Gene Symbol`, ID) %>% 
  rename(ProbeID = ID,
         Gene = `Gene Symbol`) %>% 
  mutate(Gene = str_split_i(Gene, " /// ", 1)) 

sig_df_gene_symbol <- proba_name_133A_df %>% 
  left_join(sig_df_probes, by = "ProbeID") %>% 
  filter(!is.na(Gene)) %>% 
  dplyr::select(-ProbeID) 

# all.equal(proba_name_133A_df$ProbeID, sig_df_probes$ProbeID) # TRUE
# sig_df_gene_symbol %>%
#   group_by(Gene) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))

sig_df_unique <- get_unique_gene(sig_df_gene_symbol)
write_csv(sig_df_unique, "../dSVA_datasets/IRIS_full_annotated.csv")

# Common Label ------------------------------------------------------------
common_genes <- intersect(bulk_df_unique$Gene, sig_df_unique$Gene) #12,174 genes
bulk_mat <- data.matrix(bulk_df_unique[, -1])
rownames(bulk_mat) <- bulk_df_unique$Gene
bulk_mat_common <- bulk_mat[common_genes, ]

sig_mat_full <- data.matrix(sig_df_unique[, -1])
rownames(sig_mat_full) <- sig_df_unique$Gene
sig_mat_common_all <- sig_mat_full[common_genes, ]


# Create GEP --------------------------------------------------------------
col_data <- read_csv("../dSVA_datasets/GSE22886_GPL96_coldata.csv")
cell_types <- unique(col_data$cell_type)

## CD8+ T cells
CD8T_id <- col_data %>% 
  filter(cell_type == cell_types[1]) %>% 
  dplyr::select(subject_id) %>% 
  unlist()
CD8T_GEP <- rowMeans(sig_mat_common_all[, CD8T_id])

## CD4+ T cells
CD4T_id <- col_data %>% 
  filter(cell_type == cell_types[2] & time == "control") %>% 
  dplyr::select(subject_id) %>% 
  unlist()
CD4T_GEP <- rowMeans(sig_mat_common_all[, CD4T_id])

## T cells
T_id <- c(CD8T_id, CD4T_id)
T_GEP <- rowMeans(sig_mat_common_all[, T_id])

## NK cells
NK_id <- col_data %>% 
  filter(cell_type == cell_types[4] & is.na(treatment_agent)) %>% 
  dplyr::select(subject_id) %>% 
  unlist()
NK_GEP <- rowMeans(sig_mat_common_all[, NK_id])

## B cells
B_id <- col_data %>% 
  filter(cell_type == cell_types[5]) %>% 
  dplyr::select(subject_id) %>% 
  unlist()
B_GEP <- rowMeans(sig_mat_common_all[, B_id])

## Plasma
Plasma_id <- col_data %>% 
  filter(cell_type == cell_types[8] & tissue == "peripheral blood") %>% 
  dplyr::select(subject_id) %>% 
  unlist()
Plasma_GEP <- rowMeans(sig_mat_common_all[, Plasma_id])

## Monocytes
Mono_id <- col_data %>% 
  filter(cell_type == cell_types[9]) %>% 
  dplyr::select(subject_id) %>% 
  unlist()
Mono_GEP <- rowMeans(sig_mat_common_all[, Mono_id])

## Neureophils
Neu_id <- col_data %>% 
  filter(cell_type == cell_types[11]) %>% 
  dplyr::select(subject_id) %>% 
  unlist()
Neu_GEP <- rowMeans(sig_mat_common_all[, Neu_id])

## final GEP
# sig_mat_common <- cbind(CD4T_GEP, CD8T_GEP, B_GEP, NK_GEP, Mono_GEP, Plasma_GEP, Neu_GEP)
# colnames(sig_mat_common) <- c("CD4_T", "CD8_T", "B", "NK", "Mono", "Plasma", "Neu")

sig_mat_common <- cbind(T_GEP, B_GEP, NK_GEP, Mono_GEP, Plasma_GEP, Neu_GEP)
colnames(sig_mat_common) <- c("T", "B", "NK", "Mono", "Plasma", "Neu")


# Metadata ----------------------------------------------------------------

meta_df <- read_csv("../dSVA_datasets/GSE52428_Metadata.csv") 

## clean meta_df
meta_clean <- meta_df %>% 
  dplyr::select(geo_accession, ID, pathogen, time_hr, symptom, symptom_score, 
                gender, age, race, infected_ID, infected_sample)
meta_clean <- meta_clean %>% 
  filter(infected_sample == "Y" | infected_sample == "N" | time_hr == 0) %>% 
  filter(ID != "H1N1_002" & ID != "H1N1_003" & ID != "H1N1_007")
include_id <- meta_clean$geo_accession # include only samples with clear inoculation designation
bulk_mat_inc <- bulk_mat_common[, include_id]
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


# Combine Data And Save ---------------------------------------------------
IRIS_data_ls <- list(bulk = bulk_mat_inc, sig = sig_mat_common, meta = meta_clean)
saveRDS(IRIS_data_ls, "../dSVA_datasets/Longitudinal_viral_challenge_data_IRIS.rds")
