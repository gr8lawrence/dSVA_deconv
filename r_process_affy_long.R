## This file converts the probe IDs of the affymetrix file to gene symbols
## and subset bulk data and gene signatures to common genes
library(tidyverse)

## read the bulk data and gene signatures
bulk_df_probes <- read_tsv("../dSVA_datasets/GSE52428_bulk_raw.txt") #22,277 probes each
sig_df_probes <- read_tsv("../dSVA_datasets/GSE49910_GSE63038_sig_raw.txt") %>% 
  filter(!is.na(`NK cell`)) #21,866 probes
bulk_df_probes_filt <- right_join(bulk_df_probes, sig_df_probes %>% dplyr::select(ProbeID))
# all.equal(bulk_df_probes$ProbeID, sig_df_probes$ProbeID) # TRUE

colSums(bulk_df_probes_filt[, -1])
colSums(sig_df_probes[, -1])

## annotation file
probe_name_df <- read_csv("../dSVA_datasets/Affy_probe_name.csv") %>% 
  dplyr::rename(ProbeID = probe)
probe_convert_df <- tibble(ProbeID = bulk_df_probes_filt$ProbeID) %>% 
  inner_join(probe_name_df, by = "ProbeID") %>% 
  filter(!is.na(name))
## subset to the first gene
probe_convert_unique <- probe_convert_df %>% 
  group_by(ProbeID) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  dplyr::rename(Gene = name) #20,829 annotated probes


## convert probes to gene symbols
bulk_df_gene_filt <- probe_convert_unique %>% 
  left_join(bulk_df_probes_filt, by = "ProbeID") %>% 
  dplyr::select(-ProbeID)
sig_df_gene <- probe_convert_unique %>% 
  left_join(sig_df_probes, by = "ProbeID") %>% 
  dplyr::select(-ProbeID)

bulk_df_gene_filt %>% 
  group_by(Gene) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) # only 13,112 unique genes

## subset to unique genes by total expression
# df <- bulk_df_gene_filt
get_unique_gene <- function(df) {
  total_exp <- df[, -1] %>% 
    data.matrix() %>% 
    rowSums()
  df$total_exp <- total_exp
  df_unique_gene <- df %>% 
    group_by(Gene) %>% 
    arrange(desc(total_exp)) %>% 
    filter(row_number() == 1) %>% 
    ungroup()
  df_unique_gene
}

## the filtered data
bulk_df_gene_final <- get_unique_gene(bulk_df_gene_filt) %>% 
  dplyr::select(-total_exp)
sig_df_gene_final <- get_unique_gene(sig_df_gene) %>% 
  dplyr::select(-total_exp)

write_tsv(bulk_df_gene_final, "../dSVA_datasets/Longitudinal_bulk_final.txt")
write_tsv(sig_df_gene_final, "../dSVA_datasets/Longitudinal_sig_final.txt")




