## This file process the GSE73032 data from Yu
library(SummarizedExperiment)
library(tidyverse)

## Making a GEP from IRIS data
sig_df <- read_csv("../dSVA_datasets/IRIS_full_annotated.csv")
sig_mat <- data.matrix(sig_df[, -1])
rownames(sig_mat) <- sig_df$Gene

col_data <- read_csv("../dSVA_datasets/GSE22886_GPL96_coldata.csv")
cell_types <- unique(col_data$cell_type)

## CD8+ T cells
CD8T_id <- col_data %>% 
  filter(cell_type == cell_types[1]) %>% 
  dplyr::select(subject_id) %>% 
  unlist()
CD8T_GEP <- rowMeans(sig_mat[, CD8T_id])

## CD4+ T cells
CD4T_id <- col_data %>% 
  filter(cell_type == cell_types[2] & time == "control") %>% 
  dplyr::select(subject_id) %>% 
  unlist()
CD4T_GEP <- rowMeans(sig_mat[, CD4T_id])

## T cells
T_id <- c(CD8T_id, CD4T_id)
T_GEP <- rowMeans(sig_mat[, T_id])

## NK cells
NK_id <- col_data %>% 
  filter(cell_type == cell_types[4] & is.na(treatment_agent)) %>% 
  dplyr::select(subject_id) %>% 
  unlist()
NK_GEP <- rowMeans(sig_mat[, NK_id])

## B cells
B_id <- col_data %>% 
  filter(cell_type == cell_types[5]) %>% 
  dplyr::select(subject_id) %>% 
  unlist()
B_GEP <- rowMeans(sig_mat[, B_id])

## Plasma
Plasma_id <- col_data %>% 
  filter(cell_type == cell_types[8] & tissue == "peripheral blood") %>% 
  dplyr::select(subject_id) %>% 
  unlist()
Plasma_GEP <- rowMeans(sig_mat[, Plasma_id])

## Monocytes
Mono_id <- col_data %>% 
  filter(cell_type == cell_types[9]) %>% 
  dplyr::select(subject_id) %>% 
  unlist()
Mono_GEP <- rowMeans(sig_mat[, Mono_id])

## Neureophils
Neu_id <- col_data %>% 
  filter(cell_type == cell_types[11]) %>% 
  dplyr::select(subject_id) %>% 
  unlist()
Neu_GEP <- rowMeans(sig_mat[, Neu_id])

## final GEP
# sig_mat_common <- cbind(CD4T_GEP, CD8T_GEP, B_GEP, NK_GEP, Mono_GEP, Plasma_GEP, Neu_GEP)
# colnames(sig_mat_common) <- c("CD4_T", "CD8_T", "B", "NK", "Mono", "Plasma", "Neu")

sig_mat_common <- cbind(T_GEP, B_GEP, NK_GEP, Mono_GEP, Plasma_GEP, Neu_GEP)
colnames(sig_mat_common) <- c("T", "B", "NK", "Mono", "Plasma", "Neu")
saveRDS(sig_mat_common, "../dSVA_datasets/IRIS_sig_mat_all_genes.rds")


# Bulk matrix after batch effect removals ---------------------------------


## extract the matrices and column data and subset them with the IRIS data
## get the count matrices and the column data
bulk_raw <- readRDS("../dSVA_datasets/GSE73072_list_batchCorrect_unique.rds")
matrix_ls <- lapply(bulk_raw, assay) # list of matrices
coldata_ls <- lapply(bulk_raw, function(x) x @ colData) # list of column data
# prod(rownames(bulk_raw[[1]] @ colData) == colnames(matrix_ls[[1]])) ## should be 1

## get common genes
# get_common_genes <- function(bulk_mat, sig_mat) {
#   common_genes <- intersect(rownames(bulk_mat), rownames(sig_mat))
# }

# prod(rownames(matrix_ls[[1]]) == rownames(matrix_ls[[2]]))
# prod(rownames(matrix_ls[[1]]) == rownames(matrix_ls[[3]]))
# prod(rownames(matrix_ls[[1]]) == rownames(matrix_ls[[4]])) # should all be 1

common_genes <- intersect(rownames(matrix_ls[[1]]), rownames(sig_mat_common))

bulk_common_ls <- lapply(matrix_ls, function(x) x[common_genes, ])
sig_mat_common2 <- sig_mat_common[common_genes, ]

saveRDS(list(bulk_ls = bulk_common_ls,
             sig = sig_mat_common2,
             col_data = coldata_ls),
        "../dSVA_datasets/GSE73072_processed.rds")


# Bulk matrices before batch effect removals ------------------------------

bulk_raw_batch <- readRDS("../dSVA_datasets/GSE73072_list_unique.rds")
matrix_ls_batch <- lapply(bulk_raw_batch, assay) # list of matrices
coldata_ls_batch <- lapply(bulk_raw_batch, function(x) x @ colData) # list of column data
# prod(rownames(bulk_raw[[1]] @ colData) == colnames(matrix_ls[[1]])) ## should be 1

## get common genes
# get_common_genes <- function(bulk_mat, sig_mat) {
#   common_genes <- intersect(rownames(bulk_mat), rownames(sig_mat))
# }

# prod(rownames(matrix_ls_batch[[1]]) == rownames(matrix_ls_batch[[2]]))
# prod(rownames(matrix_ls_batch[[1]]) == rownames(matrix_ls_batch[[3]]))
# prod(rownames(matrix_ls_batch[[1]]) == rownames(matrix_ls_batch[[4]])) # should all be 1

common_genes_batch <- intersect(rownames(matrix_ls_batch[[1]]), rownames(sig_mat_common))

bulk_common_ls_batch <- lapply(matrix_ls_batch, function(x) x[common_genes_batch, ])
sig_mat_common2_batch <- sig_mat_common[common_genes_batch, ]

saveRDS(list(bulk_ls = bulk_common_ls_batch,
             sig = sig_mat_common2_batch,
             col_data = coldata_ls_batch),
        "../dSVA_datasets/GSE73072_processed_batch.rds")

