## This file process the microarray data from GSE49910
library(GEOquery)
library(GEOmetadb)
library(tidyverse)

## download the bulk matrix (Affymetrix HG-U133_Plus_2): 54,675 probes
gse49910_tbl <- getGEO("GSE49910", GSEMatrix = TRUE)
gse49910_tbl2 <- gse49910_tbl$GSE49910_series_matrix.txt.gz
gse49910_bulk_mat <- Biobase::exprs(gse49910_tbl2)
saveRDS(gse49910_tbl2, "../dSVA_datasets/GSE49910_ExpressionSet.rds")
gse49910_bulk_df <- gse49910_bulk_mat %>% 
  as_tibble(rownames = "ProbeID")
write_tsv(gse49910_bulk_df, "../dSVA_datasets/GSE49910_bulk_raw.txt")

## process the matrix into a signature matrix
## build the cell-type metadata
gse49910_meta <- tibble(
  subject_id = gse49910_tbl2 @ phenoData @ data $geo_accession,
  cell_type = gse49910_tbl2 @ phenoData @ data $`cell type:ch1`,
  cell_type_subset = gse49910_tbl2 @ phenoData @ data $title
)
gse49910_meta <- gse49910_meta %>% 
  mutate(cell_type_correct = case_when(
    cell_type == "T cells" ~ "T cell",
    .default = cell_type
  ))
# all.equal(colnames(gse49910_bulk_mat), gse49910_meta$subject_id) # TRUE

## signature matrix 1: based on the umbrella cell_types
GEP_mat <- matrix(nrow = nrow(gse49910_bulk_mat))
CTs <- unique(gse49910_meta$cell_type_correct)
for (i in seq(length(CTs) - 1)) {
  # print(i)
  ct <- CTs[i]
  samp_id <- gse49910_meta$subject_id[which(gse49910_meta$cell_type_correct == ct)]
  bulk_sub <- gse49910_bulk_mat[, samp_id]
  GEP_mat <- cbind(GEP_mat, rowMeans(bulk_sub))
} 
GEP_mat <- GEP_mat[, -1]
colnames(GEP_mat) <- CTs[1:5]

## download the gene expression of NK cells (Affymetrix HG-U133A): 22,283 probes
gse63038_tbl <- getGEO("GSE63038", GSEMatrix = TRUE)
gse63038_tbl2 <- gse63038_tbl$GSE63038_series_matrix.txt.gz
gse63038_bulk_mat <- Biobase::exprs(gse63038_tbl2)
saveRDS(gse63038_tbl2, "../dSVA_datasets/GSE63038_ExpressionSet.rds")
gse63038_bulk_df <- gse63038_bulk_mat %>% 
  as_tibble(rownames = "ProbeID")
write_tsv(gse63038_bulk_df, "../dSVA_datasets/GSE49910_bulk_raw.txt")

## extract the metadata
gse63038_meta <- tibble(
  subject_id = gse63038_tbl2 @ phenoData @ data $geo_accession,
  character = gse63038_tbl2 @ phenoData @ data $source_name_ch1
)
# all.equal(colnames(gse63038_bulk_mat), gse63038_meta$subject_id) # TRUE

control_NK_id <- grepl("medium control", gse63038_meta$character)
control_NK_mat <- gse63038_bulk_mat[, control_NK_id]
NK_mean <- rowMeans(control_NK_mat)

## combine the signature matrix to get the 
common_probe_id <- intersect(rownames(GEP_mat), names(NK_mean)) ## 22,277 common probes
GEP_combined <- cbind(GEP_mat[common_probe_id, ], NK_mean[common_probe_id])
colnames(GEP_combined)[6] <- "NK cell"
colnames(GEP_combined)[2] <- "Neutrophil"
GEP_combined_final <- GEP_combined[, -5]

## save the final combined gene signatures
saveRDS(GEP_combined_final, "../dSVA_datasets/GEP_microarray_combined.rds")
GEP_combined_df <- GEP_combined_final %>% 
  as_tibble(rownames = "ProbeID")
write_tsv(GEP_combined_df, "../dSVA_datasets/GSE49910_GSE63038_sig_raw.txt")
