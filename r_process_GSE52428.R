## This file process the microarray data from GSE 52428
library(GEOmetadb)
library(tidyverse)

##  download the bulk matrix
gse52428_tbl <-  getGEO("GSE52428", GSEMatrix = TRUE)
gse52428_tbl2 <- gse52428_tbl$GSE52428_series_matrix.txt.gz
gse52428_bulk_mat <- Biobase::exprs(gse52428_tbl2)
saveRDS(gse52428_tbl2, "../dSVA_datasets/GSE52428_ExpressionSet.rds")
gse52428_bulk_df <- gse52428_bulk_mat %>% 
  as_tibble(rownames = "ProbeID")
write_tsv(gse52428_bulk_df, "../dSVA_datasets/GSE52428_bulk_raw.txt")

