library(GEOquery)
library(GEOmetadb)
library(Biobase)
library(tidyverse)

### This file processes the gene expression for gene signatures from IRIS
## download the series matrix
gse22886_tbl <- getGEO("GSE22886", GSEMatrix = TRUE)
gse22886_GPL96_tbl <- gse22886_tbl$`GSE22886-GPL96_series_matrix.txt.gz`
gse22886_GPL97_tbl <- gse22886_tbl$`GSE22886-GPL97_series_matrix.txt.gz`

## download the expression data
gse22886_raw_mat <- getGEO("GSE22886", GSEMatrix = FALSE)

## get the platforms
gsmplatforms <- lapply(GSMList(gse22886_raw_mat), function(x) Meta(x)$platform_id) # should have both GPL96 and GPL97


# GPL96 -------------------------------------------------------------------

gsmlist <- Filter(function(gsm) Meta(gsm)$platform_id == 'GPL96', GSMList(gse22886_raw_mat))

probesets <- Table(GPLList(gse22886_raw_mat)[[1]])$ID #22,283 probes
data_matrix_GPL96 <- do.call('cbind', lapply(gsmlist, function(x) {
  tab <- Table(x)
  mymatch <- match(probesets, tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
data_matrix_GPL96 <- apply(data_matrix_GPL96, 2, function(x) as.numeric(as.character(x)))
rownames(data_matrix_GPL96) <- probesets

write_csv(data_matrix_GPL96 %>% as_tibble(rownames = "ProbeID"), "../dSVA_datasets/GSE22886_GPL96.csv")
# GPL97 -------------------------------------------------------------------

gsmlist2 <- Filter(function(gsm) Meta(gsm)$platform_id == 'GPL97', GSMList(gse22886_raw_mat))

probesets2 <- Table(GPLList(gse22886_raw_mat)[[2]])$ID #22,645 probes
data_matrix_GPL97 <- do.call('cbind', lapply(gsmlist2, function(x) {
  tab <- Table(x)
  mymatch <- match(probesets2, tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
data_matrix_GPL97 <- apply(data_matrix_GPL97, 2, function(x) as.numeric(as.character(x)))
rownames(data_matrix_GPL97) <- probesets2

write_csv(data_matrix_GPL97 %>% as_tibble(rownames = "ProbeID"), "../dSVA_datasets/GSE22886_GPL97.csv")

# Combine data ------------------------------------------------------------

## sort out the column data
col_data_GPL96 <- tibble(
  subject_id = gse22886_GPL96_tbl$geo_accession,
  cell_type = gse22886_GPL96_tbl$`cell type:ch1`,
  tissue = gse22886_GPL96_tbl$`tissue:ch1`,
  time = gse22886_GPL96_tbl$`time post-treatment:ch1`,
  treatment_agent = gse22886_GPL96_tbl$`treatment agent:ch1`
)
write_csv(col_data_GPL96, "../dSVA_datasets/GSE22886_GPL96_coldata.csv")

col_data_GPL97 <- tibble(
  subject_id = gse22886_GPL97_tbl$geo_accession,
  cell_type = gse22886_GPL97_tbl$`cell type:ch1`,
  tissue = gse22886_GPL97_tbl$`tissue:ch1`,
  time = gse22886_GPL97_tbl$`time post-treatment:ch1`,
  treatment_agent = gse22886_GPL97_tbl$`treatment agent:ch1`
)
write_csv(col_data_GPL97, "../dSVA_datasets/GSE22886_GPL97_coldata.csv")
## we can combine the two chips to get the HG133-Plus 2 probes
# both_probes <- intersect(rownames(data_matrix_GPL96), rownames(data_matrix_GPL97))
# cor(data_matrix_GPL96[both_probes, ], data_matrix_GPL97[both_probes, ])






