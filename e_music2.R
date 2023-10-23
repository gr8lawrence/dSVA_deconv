library(MuSiC)
library(devtools)

## read the bulk data
fadista_bulk <- readRDS("../dSVA_datasets/bulk_Fadista.rds")
fadista_bulk <- fadista_bulk[, !is.na(fadista_bulk$hba1c)]
fadista_bulk$group <- ifelse(fadista_bulk$hba1c > 5.7, "T2D", "Healthy") # use 5.7 as in Fan et al.

## read the healthy sc data
cell_types <- c("acinar", "alpha", "beta", "delta", "ductal", "gamma")
healthy_sc <- readRDS("../dSVA_datasets/sc_healthy_Segerstolpe.rds")
healthy_sc <- healthy_sc[, healthy_sc$cellType %in% cell_types]

## make a simple signature matrix (averaging over cell types)
avg_exp <- matrix(NA, nrow(healthy_sc), length(cell_types))
rownames(avg_exp) <- rownames(healthy_sc)
colnames(avg_exp) <- cell_types
for (i in seq(1, length(cell_types))) {
  # i = 1
  exp_sub <- healthy_sc[, healthy_sc$cellType == cell_types[i]]
  avg_exp[, i] <- rowMeans(exprs(exp_sub))
}
avg_exp <- avg_exp[rowSums(avg_exp) > 10, ]

## does some feature selection

  