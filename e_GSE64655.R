library(tidyverse)
HD30_ls <- readRDS('../dSVA_datasets/HD30_keep_TPM.rds')
Y30 <- HD30_ls$bulk
Ref30 <- HD30_ls$sig
HD31_ls <- readRDS('../dSVA_datasets/HD31_keep_TPM.rds')
Y31 <- HD31_ls$bulk
Ref31 <- HD31_ls$sig
Y <- cbind(Y30, Y31)
Ref_mean <- (Ref30 + Ref31)/2

# Ref_select <- Ref_mean[rowSums(Ref_mean) > 0, ]
# ct_names <- colnames(Ref_select)
# ct <- apply(Ref_select, 1, function(x) ct_names[which(x == max(x))])
# ref_genes <- tibble(gene = rownames(Ref_select)) %>% 
#   mutate(cell_type = ct)
# write_csv(ref_genes, "../dSVA_datasets/ref_genes.csv")

col_data <- tibble(
  subject = stringr::str_split_i(colnames(Y), "_", 1),
  day = stringr::str_split_i(colnames(Y), "_", 3)
)
X <- model.matrix(~1 + Ref_mean)
B <- solve(t(X) %*% (X)) %*% t(X) %*% Y  
R <- Y - X %*% B 
pca <- prcomp(t(R), scale. = TRUE, center = TRUE)
autoplot(pca, data = col_data, colour = "subject")
autoplot(pca, data = col_data, colour = "day")

s