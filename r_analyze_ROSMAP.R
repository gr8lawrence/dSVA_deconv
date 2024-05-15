library(tidyverse)
library(ggpubr)
source("a_dSVA_functions.R")
source("r_real_data_funs.R")
library(edgeR)
library(limma)

## This data is not suitable as the elbow plot has no clear cutoffs

## read bulk data
bulk_RM <- read.delim("../dSVA_datasets/ROSMAP_data/geneExpr.txt") 
Y <- 2^bulk_RM

## read the single cell data
load("../dSVA_datasets/ROSMAP_data/GSE67835.RData") #Darmanis/DarmanisCells

## process the signature matrix
int <- intersect(rownames(bulk_RM),rownames(Darmanis))
y <- DGEList(counts= Darmanis[int,])
y <- calcNormFactors(y,method = 'TMM', Acutoff =  quantile(log2(Darmanis[,1]/sum(Darmanis[,1])),0.75))
Darmanis2 <- cpm(y, log = FALSE)
DarmanisMean <- sapply(unique(DarmanisCells),
                      function(x) rowMeans(Darmanis2[,names(which(DarmanisCells==x))]))

DarmanisSize <- sapply(unique(DarmanisCells),
                       function(x) rowMeans(Darmanis[,names(which(DarmanisCells==x))])) %>% colSums()

## process the bulk data for deconvolution
y <- cbind(DarmanisMean, Y[int, ])
y <- DGEList(counts= y)
y <- calcNormFactors(y,method = 'TMM', Acutoff =  quantile(log2(Darmanis[,1]/sum(Darmanis[,1])),0.75))
y <- cpm(y, log=FALSE)

## final data
Theta <- y$counts[, 1:5]
Y2 <- y$counts[, -(1:5)]

## load the real counts
astro_IHC <- read.delim("../dSVA_datasets/ROSMAP_data/IHC.astro.txt") 
endo_IHC <- read.delim("../dSVA_datasets/ROSMAP_data/IHC.endo.txt") 
microglia_IHC <- read.delim("../dSVA_datasets/ROSMAP_data/IHC.microglia.txt") 
neuro_IHC <- read.delim("../dSVA_datasets/ROSMAP_data/IHC.neuro.txt") 
oligo_IHC <- read.delim("../dSVA_datasets/ROSMAP_data/IHC.oligo.txt") 

all_samps <- intersect(names(oligo_IHC), 
                       intersect(intersect(names(astro_IHC), names(endo_IHC)), 
                                 intersect(names(microglia_IHC), names(neuro_IHC)))
)

P_IHC <- rbind(astro_IHC[, all_samps], 
               endo_IHC[, all_samps], 
               microglia_IHC[, all_samps], 
               neuro_IHC[, all_samps], 
               oligo_IHC[, all_samps])
rownames(P_IHC) <- colnames(Theta)
P_true <- apply(P_IHC, 2, function(x) x/sum(x))


P_nnls <- run_dSVA(Y2, Theta, intercept = FALSE, alg = "nnls", solver = "lsei", method = "none") 
P_nnls2 <- run_dSVA(Y2, Theta, intercept = TRUE, alg = "nnls", solver = "lsei", method = "cutoff")
P_nnls3 <- run_dSVA(Y2, Theta, intercept = TRUE, alg = "nnls", solver = "lsei", method = "be")
P_nnls4 <- run_dSVA(Y2, Theta, intercept = TRUE, alg = "nnls", solver = "lsei", method = "tw")


## compare the results
# P_est = P_nnls2[, all_samps]
B <- solve(t(Theta) %*% (Theta)) %*% t(Theta) %*% Y2 
R <- Y2 - Theta %*% B 
pca <- prcomp(t(R), scale. = TRUE, center = TRUE)
pca$sdev[1:100]

log(pca$sdev)[1:99] - log(pca$sdev)[2:100]
q = 2

P_nnls2 <- run_dSVA(Y2, Theta, intercept = FALSE, alg = "nnls", solver = "lsei", method = "cutoff")


comp_results <- function(P_est, P_true) {
  P_est_long <- P_est %>% 
    as_tibble(rownames = "cell_type") %>% 
    pivot_longer(-cell_type, names_to = "sample", values_to = "P_est")
  P_true_long <- P_true %>% 
    as_tibble(rownames = "cell_type") %>% 
    pivot_longer(-cell_type, names_to = "sample", values_to = "P_true")
  P_all <- left_join(P_est_long, P_true_long, by = c("sample", "cell_type"))
  P_all 
}

res_df <- comp_results(P_nnls[, all_samps], P_true)
res_df %>% 
  ggplot(aes(x = P_true, y = P_est, col = cell_type)) +
  geom_point() +
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2) +
  coord_fixed() +
  facet_wrap(~cell_type, nrow = 2)

res_df2 <- comp_results(P_nnls2[, all_samps], P_true)
res_df2 %>% 
  ggplot(aes(x = P_true, y = P_est, col = cell_type)) +
  geom_point() +
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2) +
  coord_fixed() +
  facet_wrap(~cell_type, nrow = 2)

my_mae(P_true, P_nnls2[, all_samps])
my_mae(P_true, P_nnls3[, all_samps])
my_mae(P_true, P_nnls4[, all_samps])
my_ccc(P_true, P_nnls[, all_samps])

my_ccc(P_true, P_nnls2[, all_samps])
my_mae(P_true, P_nnls[, all_samps])


# res_df3 <- comp_results(P_nnls3[, all_samps], P_true)
# res_df3 %>% 
#   ggplot(aes(x = P_true, y = P_est, col = cell_type)) +
#   geom_point() +
#   xlim(0, 1) +
#   ylim(0, 1) +
#   geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2) +
#   coord_fixed()
# 
# res_df4 <- comp_results(P_nnls4[, all_samps], P_true)
# res_df4 %>% 
#   ggplot(aes(x = P_true, y = P_est, col = cell_type)) +
#   geom_point() +
#   xlim(0, 1) +
#   ylim(0, 1) +
#   geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2) +
#   coord_fixed()
