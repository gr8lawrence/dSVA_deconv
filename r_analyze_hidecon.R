## This is the analysis of HiDecon data

## bulk data
load("../dSVA_datasets/HiDeConv/bulk_dat.rda") #bulk.dat

## ref data (each column is a purified bulk?)
load("../dSVA_datasets/HiDeConv/ref_dat.rda") # ref.dat
load("../dSVA_datasets/HiDeConv/ref_type.rda") # ref.type

load("../dSVA_datasets/HiDeConv/order_type.rda") # ref.type
load("../dSVA_datasets/HiDeConv/B.rda") # ref.type

pca <- prcomp(t(bulk.dat), scale. = TRUE, center = TRUE)
pca$sdev
plot(pca$sdev^2/sum(pca$sdev^2), main = "FHS Variation Explained")
