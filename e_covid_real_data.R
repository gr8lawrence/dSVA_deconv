library(GEOquery)
library(tidyverse)
library(biomaRt)
library(ggfortify)
library(ggrepel)
library(gprofiler2)
library(ggthemes)

## download the bulk and single cell data for the covid project

# Bulk data ---------------------------------------------------------------


covid_bulk_raw <- read_delim("../dSVA_datasets/Covid_Study1_RawCounts.txt") # bulk data
covid_series <- getGEO(GEO = "GSE152418", GSEMatrix = TRUE)
covid_sm <- covid_series$GSE152418_series_matrix.txt.gz
covid_bulk_pdata <- as_tibble(Biobase::pData(covid_sm))
covid_bulk_cols <- covid_bulk_pdata %>% 
  dplyr::select(geo_accession, title, matches(":ch1"))
colnames(covid_bulk_cols) <- c("acc_id" ,"subject_id", "cell_type", "post_infection", "disease_state", 
                               "gender", "city", "severity")
covid_raw_mat <- data.matrix(covid_bulk_raw[, -1])
colnames(covid_raw_mat) <- colnames(covid_bulk_raw)[-1]
rownames(covid_raw_mat) <- covid_bulk_raw$ENSEMBLID

## plot the boxplots
# seqUtils::expr.bp(covid_raw_mat, n = 34)
covid_mat_filtered <- covid_raw_mat[rowSums(covid_raw_mat) > 0, ]
covid_mat_filtered2 <- covid_mat_filtered
colnames(covid_mat_filtered2) <- str_split_i(colnames(covid_mat_filtered2), "_", 1)
seqUtils::pw.cor.heatmap(covid_mat_filtered2) +
  ggsci::scale_fill_gsea()

covid_bulk_cols$disease_state <- factor(covid_bulk_cols$disease_state, levels = c("COVID-19", "Healthy", "Convalescent"))
covid_bulk_cols$post_infection <- as.numeric(covid_bulk_cols$post_infection)
covid_bulk_cols$post_infection2 <- ifelse(covid_bulk_cols$post_infection > 20, 20, covid_bulk_cols$post_infection)
saveRDS(covid_bulk_cols, "../dSVA_datasets/covid_col.rds") 

covid_mat_filtered <- covid_mat_filtered[, covid_bulk_cols$subject_id]

## let's do a quick check (for batch effects and such)
pca <- prcomp(t(covid_mat_filtered), scale. = TRUE, center = TRUE)
samp_id <- str_split_i(colnames(covid_mat_filtered), "_", 1)
plot(pca$sdev[1:20]^2, main = "Elbow Plot \nEigenvalues of Raw Bulk Matrix (COVID-19 Study)", ylab = "Eigenvalue",
     log = "y")

## autoplot
autoplot(pca, data = covid_bulk_cols, colour = "disease_state") +
  labs(title = "PC of COVID-19 Bulk Samples", subtitle = "By Disease States") +
  geom_label_repel(aes(label = samp_id)) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Disease state")

# my_palette <- c("#ffd380", "#ff8531", "#ff6361", "#bc5090", "#8a508f", "#2c4875")
autoplot(pca, data = covid_bulk_cols, colour = "disease_state") +
  labs(title = "PCA Plot on COVID-19 Bulk Expression") +
  geom_label_repel(aes(label = samp_id)) +
  theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  ggsci::scale_color_tron(name = "Disease state")
 
autoplot(pca, data = covid_bulk_cols, colour = "disease_state") +
  labs(title = "By Disease Status") 

## try to plot a few loadings
# Extract PC axes for plotting
PCAvalues <- data.frame(sample = colnames(covid_mat_filtered), pca$x) %>% 
  as_tibble()

# Extract loadings of the variables
PCAloadings <- data.frame(variables = rownames(pca$rotation), pca$rotation) %>% 
  as_tibble()
PCAloadings_sub <- PCAloadings %>% 
  arrange(desc(abs(PC1))) %>% 
  filter(row_number() <= 20)
PCAloadings_sub

gene_symbol2 <- gconvert(PCAloadings_sub$variables,
                         organism = "hsapiens", 
                         target = "ENTREZGENE",
                         filter_na = FALSE)
PCAloadings_sub$gene_symbols <- gene_symbol2$target
PCAloadings_sub %>% 
  select(gene_symbols, PC1:PC34)

# Plot
# ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = sample)) +
#   geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*5),
#                                        yend = (PC2*5)), arrow = arrow(length = unit(1/2, "picas")),
#                color = "black") +
#   geom_point(size = 3) +
#   annotate("text", x = (PCAloadings_sub$PC1*5), y = (PCAloadings_sub$PC2*5),
#            label = PCAloadings_sub$variables)
# autoplot(pca, data = covid_bulk_cols, colour = "disease_state") +
#   labs(title = "PC of COVID-19 Bulk Samples", subtitle = "By Disease States") +
#   geom_label_repel(aes(label = samp_id)) +
#   ggpubr::theme_pubr() +
#   ggsci::scale_color_tron(name = "Disease state") + 
#   geom_segment(data = PCAloadings_sub, aes(x = 0, y = 0, xend = PC1 * 5,
#                                        yend = PC2 * 5), arrow = arrow(length = unit(1/2, "picas")), color = "black") +
#   annotate("text", x = (PCAloadings_sub$PC1*5), y = (PCAloadings_sub$PC2*5),
#          label = PCAloadings_sub$variables)

## check S175 - any gene dominates?
covid_mat_filtered2 %>% 
  as_tibble(rownames = "gene") %>% 
  arrange(desc(S175)) %>% 
  dplyr::select(gene, S175, S145, S147, S149, S150, S151, S155, S179)

covid_mat_filtered2 %>% 
  as_tibble(rownames = "gene") %>% 
  arrange(desc(S179)) %>% 
  dplyr::select(gene, S175, S145, S147, S149, S150, S151, S155, S179)

autoplot(pca, data = covid_bulk_cols, colour = "gender") +
  labs(title = "By Gender")

autoplot(pca, data = covid_bulk_cols, colour = "post_infection2") +
  labs(title = "By Days Post-Infection")

autoplot(pca, data = covid_bulk_cols, colour = "city") +
  labs(title = "By City")

autoplot(pca, data = covid_bulk_cols, colour = "severity") +
  labs(title = "By Disease Severity") +
  geom_label_repel(aes(label = samp_id)) +
  ggpubr::theme_pubr() +
  ggsci::scale_color_tron(name = "Disease state")

## translate ENSEMBL ID to gene symbols
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# G_list <- getBM(filters= "ensembl_gene_id", 
#                 attributes= c("ensembl_gene_id", "entrezgene", "description"), 
#                 values = gene_id, 
#                 mart = mart)
gene_id <- rownames(covid_mat_filtered)
gene_symbol <- gconvert(gene_id,
                        organism = "hsapiens", 
                        target = "ENTREZGENE",
                        filter_na = FALSE)
convert_tbl <- tibble(input = gene_symbol$input,
                      output = gene_symbol$target)
input_tbl <- table(gene_symbol$input)

## inspect the IDs mapped to multiple gene symbols
multi_list <- input_tbl[input_tbl > 1]
convert_tbl %>% 
  filter(input == "ENSG00000278233")
convert_tbl %>% 
  filter(output == "A2M") # checking the symbols are included

## find the IDs only mapped to one gene symbol
singular_list <- input_tbl[input_tbl == 1]
all_targets <- convert_tbl$output
names(all_targets) <- convert_tbl$input

## subset the bulk mat and use gene symbols as row names
covid_mat_genes <- covid_mat_filtered[names(singular_list), ]
# all.equal(rownames(covid_mat_genes), names(singular_targets))
rownames(covid_mat_genes) <- all_targets[rownames(covid_mat_genes)]
covid_bulk_symb <- covid_mat_genes[!is.na(rownames(covid_mat_genes)), ] # 24,280 genes

## one last step: make all gene symbols unique - choose the symbol with higher total expression
gene_tbl <- table(rownames(covid_bulk_symb))
good_genes <- names(gene_tbl[gene_tbl == 1])
multi_genes <- names(gene_tbl[gene_tbl > 1])
covid_bulk_good <- covid_bulk_symb[good_genes, ]
covid_bulk_rep <- covid_bulk_symb[which(rownames(covid_bulk_symb) %in% multi_genes), ] %>% 
  as_tibble(rownames = "gene")

rs <- rowSums(data.matrix(covid_bulk_rep[, -1]))
covid_bulk_rep <- covid_bulk_rep %>% 
  mutate(rs = rs) 

covid_bulk_fixed <- covid_bulk_rep %>% 
  group_by(gene) %>% 
  arrange(desc(rs)) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  dplyr::select(-rs) 
covid_bulk_fixed_mat <- data.matrix(covid_bulk_fixed[, -1])
rownames(covid_bulk_fixed_mat) <- covid_bulk_fixed$gene

covid_bulk_final <- rbind(covid_bulk_good, covid_bulk_fixed_mat) # 24,259 genes
# dim(covid_bulk_final)
saveRDS(covid_bulk_final, "../dSVA_datasets/covid_bulk_count_processed.rds")

# Single-cell data --------------------------------------------------------
## GSE64655
full_data_df <- read_csv("../dSVA_datasets/GSE64655_pc_raw_counts.csv") 

# %>% filter(!is.na(gene)) 

## let's inspect the data first
sorted_lines_df <- full_data_df %>% dplyr::select(-matches("PBMC"))
sorted_lines_mat <- data.matrix(sorted_lines_df[, -1])
rownames(sorted_lines_mat) <- c("Unmapped", sorted_lines_df$gene[-1])
sorted_lines_mat <- sorted_lines_mat[rowSums(sorted_lines_mat) > 0, ] # filter out genes with 0 expression

sorted_lines_col_data <- tibble(name = colnames(sorted_lines_df)[-1]) %>% 
  mutate(subject = str_split_i(name, "_", 1),
         cell_type = str_split_i(name, "_", 2),
         time = str_split_i(name, "_", 3)) %>% 
  mutate(log10_total_count = colSums(sorted_lines_mat))

pca_sc <- prcomp(t(sorted_lines_mat), scale. = TRUE, center = TRUE)
autoplot(pca_sc, data = sorted_lines_col_data, colour = "subject") +
  labs(title = "Sorted Cell Lines By Subject")

autoplot(pca_sc, data = sorted_lines_col_data, colour = "cell_type") +
  labs(title = "Sorted Cell Lines By Cell Type")

autoplot(pca_sc, data = sorted_lines_col_data, colour = "time") +
  labs(title = "Sorted Cell Lines By Time Point")

## summarize cell size (in log10)

sorted_lines_col_data %>% 
  ggplot(aes(x = cell_type, y = log10_total_count, fill = cell_type)) +
  geom_boxplot() +
  ggpubr::theme_pubr() +
  labs(title = "Total mRNA counts among cells")

## seems reasonable, now make signature matrix by averaging them
T_sig <- rowMeans(sorted_lines_mat[, grepl("_T_", colnames(sorted_lines_mat))])
B_sig <- rowMeans(sorted_lines_mat[, grepl("_B_", colnames(sorted_lines_mat))])
NK_sig <- rowMeans(sorted_lines_mat[, grepl("_NK_", colnames(sorted_lines_mat))])
Mono_sig <- rowMeans(sorted_lines_mat[, grepl("_Mono_", colnames(sorted_lines_mat))])
DC_sig <- rowMeans(sorted_lines_mat[, grepl("_DC_", colnames(sorted_lines_mat))])
Neut_sig <- rowMeans(sorted_lines_mat[, grepl("_Neut_", colnames(sorted_lines_mat))])

sig_total_mat <- cbind(T_sig, B_sig, NK_sig, Mono_sig, DC_sig)
colnames(sig_total_mat) <- c("T", "B", "NK", "Mono", "DC")
saveRDS(sig_total_mat, "../dSVA_datasets/GSE64655_sorted_count_processed.rds")

sig_total_mat <- cbind(T_sig, B_sig, NK_sig, Mono_sig, DC_sig, Neut_sig)
colnames(sig_total_mat) <- c("T", "B", "NK", "Mono", "DC", "Neut")
saveRDS(sig_total_mat, "../dSVA_datasets/GSE64655_sorted_count_whole_blood_processed.rds")

## only use the Day 0 data
T_sig <- rowMeans(sorted_lines_mat[, grepl("_T_0d", colnames(sorted_lines_mat))])
B_sig <- rowMeans(sorted_lines_mat[, grepl("_B_0d", colnames(sorted_lines_mat))])
NK_sig <- rowMeans(sorted_lines_mat[, grepl("_NK_0d", colnames(sorted_lines_mat))])
Mono_sig <- rowMeans(sorted_lines_mat[, grepl("_Mono_0d", colnames(sorted_lines_mat))])
DC_sig <- rowMeans(sorted_lines_mat[, grepl("_DC_0d", colnames(sorted_lines_mat))])
Neut_sig <- rowMeans(sorted_lines_mat[, grepl("_Neut_0d", colnames(sorted_lines_mat))])

sig_total_mat_0d <- cbind(T_sig, B_sig, NK_sig, Mono_sig, DC_sig)
colnames(sig_total_mat_0d) <- c("T", "B", "NK", "Mono", "DC")
saveRDS(sig_total_mat_0d, "../dSVA_datasets/GSE64655_sorted_count_0d_processed.rds")

sig_total_mat <- cbind(T_sig, B_sig, NK_sig, Mono_sig, DC_sig, Neut_sig)
colnames(sig_total_mat) <- c("T", "B", "NK", "Mono", "DC", "Neut")
saveRDS(sig_total_mat, "../dSVA_datasets/GSE64655_sorted_count_whole_blood_0d_processed.rds")
