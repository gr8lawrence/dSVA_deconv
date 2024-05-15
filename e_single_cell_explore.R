
## TODO: this file explores whether marker genes in the Monaco matrix are differentially expressed

## read single cell files
library(data.table)
sc_DT <- fread("../dSVA_datasets/sc-dendritic/sc_signature_case-1024-lognormalized.txt")
sc_DT_controls <- fread("../dSVA_datasets/sc-dendritic/sc_signature_control-1024-lognormalized.txt")
# unique(colnames(sc_DT))

## read reference gene expression
# sig_full <- read.table("../dSVA_datasets/sigmatrixRNAseq.txt", header = TRUE, sep = "\t") %>% 
#   data.matrix()

covid_deconv_genes <- readRDS("../dSVA_datasets/covid_deconv_genes.rds")

## process the raw single cell matrices
sc_case_mat <- data.matrix(sc_DT[, -1])
rownames(sc_case_mat) <- sc_DT$genesymbol
colnames(sc_case_mat) <- colnames(sc_DT)[-1]

sc_control_mat <- data.matrix(sc_DT_controls[, -1])
rownames(sc_control_mat) <- sc_DT_controls$genesymbol
colnames(sc_control_mat) <- colnames(sc_DT_controls)[-1]

## subset the single cell files
gene_int_control <- intersect(covid_deconv_genes, rownames(sc_control_mat))
gene_int_case <- intersect(covid_deconv_genes, rownames(sc_case_mat))
gene_int_both <- intersect(gene_int_control, gene_int_case)

sc_control_sub <- sc_control_mat[gene_int_both, ] # 685 genes
sc_case_sub <- sc_case_mat[gene_int_both, ]

saveRDS(sc_control_sub, "../dSVA_datasets/covid_control_sub.rds")
saveRDS(sc_case_sub, "../dSVA_datasets/covid_case_sub.rds")

## find marker genes
cell_types <- unique(colnames(sc_DT))[-1]

case_sub_df <- as_tibble(t(sc_case_sub), rownames = "Cell_type") %>% 
  mutate(Cell_number = 1:ncol(sc_case_sub), Status = "Case") %>% 
  pivot_longer(-c(Cell_type, Cell_number, Status), names_to = "Gene", values_to = "Count")
case_sub_df

control_sub_df <- as_tibble(t(sc_control_sub), rownames = "Cell_type") %>% 
  mutate(Cell_number = 1:ncol(sc_control_sub), Status = "Control") %>% 
  pivot_longer(-c(Cell_type, Cell_number, Status), names_to = "Gene", values_to = "Count")

all_sub_df <- rbind(case_sub_df, control_sub_df)

## calculate mean fold change
mean_df <- all_sub_df %>% 
  group_by(Cell_type, Status, Gene) %>% 
  summarise(Mean_count = mean(Count)) %>% 
  ungroup %>% 
  pivot_wider(id_cols = c(Cell_type, Gene), names_from = Status, values_from = Mean_count) %>% 
  mutate(Log2_FC = log2(Case/Control))

ggplot(mean_df, aes(x = Log2_FC, fill = Cell_type)) +
  geom_histogram() +
  geom_vline(xintercept = 0, col = "navy", linetype = 2) +
  facet_wrap(. ~ Cell_type, scales = "free") +
  labs(title = "Cell-Type Expression Fold Changes",
       subtitle = "From Control to COVID-19 Samples",
       x = "Log2 fold change", y = "Count") +
  theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "none",
        plot.background = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5))
ggsave("plots/COVID_19/Expression_FC_20240303.pdf")

# Proportion Changes ------------------------------------------------------

case_true <- read_csv("../dSVA_datasets/case-celltype.csv") %>% 
  select(-1) %>% 
  dplyr::rename(cell_type = y,
                true_prop = prop) %>%
  left_join(meta_data, by = "sample") 
control_true <- read_csv("../dSVA_datasets/control-celltype.csv") %>%   
  select(-1) %>% 
  dplyr::rename(cell_type = y,
                true_prop = prop) %>% 
  left_join(meta_data, by = "sample")
all_true <- rbind(case_true, control_true)
all_true$cell_type <- factor(all_true$cell_type, levels = c("T cells", "B cells", "NK cells", "Monocytes"))

all_true %>% 
  ggplot(aes(x = cell_type, y = true_prop, fill = disease)) +
  geom_boxplot() +
  labs(title = "Cell-Type Proportion Changes",
       x = "Cell type", y = "Single-Cell proportion") +
  theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom",
        plot.background = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
        panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5)) +
  scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f"))
ggsave("plots/COVID_19/Cell_type_props_20240303.pdf")
