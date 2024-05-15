library(tidyverse)

## function to find the top-n marker genes (by ratio-of-expression methods)
# get_marker_genes <- function(sig_sub1, n_marker) {
#   rat_df <- tibble(gene = rownames(sig_sub1),
#                    cell_type = apply(sig_sub1, 1, function(x) colnames(sig_sub1)[which(x == max(x))[1]]))
#   max_val <- apply(sig_sub1, 1, max)
#   sec_val <- apply(sig_sub1, 1, function(x) sort(x, decreasing = TRUE)[2])
#   rat_df <- rat_df %>% 
#     mutate(ratio = max_val/sec_val)
#   marker_genes <- rat_df %>% 
#     group_by(cell_type) %>% 
#     arrange(desc(ratio)) %>% 
#     filter(row_number() <= n_marker) %>% 
#     ungroup %>% 
#     dplyr::select(gene) %>% 
#     unlist
#   marker_genes
# }

get_marker_genes <- function(sig_sub1, n_marker, max_cutoff = 5) {
  rat_df <- tibble(gene = rownames(sig_sub1),
                   cell_type = apply(sig_sub1, 1, function(x) colnames(sig_sub1)[which(x == max(x))[1]]))
  max_val <- apply(sig_sub1, 1, max)
  sec_val <- apply(sig_sub1, 1, function(x) sort(x, decreasing = TRUE)[2])
  rat_df <- rat_df %>% 
    mutate(ratio = max_val/sec_val, max_val = max_val) %>% 
    filter(max_val > max_cutoff)
  marker_genes <- rat_df %>% 
    group_by(cell_type) %>% 
    arrange(desc(ratio)) %>% 
    filter(row_number() <= n_marker) %>% 
    ungroup %>% 
    dplyr::select(gene) %>% 
    unlist
  marker_genes
}

get_all_df <- function(controls, cases) {
  con_df_sub <- as_tibble(controls, rownames = "gene") %>% 
    pivot_longer(-1, names_to = "cell_type", values_to = "count") %>% 
    mutate(disease = "Healthy")
  case_df_sub <- as_tibble(cases, rownames = "gene") %>% 
    pivot_longer(-1, names_to = "cell_type", values_to = "count") %>% 
    mutate(disease = "Covid-19")
  all_df_sub <- rbind(con_df_sub, case_df_sub) %>% 
    mutate(log_count = log(count + 1, base = 10))
  all_df_sub
}

## read the bulk data
bulk_full <- readRDS("../dSVA_datasets/covid_bulk_count_processed.rds")

## gene signature from controls
sig_controls <- readRDS("../dSVA_datasets/Chuwen_sc_sig_sum_mid_500_cells_control.rds")
con_genes <- intersect(rownames(sig_controls), rownames(bulk_full))
sig_con2 <- sig_controls[con_genes, ]

## cases
sig_cases <- readRDS("../dSVA_datasets/Chuwen_sc_sig_sum_mid_500_cells.rds")
case_genes <- intersect(rownames(sig_cases), rownames(bulk_full))
sig_case2 <- sig_cases[case_genes, ]


# Plot the distribution of all genes' expression --------------------------
con_df <- as_tibble(sig_con2, rownames = "gene") %>% 
  pivot_longer(-1, names_to = "cell_type", values_to = "count") %>% 
  mutate(disease = "Healthy")
case_df <- as_tibble(sig_case2, rownames = "gene") %>% 
  pivot_longer(-1, names_to = "cell_type", values_to = "count") %>% 
  mutate(disease = "Covid-19")
all_df <- rbind(con_df, case_df) %>% 
  mutate(log_count = log(count + 1, base = 10))

all_df %>% 
  ggplot(aes(x = cell_type, y = log_count, fill = disease)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = "Log10(expression + 1)", x = "Cell type") +
  theme(legend.position = "bottom")

all_df %>% 
  ggplot(aes(x = log_count, fill = disease)) +
  geom_density(position = "identity", alpha = .6) +
  facet_wrap(~cell_type) +
  labs(x = "Log10(expression + 1)", y = "Count") +
  theme(legend.position = "bottom")

all_df %>% 
  ggplot(aes(x = log_count, fill = disease)) +
  geom_histogram(position = "identity", alpha = .6) +
  facet_wrap(~cell_type) +
  labs(x = "Log10(expression + 1)", y = "Count") +
  theme(legend.position = "bottom")


# Subset to marker genes --------------------------------------------------

## controls
sig_con_200 <- get_marker_genes(sig_con2[, 1:4], 200) # NOT include DC (800 genes)

## cases
sig_case_200 <- get_marker_genes(sig_case2[, 1:4], 200) #800 genes

## overlap
overlap <- intersect(sig_con_200, sig_case_200) #245 genes

### non-overlap

all_df_sub <- get_all_df(sig_con2[sig_con_200, 1:4], sig_case2[sig_case_200, 1:4]) 
sig_con_200_overlap <- intersect(sig_con_200, rownames(sig_case2))
all_df_sub_h <- get_all_df(sig_con2[sig_con_200_overlap, 1:4], sig_case2[sig_con_200_overlap, 1:4]) 
sig_case_200_overlap <- intersect(sig_case_200, rownames(sig_con2))
all_df_sub_c <- get_all_df(sig_con2[sig_case_200_overlap, 1:4], sig_case2[sig_case_200_overlap, 1:4]) 



all_df_sub %>% 
  ggplot(aes(x = cell_type, y = log_count, fill = disease)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "Markers Found Separately",
       y = "Log10(expression + 1)", 
       x = "Cell type") +
  theme(legend.position = "bottom")

all_df_sub_h %>% 
  ggplot(aes(x = cell_type, y = log_count, fill = disease)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "Markers Found on Healthy Controls",
       y = "Log10(expression + 1)", 
       x = "Cell type") +
  theme(legend.position = "bottom")

all_df_sub_c %>% 
  ggplot(aes(x = cell_type, y = log_count, fill = disease)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "Markers Found on COVID-19 Patients",
       y = "Log10(expression + 1)", 
       x = "Cell type") +
  theme(legend.position = "bottom")

all_df_sub %>% 
  ggplot(aes(x = log_count, fill = disease)) +
  geom_histogram(position = "identity", alpha = .6) +
  facet_wrap(~cell_type) +
  labs(title = "Markers Found Separately", x = "Log10(expression + 1)", y = "Count") +
  theme(legend.position = "bottom")

all_df_sub_h %>% 
  ggplot(aes(x = log_count, fill = disease)) +
  geom_histogram(position = "identity", alpha = .6) +
  facet_wrap(~cell_type) +
  labs(title = "Markers Found on Healthy Controls", x = "Log10(expression + 1)", y = "Count") +
  theme(legend.position = "bottom")

all_df_sub_c %>% 
  ggplot(aes(x = log_count, fill = disease)) +
  geom_histogram(position = "identity", alpha = .6) +
  facet_wrap(~cell_type) +
  labs(title = "Markers Found on COVID-19 Patients", x = "Log10(expression + 1)", y = "Count") +
  theme(legend.position = "bottom")



# Bulk data ---------------------------------------------------------------

bulk_cols <- readRDS("../dSVA_datasets/covid_col.rds")
bulk_cols$severity <- factor(bulk_cols$severity, levels = c("ICU",
                                                            "Severe",
                                                            "Moderate",
                                                            "Convalescent",
                                                            "Healthy"))
bulk_cols$disease_state_bin <- case_when(
  bulk_cols$disease_state == "Healthy" | bulk_cols$disease_state == "Convalescent" ~ "Healthy/Convalescent",
  .default = "COVID-19"
)
bulk_cols$sample <- str_split_i(bulk_cols$subject_id, "_", 1)

bulk_df <- as_tibble(bulk_full, rownames = "gene") %>% 
  pivot_longer(-1, names_to = "subject_id", values_to = "count") %>% 
  left_join(bulk_cols, by = "subject_id") %>% 
  mutate(log_count = log(count + 1, base = 10))

bulk_df %>% 
  ggplot(aes(x = sample, y = log_count, fill = disease_state_bin)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = "Log10(expression + 1)", x = "Subject") +
  theme(legend.position = "bottom")

bulk_df %>% 
  ggplot(aes(x = sample, y = log_count, fill = severity)) +
  geom_boxplot() +
  coord_flip() +
  labs(y = "Log10(expression + 1)", x = "Subject") +
  theme(legend.position = "bottom")
