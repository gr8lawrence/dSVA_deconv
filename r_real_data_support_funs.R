## function to find the top-n marker genes (by ratio-of-expression methods)
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

## function to get the estimated proportions
get_P_ls <- function(bulk_marker, sig_marker, method, alg = "nnls", q = NULL, intercept = TRUE) {
  if (!is.null(q)) { 
    P_ls <- run_dSVA(bulk_marker, 
                   sig_marker, 
                   intercept = intercept, 
                   exclude = NULL, 
                   alg = alg, 
                   solver = "lsei", 
                   q = q)
  } else {
    P_ls <- run_dSVA(bulk_marker, 
                     sig_marker, 
                     intercept = intercept, 
                     exclude = NULL, 
                     alg = alg, 
                     solver = "lsei", 
                     method = method)
  }
 
  P_ls
}

## get the prop boxplots
get_prop_boxplots <- function(P_long, subtitle) {
  P_long %>% 
    ggplot(aes(x = Cell_type, y = Proportion, fill = disease_state_bin)) +
    geom_boxplot() +
    labs(title = "Deconvoluted PBMC proportions", 
         subtitle = subtitle,
         x = "Cell type",
         caption = paste("Using ABIS-seq as reference; q =", P_ls_con_300$q_hat)) +
    theme_base(base_size = 12, base_family = "Times") +
    scale_fill_manual(name = "Disease State", values = c("#ffd380", "#8a508f")) +
    theme(plot.background = element_blank(),
          legend.position = "bottom",
          legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
          panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5)) 
  
}

## translate proportion matrix to a df
get_prop_df <- function(P_hat) {
  all_prop <- P_hat %>% 
    as_tibble(rownames = "cell_type") %>% 
    pivot_longer(-cell_type, 
                 names_to = "subject_id", 
                 values_to = "est_prop") %>% 
    mutate(subject_id = str_split_i(subject_id, "_", 1))
}

## get the comparison plot (the publication style)
get_compare_plot <- function(P_df, all_true, subtitle) {
  all_df <- all_true %>% 
    left_join(P_df , by = c("subject_id", "cell_type"))
  all_df$cell_type <- factor(all_df$cell_type, levels = c("T", "B", "NK", "Mono"))
  
  p <- ggplot(all_df, aes(x = true_prop, y = est_prop, col = cell_type, shape = cell_type)) +
    geom_point(size = 1.5) +
    # geom_text(aes(label = paste0("Pearson's correlation = ", pearson)), x = 0.6, y = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, col = "blue") + 
    scale_color_manual(name = "Cell type", values = c("#D72638", "#3F88C5", "#F49D37", "#140F2D")) +
    scale_shape_discrete(name = "Cell type") +
    coord_fixed() +
    xlim(0, 1) +
    ylim(0, 1) +
    facet_wrap(~disease) +
    theme_base(base_size = 12, base_family = "Times") +
    # theme_linedraw() +
    labs(title = "Estimated vs. Single Cell Proportions",
         subtitle = subtitle,
         x = "Single cell proportion", y = "Estimated proportion") +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "vertical",
          plot.background = element_blank(),
          legend.background = element_rect(fill = "white", color = "black", linewidth = .6),
          panel.grid = element_line(color = "lightgrey", linetype = 2, linewidth = .5))
  
  p
}

get_compare_stat <- function(P_df, all_true) {
  all_df <- all_true %>% 
    left_join(P_df , by = c("subject_id", "cell_type"))
  all_df$cell_type <- factor(all_df$cell_type, levels = c("T", "B", "NK", "Mono"))
  
  ## correlation by cell type - NNLS
  stat_by_CT <- all_df %>% 
    group_by(cell_type, disease) %>% 
    summarise(CCC = DescTools::CCC(true_prop, est_prop, na.rm = TRUE)$rho.c[1],
              MAE = mean(abs(true_prop - est_prop), na.rm = TRUE))
  
  ## correlation by sbuject - NNLS
  stat_by_prop <- all_df %>% 
    group_by(subject_id, disease) %>% 
    summarise(CCC = DescTools::CCC(true_prop, est_prop, na.rm = TRUE)$rho.c[1],
              # MSE = mean((true_prop - est_prop)^2, na.rm = TRUE),
              MAE = mean(abs(true_prop - est_prop), na.rm = TRUE)) 
  
  list(CT = stat_by_CT, subject = stat_by_prop)
}

## get the residuals
get_R_hat <- function(Y, Theta, intercept = TRUE) {
  
  ## add an intercept 
  if (intercept) {
    X <- model.matrix(~1 + Theta)
  } else {
    X <- Theta
  }
  
  ## perform first pass regression
  Bhat <- solve(t(X) %*% X) %*% t(X) %*% Y  
  R <- Y - X %*% Bhat
  
  R
}
