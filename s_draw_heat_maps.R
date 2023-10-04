## this script is for generating the heatmap
source("s_sources.R")
set.seed(100)

n <- 20
m <- 1000
K <- 5
q <- 1
err <- TRUE
p_sig <- c(0.25, 0.5, 0.75)
lambda <- 5
# gamma <- 3
gamma_seq <- c(1/4, 1, 4)
# n_sv <- 6 # the number of sv we want to plot
first_effect <- "bin"
second_effect <- "con"
B <- 100

if (q == 1) {
  effect_name <- case_when(first_effect == "bin" ~ "Binary Structure",
                           first_effect == "con" ~ "Continuous Confounding",
                           first_effect == "flat" ~ "Flat Effects",
                           first_effect == "cc" ~ "Case/control Structure",
                           first_effect == "miss" ~ "Missing Data Latent Structure")
} else if (q == 2) {
  effect_name <- case_when(first_effect == "bin" & second_effect == "bin" ~ "Binary + Binary Structure",
                           first_effect == "bin" & second_effect == "con" ~ "Binary Structure + Continuous Confounding",
                           first_effect == "con" & second_effect == "con" ~ "Continuous + Continuous Confounding",
                           first_effect == "cc" & second_effect == "bin" ~ "Case/control + Binary Structure",
                           first_effect == "cc" & second_effect == "con" ~ "Case/control + Continuous Structure")
}


ht_opt(heatmap_column_names_gp = gpar(fontface = "italic"), 
       heatmap_column_title_gp = gpar(fontsize = 10),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE
)

if (q == 1) {
  # pdf(paste0("plots/matrix_heatmap_", first_effect, ".pdf"))
  pdf(paste0("plots/matrix_heatmap_", first_effect, "_20230923.pdf"))
  
} else if (q == 2) {
  # pdf(paste0("plots/matrix_heatmap_q", q, "_", first_effect, "_", second_effect, ".pdf"))
  pdf(paste0("plots/matrix_heatmap_q", q, "_", first_effect, "_", second_effect, "_20230923.pdf"))
}
for (p in p_sig) {
  for (gamma in gamma_seq) {
    # p <- 1
    # gamma <- 4
    
    ## true data
    true_data <- dSVA_model_sim_intercept(m, n, K, q, p, lambda, gamma, err = err, first_effect = first_effect, second_effect = second_effect)
    
    ## draw heat maps of Y and Y_lat
    H_Y <- Heatmap(true_data$Y, column_title = "Bulk Expression", name = "value", cluster_rows = FALSE, cluster_columns = FALSE)
    H_Y_lat <- Heatmap(true_data$Y_lat, column_title = "Latent Structure", name = "value", cluster_rows = FALSE, cluster_columns = FALSE)
    
    ## put them into the same list
    H_both <- H_Y + H_Y_lat
    draw(H_both, 
         # column_title = paste0(effect_name, ", gamma=", gamma, ", p_sig=", p), 
         column_title = paste0(effect_name, " (q = ", q, ") example"), 
         column_title_gp = gpar(fontsize = 16))
  }
}
dev.off()





# H_X <- Heatmap(true_data$X, column_title = "Signature Matrix", name = "value", cluster_rows = FALSE, cluster_columns = FALSE)
# H_X
