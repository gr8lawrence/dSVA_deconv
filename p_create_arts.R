## Draw the illustration art of the dSVA algorithm

source("s_sources.R")
library(ggfortify)
set.seed(100)

n <- 40
m <- 2000
K <- 10
q <- 2
err <- TRUE
p_sig <- 0.25
lambda <- 5
gamma <- 1
# n_sv <- 6 # the number of sv we want to plot
first_effect <- "bin"
second_effect <- "bin"
# B <- 100

## Generate the true data
true_data <- dSVA_model_sim_intercept(m, n, K, q, p_sig, lambda, gamma, err = err, first_effect = first_effect, second_effect = second_effect)

Heatmap(true_data$D)

pca <- prcomp(t(true_data$Y), scale. = TRUE, center = TRUE)

plot(pca$sdev[1:20]^2, main = "Elbow Plot", ylab = "Eigenvalue",
     log = "y")

## draw illustrations
ht_opt(heatmap_column_names_gp = gpar(fontface = "italic"), 
       heatmap_column_title_gp = gpar(fontsize = 10),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE
)

## Y
jpeg("images/illustrations/Y_strut.jpg")
p <- Heatmap(true_data$Y, column_title = NULL, name = "Expression", cluster_rows = FALSE, cluster_columns = FALSE)
draw(p)
dev.off()


## Theta
Heatmap(true_data$X, column_title = "Bulk Expression", name = "value", cluster_rows = FALSE, cluster_columns = FALSE)

## Proportions
Heatmap(true_data$P_star, column_title = "Bulk Expression", name = "value", cluster_rows = FALSE, cluster_columns = FALSE)

## Latent factor
jpeg("images/illustrations/D_strut.jpg", width = 1080, height = 240)
Heatmap(true_data$D, 
        column_title = NULL, 
        name = "value", 
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        rect_gp = gpar(col = "white", lwd = 1)
        )
dev.off()
