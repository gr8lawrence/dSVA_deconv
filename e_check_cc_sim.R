## check the simulated case-control study
## source files and libraries
source("s_sources.R")
set.seed(100)
true_data <- dSVA_model_sim_intercept(m = 2000, 
                                      n = 40, 
                                      K = 10, 
                                      q = 1, 
                                      p_sig = .1, 
                                      lambda = 5, 
                                      gamma = 1, 
                                      err = TRUE, 
                                      first_effect = "cc", 
                                      second_effect = "bin",
                                      p_pert = .1)

## column data
sim_col <- tibble(
                  cc = c(rep("Control", 20), rep("Case", 20)),
                  batch = rep(c(rep("Batch1", 10), rep("Batch2", 10)), 2))

## check the heatmap of the proportions
P_star <- true_data$P_star
rownames(P_star) <- paste0("CT_", seq(nrow(P_star)))
colnames(P_star) <- paste0("sample_", seq(ncol(P_star)))
P_star_df <- as_tibble(P_star, rownames = "Cell_type")
P_star

cor(true_data$P_star)

## check the PCA of the residuals
Y <- true_data$Y
Theta <- true_data$X

## add an intercept 
# if (intercept) {
#   X <- model.matrix(~1 + Theta)
# } else {
#   X <- Theta
# }

## perform first pass regression
X <- model.matrix(~1 + Theta)
Bhat <- solve(t(X) %*% X) %*% t(X) %*% Y  
R <- Y - X %*% Bhat

## plot pca
pca2 <- prcomp(t(R), scale. = TRUE, center = TRUE)
# pca$sdev
plot(pca2$sdev^2/sum(pca2$sdev^2), main = "Bulk Residuals PCA Variation Explained",
     sub = "Simulated Case-Control Structure (q = 1)",
     ylab = "Variance explained by PC", xlab = "PC")

p2 <- autoplot(pca2, data = sim_col, colour = "cc") +
  labs(title = "PCA Plot on SLE Bulk Residuals",
       caption = "Simulated Case-Control Structure (q = 1)") +
  # geom_label_repel(aes(label = sample)) +
  ggthemes::theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Disease status", values = c("#ffd380", "#8a508f")) 
p2


true_data2 <- dSVA_model_sim_intercept(m = 2000, 
                                      n = 40, 
                                      K = 10, 
                                      q = 2, 
                                      p_sig = .1, 
                                      lambda = 5, 
                                      gamma = 1, 
                                      err = TRUE, 
                                      first_effect = "cc", 
                                      second_effect = "bin",
                                      p_pert = .1)
Y <- true_data2$Y
Theta <- true_data2$X

## add an intercept 
# if (intercept) {
#   X <- model.matrix(~1 + Theta)
# } else {
#   X <- Theta
# }

## perform first pass regression
X <- model.matrix(~1 + Theta)
Bhat <- solve(t(X) %*% X) %*% t(X) %*% Y  
R <- Y - X %*% Bhat


## plot pca
pca3 <- prcomp(t(R), scale. = TRUE, center = TRUE)
# pca$sdev
plot(pca3$sdev^2/sum(pca3$sdev^2), main = "Bulk Residuals PCA Variation Explained",
     sub = "Simulated Case-Control Structure (q = 2)",
     ylab = "Variance explained by PC", xlab = "PC")

p3 <- autoplot(pca3, data = sim_col, colour = "cc") +
  labs(title = "PCA Plot on SLE Bulk Residuals",
       caption = "Simulated Case-Control Structure (q = 2)") +
  # geom_label_repel(aes(label = sample)) +
  ggthemes::theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Disease status", values = c("#ffd380", "#8a508f")) 
p3

p4 <- autoplot(pca3, data = sim_col, colour = "batch") +
  labs(title = "PCA Plot on SLE Bulk Residuals",
       caption = "Simulated Case-Control Structure (q = 2)") +
  # geom_label_repel(aes(label = sample)) +
  ggthemes::theme_base(base_size = 12, base_family = "Times") +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.background=element_blank()) +
  scale_color_manual(name = "Batch", values = c("#ffd380", "#8a508f")) 
p4

## Now let's try when only the proportion changes

