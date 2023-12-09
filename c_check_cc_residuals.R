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
gamma_seq <- c(1/4, 1/2, 1, 2, 4)
n_sv <- 6 # the number of sv we want to plot
first_effect <- "cc"
B <- 10
plot_sv <- FALSE # whether to plot the singular/eigenvalues
true_data <- dSVA_model_sim_intercept(m, n, K, q, p, lambda, gamma, err = err, first_effect = first_effect, second_effect = "bin")
col_data <- tibble(subject = paste("sub", 1:ncol(Y), sep = "_"),
                   group = as.character(true_data$D))
Y <- true_data$Y
colnames(Y) <- col_data$subject
X <- model.matrix(~ 1 + true_data$X)
B <- solve(t(X) %*% (X)) %*% t(X) %*% Y  
R <- Y - X %*% B 
pca <- prcomp(t(R), scale. = TRUE, center = TRUE)
autoplot(pca, data = col_data, colour = "group")

