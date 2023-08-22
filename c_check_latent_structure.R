## this file is for checking the estimated latent structuress when dSVA performs badly
source("s_sources.R")
set.seed(100)
n <- 20
m <- 1000
K <- 5
q <- 1
err <- TRUE
p <- 0.75
lambda <- 5
gamma <- 4
n_sv <- 6 # the number of sv we want to plot
first_effect <- "con"
second_effect <- "con"

## simulate a truth
true_data <- dSVA_model_sim_intercept(m, n, K, q, p, lambda, gamma, err = err, first_effect = first_effect, second_effect = second_effect)

## check dSVA latent factors
sim_ls <- dsva_for_sim(Y = true_data$Y, Theta = true_data$X, n_comp = 1, intercept = TRUE, test = TRUE)
sim_ls$Y_lat_hat - true_data$Y_lat
my_mae(sim_ls$P_hat, true_data$P_star)
