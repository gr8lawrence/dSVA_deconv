## This file contains all the official simulation functions
## performance metrics
my_cor <- function(P1, P2) mean(diag(cor(P1, P2)))
my_ccc <- function(P1, P2) {
  ccc <- vector(mode = "double", length = ncol(P1))
  for (i in 1:ncol(P1)) ccc[i] <- DescTools::CCC(x = P1[, i], y = P2[, i])$rho.c[1, 1]
  mean(ccc)
}
my_mse <- function(P1, P2) mean((P1 - P2)^2)
my_mae <- function(P1, P2) mean(abs(P1 - P2))

## z-score normalize an array
my_normal <- function(v) {
  (v - mean(v))/sd(v)
} 

## the functions for the case-control scenario
disturb <- function(alpha) {
  K <- length(alpha)
  alpha2 <- vector(mode = "double", length = K)
  alpha2[seq(K - 1)] <- alpha[seq(K - 1)] + runif(K - 1, 0, mean(alpha)/(K - 1))
  alpha2[K] <- alpha[K] - runif(1, 0, mean(alpha))
  alpha2
}

softmax <- function(x) exp(x)/sum(exp(x))
row_center <- function(X) t(apply(X, 1, function(x) x - mean(x)))

## get Theta (from project 2)
## add a rotation matrix to record the rotation
get_signatures <- function(m, n, K, p_sig = 0.5, lambda = 3, chi_df = 200) {
  W <- extraDistr::rbern(m, p_sig)
  D_exp <- diag(rchisq(n = m, df = chi_df)) # mean expression value of each gene
  expr_list <- lapply(1:m, function(i) {
    expr <- W[i] * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - W[i]) * extraDistr::rdirichlet(1, rep(1, K))
    names(expr) <- as.character(1:K)
    sample(expr)
  })
  X <- D_exp %*% t(simplify2array(expr_list))
  rho <- t(sapply(expr_list, function(expr) as.integer(names(expr[order(expr, decreasing = TRUE)]))))
  list(X = X, W = W, rho = rho)
}

## simulate a second set of signatures using the marker genes from gene_signatures
get_second_signatures <- function(m, n, K, W, rho, lambda = 3, chi_df = 200) {
  X2 <- matrix(NA_real_, nrow = m, ncol = K)
  D_exp <- diag(rchisq(n = m, df = chi_df)) # mean expression value of each gene
  for (i in 1:m) {
    expr2 <- W[i] * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - W[i]) * extraDistr::rdirichlet(1, rep(1, K))
    X2[i, ] <- expr2[rho[i, ]]
  }
  X2 <- D_exp %*% X2
  list(X = X2, W = W)
}

## build a new model with an intercept
dSVA_model_sim_intercept <- function(m, n, K, q = 1, p_sig = 0.5, lambda = 3, gamma = 2, chi_df = 200, err = FALSE, 
                                     first_effect = c("bin", "con", "small", "cc"), second_effect = c("con", "bin"), 
                                     p_w2 = 0.5, p_pert = 0.5) {
  sig_ls <- get_signatures(m, n, K, p_sig, lambda, chi_df)
  X <- sig_ls$X
  W <- sig_ls$W
  P_star <- matrix(nrow = K, ncol = n)
  
  if (first_effect != "cc") {
    P_star <- t(extraDistr::rdirichlet(n = n, alpha = rep(1, K)))
    Y_reg <- X %*% P_star
  } else {
    Alpha <- Alpha2 <- t(extraDistr::rdirichlet(n, 2^seq(1, K))) # 20240322 changed to 2^1:K from 1:K
    Alpha2[, seq(n/2 + 1, n)] <- apply(Alpha2[, seq(n/2 + 1, n)], 2, disturb)
    P_star1 <- apply(Alpha, 2, softmax)
    P1 <- apply(Alpha2, 2, softmax)
    Dmat <- matrix(c(rep(0, floor(n/2)), rep(1, ceiling(n/2))), nrow = q, ncol = n)
    rho <- sig_ls$rho
    sig_ls2 <- get_second_signatures(m, n, K, W, rho, lambda, 1/2 * chi_df)
    X2 <- extraDistr::rbern(m, p_pert) * sig_ls2$X
    Y <- cbind(X %*% P_star1[, Dmat[1, ] == 0], (X + X2) %*% P1[, Dmat[1, ] == 1])
    
    ## 20240229: fix the wrong P_star for cc simulations
    P_star <- cbind(P_star1[, Dmat[1, ] == 0], P1[, Dmat[1, ] == 1])
  }
  
  Dmat <- if (first_effect != "cc" | q >= 2) matrix(0, nrow = q, ncol = n) else if (first_effect == "cc" | q == 1) {
    matrix(c(rep(0, floor(n/2)), rep(1, ceiling(n/2))), nrow = q, ncol = n)
  }
  if (first_effect == "bin") Dmat[1, ] <- c(rep(0, floor(n/2)), rep(1, ceiling(n/2)))
  if (first_effect == "con") Dmat[1, ] <- rchisq(n, df = chi_df/(5 * K))
  
  if (q >= 2) {
    if (second_effect == "con") {
      Dmat[2, ] <- rchisq(n, df = chi_df/K)
    } else if (second_effect == "bin") {
      Dmat[2, ] <- rep(c(rep(0, floor(n/4)), rep(1, ceiling(n/4))), 2)
    }
  }  
  
  if (first_effect == "small") {
    Y_lat <- get_Y_lat_small_effects(m, n, chi_df, gamma, W)
  } else if (first_effect == "bin") {
    Y_lat <- get_Y_lat_binary(m, n, chi_df/2, gamma, W)
  } else if (first_effect == "con") {
    Y_lat <- get_Y_lat_continuous(m, n, chi_df, gamma, W, Dmat[1, ], p_w2)
  } else if (first_effect == "cc") {
    Y_reg <- X %*% P_star 
    Y_lat <- X2 %*% (P1 - P_star) 
  } 
  
  if (q >= 2) {
    Y_lat2 <- matrix(0, nrow = m, ncol = n)
    if (second_effect == "con") {
      W2 <- extraDistr::rbern(m, 0.5)
      for (i in 1:n) {
        Y_lat2[, i] <- W2 * (W * rnorm(m, gamma * chi_df/4, sqrt(gamma * chi_df/3)) + (1 - W) * rnorm(m, chi_df/4, sqrt(chi_df/3))) * Dmat[2, i]
      }
    } else if (second_effect == "bin") {
      Y_lat2 <- get_Y_lat_binary2(m, n, chi_df/2, gamma, Dmat[2, ], W)
    }
    Y_lat <- Y_lat + Y_lat2
  }
  
  if (first_effect != "cc") {
    Y <- Y_reg + Y_lat
  } 
  
  if (err) {
    E <- matrix(0, nrow = m, ncol = n)
    for (i in 1:m) {
      E[i, ] <- rnorm(n, mean(Y_reg[i, ]), sqrt(mean(Y_reg[i, ])/5))
    }
    Y <- Y + E
  } else {
    E <- NA
  }
  
  p_neg <- sum(Y < 0)/(m * n)
  Y[Y < 0] <- 0
  
  ls <- list(X = X, P_star = P_star, D = Dmat, Y = Y, E = E, Y_lat = Y_lat, W = W, p_neg = p_neg)
  if (first_effect == "cc") {
    ls$X2 <- X2 
    ls$P_circ <- P_star1
    ls$P1 <- P1
  }
  return(ls)
}


## first binary effect
get_Y_lat_binary <- function(m, n, chi_df, gamma, W) {
  Y_lat <- matrix(0, m, n)
  
  ## for the group where Dmat = 0, select m/4 genes that increase in expression
  m1 <- ceiling(m/4)
  n1 <- floor(n/2)
  W1 <- W[seq(m1)]
  Y_lat[seq(m1), seq(n1)] <- W1 * rchisq(n = m1 * n1, df = gamma * chi_df) + (1 - W1) * rchisq(n = m1 * n1, df = chi_df)
  
  ## for the group where Dmat = 1, do the same
  m2 <- floor(m/4)
  n2 <- ceiling(n/2)
  W2 <- W[seq(m - m2 + 1, m)]
  Y_lat[seq(m - m2 + 1, m), seq(n - n2 + 1, n)] <- W2 * rchisq(n = m2 * n2, df = gamma * chi_df) + (1 - W2) * rchisq(n = m2 * n2, df = chi_df)
  
  ## return Y_lat
  Y_lat
}

## second binary effect
get_Y_lat_binary2 <- function(m, n, chi_df, gamma, d, W) {
  Y_lat <- matrix(0, m, n)
  
  ## for the group where d = 0, select m/4 genes that increase in expression
  m1 <- ceiling(m/4)
  n1 <- sum(d == 0)
  W1 <- W[seq(m1 + 1, ceiling(m/2))]
  Y_lat[seq(m1 + 1, ceiling(m/2)), (d == 0)] <- W1 * rchisq(n = m1 * n1, df = gamma * chi_df) + (1 - W1) * rchisq(n = m1 * n1, df = chi_df)
  
  ## for the group where d = 1, do the same
  m2 <- floor(m/4)
  n2 <- sum(d == 1)
  W2 <- W[seq(ceiling(m/2) + 1, ceiling(m/2) + m2)]
  Y_lat[seq(ceiling(m/2) + 1, ceiling(m/2) + m2), (d == 1)] <- W2 * rchisq(n = m2 * n2, df = gamma * chi_df) + (1 - W2) * rchisq(n = m2 * n2, df = chi_df)
  
  ## return Y_lat
  Y_lat
}

get_Y_lat_small_effects <- function(m, n, chi_df, gamma, W) {
  Y_lat <- matrix(0, m, n)
  
  ## for the group where Dmat = 0, select m/4 genes that increase in expression
  m1 <- ceiling(m/4)
  n1 <- floor(n/2)
  W1 <- W[seq(m1)]
  Y_lat[seq(m1), seq(n1)] <- W1 * rnorm(n = m1 * n1, mean = gamma * chi_df, sd = sqrt(chi_df/3)) + (1 - W1) * rnorm(n = m1 * n1, mean = chi_df, sd = sqrt(chi_df/3))
  
  ## for the group where Dmat = 1, do the same
  m2 <- floor(m/4)
  n2 <- ceiling(n/2)
  W2 <- W[seq(m - m2 + 1, m)]
  Y_lat[seq(m - m2 + 1, m), seq(n - n2 + 1, n)] <- W2 * rnorm(n = m2 * n2, mean = gamma * chi_df, sd = sqrt(chi_df/3)) + (1 - W2) * rnorm(n = m2 * n2, mean = chi_df, sd = sqrt(chi_df/3))
  
  ## return Y_lat
  Y_lat
}

## get continuous latent structures
get_Y_lat_continuous <- function(m, n, chi_df, gamma, W, d, p_w2 = 0.5) {
  Y_lat2 <- matrix(0, nrow = m, ncol = n)
  W2 <- extraDistr::rbern(m, p_w2)
  for (i in 1:n) Y_lat2[, i] <- W2 * (W * rnorm(m, gamma * chi_df/6, sqrt(gamma * chi_df/3)) + (1 - W) * rnorm(m, chi_df/6, sqrt(chi_df/3))) * d[i]
  Y_lat2
}

## get the singular values based on the sample-wise PCA
get_sing_vals <- function(true_data, intercept = TRUE, n_sv = 10) {
  ## see the distribution of top 10 eigenvalues for Y, Y_lat, R
  Y <- true_data$Y
  Y_lat <- true_data$Y_lat
  
  if (intercept) {
    X <- model.matrix(~1 + true_data$X)
  } else {
    X <- true_data$X
  }
  
  svd_X <- svd(X)
  U_x <- svd_X$u
  Sigma_x <- svd_X$d
  V_x <- svd_X$v
  B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y   
  R <- Y - X %*% B_star_hat
  
  d_Y <- prcomp(t(Y))$sdev^2
  d_Y_lat <- prcomp(t(Y_lat))$sdev^2
  d_R <- prcomp(t(R))$sdev^2
  
  mat <- rbind(d_Y[seq(n_sv)], d_Y_lat[seq(n_sv)], d_R[seq(n_sv)])
  colnames(mat) <- as.character(seq(n_sv))
  mat
}

## a function to produce bi-plots based on PCA results
get_bi_plot <- function(true_data, orient = c("gene", "sample"), intercept = TRUE) {
  ## extract Y and Y_lat
  Y <- true_data$Y
  Y_lat <- true_data$Y_lat

  ## calculate the OLS residuals
  if (intercept) {
    X <- model.matrix(~1 + true_data$X)
  } else {
    X <- true_data$X
  }
  svd_X <- svd(X)
  U_x <- svd_X$u
  Sigma_x <- svd_X$d
  V_x <- svd_X$v
  B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y
  R <- Y - X %*% B_star_hat
  
  ## name the matrices (or the biplot function would call an error)
  rownames(Y) <- rownames(Y_lat) <- rownames(R) <- paste("gene", as.character(1:m))
  colnames(Y) <- colnames(Y_lat) <- colnames(R) <- paste("sample", as.character(1:n)) 
  
  ## set the color keys
  ck <- c('0' = 'forestgreen', '1' ='gold')
  
  ## proceed with the orientation
  if (orient == "sample") {
    ## metadata (samples)
    md <- matrix(as.factor(true_data$D), n, 1)
    rownames(md) <- colnames(Y)
    colnames(md) <- "latent group"
    
    ## get the PCA objects (sample-wise)
    pca_Y <- pca(Y, metadata = md)
    pca_Y_lat <- pca(Y_lat, metadata = md)
    pca_R <- pca(R, metadata = md)
    
    ## save the plotted bi-plots
    rt_ls <- list(bi_Y = biplot(pca_Y, colby = "latent group") + labs(title = "Y", subtitle = "Sample-wise"),
                  bi_Y_lat = biplot(pca_Y_lat, colby = "latent group") + labs(title = "Latent Factors", subtitle = "Sample-wise"),
                  bi_R = biplot(pca_R, colby = "latent group") + labs(title = "Fitted Residuals", subtitle = "Sample-wise"))
    
  } else if (orient == "gene") {
    ## metadata (genes)
    md <- matrix(as.factor(true_data$W), m, 1)
    rownames(md) <- rownames(Y)
    colnames(md) <- "W"
    
    ## get the PCA objects (gene-wise)
    pca_Y <- pca(t(Y), metadata = md)
    pca_Y_lat <- pca(t(Y_lat), metadata = md)
    pca_R <- pca(t(R), metadata = md)
    
    ## save the plotted bi-plots
    rt_ls <- list(bi_Y = biplot(pca_Y, colby = "W", colkey = ck) + labs(title = "Y", subtitle = "Gene-wise"),
                  bi_Y_lat = biplot(pca_Y_lat, colby = "W", colkey = ck) + labs(title = "Latent Factors", subtitle = "Gene-wise"),
                  bi_R = biplot(pca_R, colby = "W", colkey = ck) + labs(title = "Fitted Residuals", subtitle = "Gene-wise"))
  }

  rt_ls
}

## 20240306: write a function that simulate the perturbation of the signature
dSVA_model_sim_pert <- function(m, n, K, q = 1, p_sig = 0.5, lambda = 3, chi_df = 200, err = FALSE, p_pert = 0.5) {
  # set.seed(100) # simulate a data for plotting

  sig_ls <- get_signatures(m, n, K, p_sig, lambda, chi_df)
  X <- sig_ls$X
  W <- sig_ls$W
  # P_star <- matrix(nrow = K, ncol = n)
  
  Alpha <- Alpha2 <- t(extraDistr::rdirichlet(n, 1:K))
  # Alpha2[, seq(n/2 + 1, n)] <- apply(Alpha2[, seq(n/2 + 1, n)], 2, disturb)
  P_star1 <- apply(Alpha, 2, softmax)
  P1 <- P_star1 
  # scales_P <- rlnorm(n, 0, 1)
  P1[nrow(P_star1), ] <- pmin(P_star1[nrow(P_star1), ] * rlnorm(n, 0, 1), .99)
  P1 <- apply(P1, 2, function(x) x/sum(x))
  # P1[seq(nrow(P1) - 1), ] <- P1[seq(nrow(P1) - 1), ] %*% diag(1 - P1[nrow(P1), ]) 
  Dmat <- matrix(c(rep(0, floor(n/2)), rep(1, ceiling(n/2))), nrow = q, ncol = n)
  rho <- sig_ls$rho
  # sig_ls2 <- get_second_signatures(m, n, K, W, rho, lambda, 1/2 * chi_df)
  v_change <- extraDistr::rbern(m, p_pert) 
  X_change <- (1 - v_change) * X[, ncol(X)] + v_change * X[, ncol(X)] * rlnorm(m, 0, 1)
  X2 <- cbind(X[, seq(ncol(X) - 1)],  X_change)
  Y <- cbind(X %*% P_star1[, Dmat[1, ] == 0], X2 %*% P1[, Dmat[1, ] == 1])
  
  ## 20240229: fix the wrong P_star for cc simulations
  P_star <- cbind(P_star1[, Dmat[1, ] == 0], P1[, Dmat[1, ] == 1])
  Y_lat <- (X2 - X) %*% P_star
  
  ## add errors
  if (err) {
    E <- matrix(0, nrow = m, ncol = n)
    for (i in 1:m) {
      E[i, ] <- rnorm(n, mean(Y[i, ]), sqrt(mean(Y[i, ])/5))
    }
    Y <- Y + E
  } else {
    E <- NA
  }
  
  p_neg <- sum(Y < 0)/(m * n)
  Y[Y < 0] <- 0
  
  ls <- list(X = X, P_star = P_star, D = Dmat, Y = Y, E = E, Y_lat = Y_lat, W = W, p_neg = p_neg)
  ls$X2 <- X2 
  ls$P_circ <- P_star1
  ls$P1 <- P1
  return(ls)
}

## plot the expression and proportion changes on the simulated data
## m = 2,000, n = 40, K = 10, p_sig = .25, lambda = 5, chi_df = 200, err = TRUE, p_pert = .75, seed = 100
# p_sig = .25
# p_pert = .75
# lambda = 5
# chi_df = 200
# err = TRUE
# 
# pdf("plots/Sim_COVAR_log_fc.pdf")
# hist(log(X_change/X[, ncol(X)]), breaks = 30,
#      main = "Histogram of Log FC of Simulated Expression \n Between Cases and Controls",
#      xlab = "Log fold change")
# dev.off()
# 
# prop_df <- tibble(Prop = c(P_star1[nrow(P_star1), Dmat[1, ] == 0], P1[nrow(P1), Dmat[1, ] == 1]),
#                   Disease = c(rep("Healthy Controls", sum(Dmat[1, ] == 0)), rep("Cases", sum(Dmat[1, ] == 1))))
# pdf("plots/Sim_COVAR_prop.pdf")
# boxplot(Prop ~ Disease, prop_df,
#         main = "Boxplot of Simulated Proportions",
#         ylab = "Proportion",
#         col = c("#8a508f", "#ffd380"))
# dev.off()
# 
# pdf("plots/Sim_COVAR_plots.pdf", width = 10.5, height = 5.5)
# par(mfrow = c(1, 2))
# hist(log(X_change/X[, ncol(X)]), breaks = 30,
#      main = "Histogram of Log FC of Simulated Expression \n Between Cases and Controls",
#      xlab = "Log fold change")
# boxplot(Prop ~ Disease, prop_df,
#         main = "Boxplot of Simulated Proportions",
#         ylab = "Proportion",
#         col = c("#8a508f", "#ffd380"))
# dev.off()

