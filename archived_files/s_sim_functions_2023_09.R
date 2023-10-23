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
get_signatures <- function(m, n, K, p_sig = 0.5, lambda = 3, chi_df = 200) {
  X <- matrix(NA_real_, nrow = m, ncol = K)
  W <- extraDistr::rbern(m, p_sig)
  D_exp <- diag(rchisq(n = m, df = chi_df)) # mean expression value of each gene
  for (i in 1:m) X[i, ] <- sample(W[i] * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - W[i]) * extraDistr::rdirichlet(1, rep(1, K)))
  X <- D_exp %*% X
  list(X = X, W = W)
}

## simulate a second set of signatures using the marker genes from gene_signatures
get_second_signatures <- function(m, n, K, W, lambda = 3, chi_df = 200) {
  X <- matrix(NA_real_, nrow = m, ncol = K)
  D_exp <- diag(rchisq(n = m, df = chi_df)) # mean expression value of each gene
  for (i in 1:m) X[i, ] <- sample(W[i] * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - W[i]) * extraDistr::rdirichlet(1, rep(1, K)))
  X <- D_exp %*% X
  list(X = X, W = W)
}

## TODO: this is generating signature matrix with measurement errors
get_signatures_bx <- function(m, n, K, p_sig = 0.5, lambda = 3, chi_df = 200, d = 0.15) {
  X <- matrix(NA_real_, nrow = m, ncol = K)
  W <- extraDistr::rbern(m, p_sig)
  D_exp <- diag(rchisq(n = m, df = chi_df)) # mean expression value of each gene
  for (i in 1:m) X[i, ] <- sample(W[i] * extraDistr::rdirichlet(1, c(lambda, rep(1, K - 1))) + (1 - W[i]) * extraDistr::rdirichlet(1, rep(1, K)))
  X <- D_exp %*% X
  X1 <- X
  X1[X1 < quantile(X1, d)] <- 0 # undetected reads (drop-outs)
  list(X1 = X1, X = X, W = W)
}

## TODO: simulate when q = 2, 3, 4 (1 discrete + 1 continuous, 1 discrete + 2 continuous, 2 discrete + 2 continuous)
## TODO: simulate continuous functions only when q = 1.
## allow errors (err = TRUE)
## build a new model with an intercept
dSVA_model_sim_intercept <- function(m, n, K, q = 1:4, p_sig = 0.5, lambda = 3, gamma = 2, chi_df = 200, err = FALSE, 
                                     first_effect = c("bin", "con", "small", "flat", "cc", "miss"), second_effect = c("con", "bin"), 
                                     p_w2 = 0.5, p_pert = 0.5) {
  
  ## generate signature matrices
  sig_ls <- get_signatures(m, n, K, p_sig, lambda, chi_df)
  X <- sig_ls$X
  W <- sig_ls$W
  
  ## generate proportions
  P_star <- matrix(nrow = K, ncol = n)
  
  ## generating the regression part and the latent variable part
  if (first_effect != "cc") {
    ## if in a case_control structrue: half of the samples are affected by a different CT abundances and a different GEP
    P_star <- t(extraDistr::rdirichlet(n = n, alpha = rep(1, K)))
    # for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = rep(1, K))
    Y_reg <- X %*% P_star
  } else {
    ## case-control!
    ## use the improved case-control structure in 2023/09
    ## making the same group assignment as the binary grouping case
    Alpha <- Alpha2 <- t(extraDistr::rdirichlet(n, 1:K))
    Alpha2[, seq(n/2 + 1, n)] <- apply(Alpha2[, seq(n/2 + 1, n)], 2, disturb)
    P_star <- apply(Alpha, 2, softmax)
    P1 <- apply(Alpha2, 2, softmax)
    # Dmat <- P1 - P_star
    Dmat <- matrix(c(rep(0, floor(n/2)), rep(1, ceiling(n/2))), nrow = q, ncol = n)
    # P_star[,  Dmat[1, ] == 0] <- t(extraDistr::rdirichlet(n = floor(n/2), alpha = rep(1, K)))
    # P_star[,  Dmat[1, ] == 1] <- t(extraDistr::rdirichlet(n = ceiling(n/2), alpha = K * 2^seq(K)/sum(2^seq(K))))
    # for (i in 1:floor(n/2)) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = rep(1, K))
    # for (i in seq(floor(n/2) + 1, n)) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = K * 2^seq(K)/sum(2^seq(K)))
    
    ## if the second kind of cell type exhibite different profiles based on a binary assignment
    sig_ls2 <- get_second_signatures(m, n, K, W, 2 * lambda, 2 * chi_df)
    X2 <- extraDistr::rbern(m, p_pert) * sig_ls2$X # only a subset of genes are perturbed 
    Y <- cbind(X %*% P_star[, Dmat[1, ] == 0], (X + X2) %*% P1[, Dmat[1, ] == 1])
  }
 
  ## generate the categorical and continuous latent variables in the D matrix
  # Z <- matrix(nrow = m, ncol = q)
  if (first_effect != "cc") Dmat <- matrix(0, nrow = q, ncol = n)
  if (first_effect == "bin" | first_effect == "miss") Dmat[1, ] <- c(rep(0, floor(n/2)), rep(1, ceiling(n/2)))
  if (first_effect == "con") Dmat[1, ] <- rchisq(n, df = chi_df/K)
  if (q >= 2) {
    ## add another continuous feature (purely positive)
    if (second_effect == "con") {
      Dmat[2, ] <- rchisq(n, df = chi_df/K)
    } else if (second_effect == "bin") {
      Dmat[2, ] <- rep(c(rep(0, floor(n/4)), rep(1, floor(n/4))), 2)
    }
  }  
  if (q >= 3) {
    ## add another continuous feature (could be negative)
    # Dmat[3, ] <- 
  }
  if (q == 4) {
    ## add another discrete feature
    # Dmat[4, ] <- 
  }
  ## choose half of the 0, 1 group to be up regulated
  
  ## encode the latent effects of the binary groups directly instead of using Z * Dmat
  if (first_effect == "small") {
    Y_lat <- get_Y_lat_small_effects(m, n, chi_df, gamma, W)
  } else if (first_effect == "bin") {
    Y_lat <- get_Y_lat_binary(m, n, chi_df/2, gamma, W)
  } else if (first_effect == "con") {
    Y_lat <- get_Y_lat_continuous(m, n, chi_df, gamma, W, Dmat[1, ], p_w2)
  } else if (first_effect == "flat") {
    Y_lat <- get_Y_lat_flat(m, n, chi_df, gamma, W)
  } else if (first_effect == "cc") {
    Y_reg <- X %*% P_star 
    Y_lat <- X2 %*% (P1 - P_star) 
  } else if (first_effect == "miss") {
    Y <- create_miss(Y_reg, Dmat[1, ])
    Y_lat <- Y - Y_reg
  }

  ## get the second latent effect
  if (q >= 2) {
    Y_lat2 <- matrix(0, nrow = m, ncol = n)
    if (second_effect == "con") {
      W2 <- extraDistr::rbern(m, 0.5)
      for (i in 1:n) {
        # Y_lat2[, i] <- W2 * (W * rchisq(m, df = gamma * chi_df/4) + (1 - W) * rchisq(m, df = chi_df/4)) * Dmat[2, i]
        # Y_lat2[, i] <- W2 * (W * rnorm(m, gamma * chi_df, sqrt(chi_df/3)) + (1 - W) * rnorm(m, chi_df, sqrt(chi_df/3))) * Dmat[2, i]
        Y_lat2[, i] <- W2 * (W * rnorm(m, gamma * chi_df/4, sqrt(gamma * chi_df/3)) + (1 - W) * rnorm(m, chi_df/4, sqrt(chi_df/3))) * Dmat[2, i]
      }
    } else if (second_effect == "bin") {
       Y_lat2 <- get_Y_lat_binary2(m, n, chi_df/2, gamma, Dmat[2, ], W)
    }
    Y_lat <- Y_lat + Y_lat2
  }
  # if (q >= 3) {
  #   Y_lat3 <- matrix(0, nrow = m, ncol = n)
  #   W3 <- extraDistr::rbern(m, 0.4)
  # }
  
  ## generate the measurement error
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
  
  ## get the proportion of Y being negative
  p_neg <- sum(Y < 0)/(m * n)
  Y[Y < 0] <- 0
  
  ## gather the results into a list
  if (first_effect == "cc") {
    ls <- list(X = X,
               # Z = Z, # Z is not important!
               X2 = X2,
               P_star = P_star,
               P1 = P1,
               D = Dmat,
               Y = Y,
               E = E,
               Y_lat = Y_lat,
               W = W,
               p_neg = p_neg)
  } else {
    ls <- list(X = X,
               # Z = Z, # Z is not important!
               P_star = P_star,
               D = Dmat,
               Y = Y,
               E = E,
               Y_lat = Y_lat,
               W = W,
               p_neg = p_neg)
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

get_Y_lat_continuous <- function(m, n, chi_df, gamma, W, d, p_w2 = 0.5) {
  Y_lat2 <- matrix(0, nrow = m, ncol = n)
  W2 <- extraDistr::rbern(m, p_w2)
  for (i in 1:n) {
    # Y_lat2[, i] <- W2 * (W * rnorm(m, gamma * chi_df/8, sqrt(gamma * chi_df/3)) + (1 - W) * rnorm(m, chi_df/4, sqrt(chi_df/3))) * d[i]
    # Y_lat2[, i] <- W2 * (W * gamma * runif(m, -1, 1) + (1 - W) * runif(m, -1, 1)) * d[i]
    # Y_lat2[, i] <- W2 * (W * gamma + (1 - W)) * d[i] # the most naive model
    Y_lat2[, i] <- W2 * (W * rnorm(m, gamma * chi_df/4, sqrt(gamma * chi_df/3)) + (1 - W) * rnorm(m, chi_df/4, sqrt(chi_df/3))) * d[i]
  }
  Y_lat2
}

## flat signals in the latent factor
get_Y_lat_flat <- function(m, n, chi_df, gamma, W) {
  Y_lat <- matrix(0, nrow = m, ncol = n)
  for (i in 1:n) {
    Y_lat[, i] <- W * rnorm(m, gamma * chi_df/4, sqrt(gamma * chi_df/3)) + (1 - W) * rnorm(m, chi_df/4, sqrt(chi_df/3))
  }
  Y_lat
}
 
## create missing
create_miss <- function(Y_reg, d) {
  m <- nrow(Y_reg)
  Y_reg[seq(ceiling(m/4)), (d == 0)] <- 0 # mark missing = 0
  Y_reg[seq(m - floor(m/4) + 1, m), (d == 1)] <- 0 # mark missing = 0
  Y_reg
}

# get_sing_vals <- function(true_data, intercept = TRUE, n_sv = 10) {
#   ## see the distribution of top 10 eigenvalues for Y, Y_lat, R
#   Y <- true_data$Y
#   Y_lat <- true_data$Y_lat
#   
#   if (intercept) {
#     X <- model.matrix(~1 + true_data$X)
#   } else {
#     X <- true_data$X
#   }
#   
#   svd_X <- svd(X)
#   U_x <- svd_X$u
#   Sigma_x <- svd_X$d
#   V_x <- svd_X$v
#   B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y   
#   R <- Y - X %*% B_star_hat
#   
#   d_Y <- svd(Y)$d[seq(n_sv)]
#   d_Y_lat <- svd(Y_lat)$d[seq(n_sv)]
#   d_R <- svd(R)$d[seq(n_sv)]
#   
#   mat <- rbind(d_Y, d_Y_lat, d_R)
#   colnames(mat) <- as.character(seq(n_sv))
#   mat
# }

## RE-WRITTEN to be of Sample-Wise PCA!!
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

## For the Bingxin simulation
get_sing_vals2 <- function(true_data, intercept = TRUE, n_sv = 10) {
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

## TODO: write a function to produce bi-plots based on PCA results
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

## use the base R to get the biplot version
# get_bi_plot <- function(true_data, intercept = TRUE) {
#   ## extract Y and Y_lat
#   Y <- true_data$Y
#   Y_lat <- true_data$Y_lat
#   
#   ## calculate the OLS residuals
#   if (intercept) {
#     X <- model.matrix(~1 + true_data$X)
#   } else {
#     X <- true_data$X
#   }
#   svd_X <- svd(X)
#   U_x <- svd_X$u
#   Sigma_x <- svd_X$d
#   V_x <- svd_X$v
#   B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y
#   R <- Y - X %*% B_star_hat
#   
#   ## name the matrices (or the biplot function would call an error)
#   rownames(Y) <- rownames(Y_lat) <- rownames(R) <- paste("gene", as.character(1:m))
#   colnames(Y) <- colnames(Y_lat) <- colnames(R) <- paste("sample", as.character(1:n)) 
#   
#   ## metadata
#   # md <- matrix(as.factor(true_data$D), n, 1)
#   # rownames(md) <- colnames(Y)
#   
#   ## get the svd objects
#   pca_Y <- prcomp(t(Y))
#   pca_Y_lat <- prcomp(t(Y_lat))
#   pca_R <- prcomp(t(R))
#   
#   ## plot the biplots
#   rt_ls <- list(bi_Y = stats::biplot(pca_Y),
#                 bi_Y_lat = stats::biplot(pca_Y_lat),
#                 bi_R = stats::biplot(pca_R))
#   
#   rt_ls
# }

## TODO: write funcitons that will perform simulations stipulated by BZ
## In this, we do not need the latent factors
get_Y_bz <- function(m, n, K, p_sig = 0.5, chi_df = 200, d = 0.15, err = TRUE) {
  ## generate signature matrices
  sig_ls <- get_signatures_bx(m, n, K, p_sig, lambda, chi_df, d)
  X <- sig_ls$X
  W <- sig_ls$W
  X1 <- sig_ls$X1
  
  ## generate proportions
  P_star <- matrix(nrow = K, ncol = n)
  for (i in 1:n) P_star[, i] <- extraDistr::rdirichlet(n = 1, alpha = rep(1, K))
  
  ## obtain the true Y
  Y <- X %*% P_star  
  
  ## getting the errors
  if (err) {
    E <- matrix(0, nrow = m, ncol = n)
    for (i in 1:m) {
      E[i, ] <- rnorm(n, mean(Y_reg[i, ]), sqrt(mean(Y_reg[i, ])/5))
    }
    Y <- Y + E
  } else {
    E <- NA
  }
  
  ## get the proportion of Y being negative
  p_neg <- sum(Y < 0)/(m * n)
  Y[Y < 0] <- 0
  
  ## create a list of returned values
  ls <- list(X = X,
             X1 = X1, # the wrong signature matrix
             # Z = Z, # Z is not important!
             P_star = P_star,
             D = NA,
             Y = Y,
             E = E,
             Y_lat = NA,
             W = W,
             p_neg = p_neg)
  return(ls)
}

## getting the number of latent factors
## TODO: finish this function
# estimate_q <- function(Y, X, intercept = TRUE) {
#   Y <- true_data$Y
#   X <- true_data$X1
#   
#   if (intercept) {
#     X <- model.matrix(~1 + X)
#   } else {
#     X <- true_data$X
#   }
#   
#   ## get the residuals
#   svd_X <- svd(X)
#   U_x <- svd_X$u
#   Sigma_x <- svd_X$d
#   V_x <- svd_X$v
#   B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y   
#   R <- Y - X %*% B_star_hat
#   
#   ## get the sequence of PCAs
#   d_R <- log10(prcomp(t(R))$sdev^2)
#   
#   mat <- rbind(d_Y[seq(n_sv)], d_Y_lat[seq(n_sv)], d_R[seq(n_sv)])
#   colnames(mat) <- as.character(seq(n_sv))
#   mat
# }
