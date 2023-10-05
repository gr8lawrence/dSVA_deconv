## extended NNLS/pNNLS
NNLS_ext <- function(Y, Theta, alg = c("nnls", "pnnls"), intercept = TRUE, centralized_residual = TRUE) {
  if (intercept) {
    ## add an intercept to the model matrix
    X <- model.matrix(~1 + Theta)
  } else {
    X <- Theta
  }
  # ncol(X) - ncol(Theta)
  
  if (centralized_residual) {
    ## centralize the residual after a least square problem
    # solve a least square problem
    svd_X <- svd(X)
    U_x <- svd_X$u
    Sigma_x <- svd_X$d
    V_x <- svd_X$v
    B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y   
    Y_hat <- X %*% B_star_hat
    R <- Y - Y_hat
    
    # centralize the residuals
    R_cen <- (diag(1, m) - matrix(1/m, m, m)) %*% R
    Y_cen <- Y_hat + R_cen
  } else {
    Y_cen <- Y
  }
  
  ## perform nnls or pnnls on Y combined with the centralized residuals
  if (alg == "nnls") {
    P_tilde <- apply(Y_cen, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(Theta), sum = NULL)$x})
    if (intercept) P_tilde <- P_tilde[-1, ]
    P_hat <- apply(P_tilde, 2, function(x) x/sum(x))
  } else if (alg == "pnnls") {
    P_hat <- apply(Y_cen, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(Theta), sum = 1)$x})
    if (intercept) P_hat <- P_hat[-1, ]
  }
  
  ## return P_hat
  P_hat
}

# sol <- NNLS_ext(Y = true_data$Y, Theta = true_data$X, alg = "pnnls")

## limma-extended (for one batch, blind to the )
limma_ext <- function(Y, Theta, batch, intercept = TRUE, ...) {
  ## adjust for batch effects on the log-expression
  log_Y <- log(Y)
  log_Y2 <- limma::removeBatchEffect(x = log_Y, batch = batch, ...)
  Y2 <- exp(log_Y2)
  
  ## get P_hat on adjusted Y the using pnnls
  P_hat <- NNLS_ext(Y = Y2, Theta = Theta, alg = "pnnls", intercept = TRUE, centralized_residual = FALSE)
  
  ## return P_hat
  P_hat
}

## RUVr-extended (control gene information needed)
RUVr_ext <- function(Y, Theta, n_comp, ctl_gene = NULL, intercept = TRUE, ...) {
  # Y <- true_data$Y
  # Theta <- true_data$X
  
  ## control genes are listed as those corresponding to hidden factor 1 unless otherwise provided
  if (is.null(ctl_gene)) {
    # m <- nrow(Y)
    # ctl <- rep(TRUE, m)
    # m1 <- ceiling(m/4)
    # m2 <- floor(m/4)
    # ctl[c(seq(m1), seq(m - m2 + 1, m))] <- FALSE
    ## every gene is a control gene
    ctl <- rep(TRUE, m)
  }

  
  ## adjust for batch effects on the original expression
  if (intercept) {
    ## add an intercept to the model matrix
    X <- model.matrix(~1 + Theta)
  } else {
    X <- Theta
  }
  
  ## obtain the canonical model residual
  svd_X <- svd(X)
  U_x <- svd_X$u
  Sigma_x <- svd_X$d
  V_x <- svd_X$v
  B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y   
  R <- Y - X %*% B_star_hat
  
  ## adjust the original Y based on the residuals
  Y2 <- RUVSeq::RUVr(x = Y, cIdx = ctl, k = n_comp, residuals = R, ...)
  
  ## use the described method in RUVr manuscript to finish the method
  ## perform an SVD on the residual
  # svd_R <- svd(R)
  # if (n_comp == 1) {
  #   W_hat <- as.matrix(svd_R$u[, 1]) * svd_R$d[1]
  # } else {
  #   W_hat <- svd_R$u[, seq(n_comp)] %*% diag(svd_R$d[seq(n_comp)])
  # }
  
  ## create a new feature matrix
  # if (intercept) {
  #   X_new <- model.matrix(~1 + cbind(W_hat, Theta))
  # } else {
  #   X_new <- cbind(W_hat, Theta)
  # }
  # B_hat <- apply(Y_cen, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = ncol(X_new) - ncol(Theta), sum = 1)$x})
  # P_hat <- B_hat[-seq(ncol(X_new) - ncol(Theta)), ]
  
  ## get P_hat on adjusted Y the using pnnls
  P_hat <- NNLS_ext(Y = Y2$normalizedCounts, Theta = Theta, alg = "pnnls", intercept = TRUE, centralized_residual = FALSE)
  ## return P_hat
  return(P_hat)
}


## RUV4-extended
## without covariates RUV4 is not useful!
# RUV4_ext <- function(Y, Theta, intercept = TRUE) {
#   if (intercept) {
#     ## add an intercept to the model matrix
#     X <- model.matrix(~1 + Theta)
#   } else {
#     X <- Theta
#   }
  
  ## start with an initial pNNLS estimate
  # current_P <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(Theta), sum = 1)$x})
  
  ## perform RUV4 (with an intercept as covariates, not ideal...)
  
  # m <- nrow(Y)
  # ctl <- rep(TRUE, m)
  # m1 <- ceiling(m/4)
  # m2 <- floor(m/4)
  # ctl[c(seq(m1), seq(m - m2 + 1, m))] <- FALSE
  # fit <- ruv::RUV4(Y = log(t(Y)), X = matrix(1, n, 1), Z = NULL, eta = NULL, ctl = ctl, k = q, include.intercept = FALSE)
  # Y_rm <- t(exp(log(t(Y)) - fit$W %*% fit$alpha))
  
  ## perform pNNLS
  # P_hat <- apply(Y_rm, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(Theta), sum = 1)$x})
  # if (intercept) P_hat <- P_hat[-1, ]
# }

## if when the true data is known
get_p_known <- function(true_data, first_effect, intercept = FALSE) {
  if (intercept) {
    ## adding an intercept also makes the known method better
    X <- model.matrix(~1 + true_data$X)
    if (first_effect == "cc") {
      X3 <- true_data$X + true_data$X2
      X2 <- model.matrix(~1 + X3)
    } 
  } else {
    X <- true_data$X
    if (first_effect == "cc") X2 <- true_data$X + true_data$X2
  }
  if (first_effect == "cc") {
    ## recovering the true proportions under separate gene signatures for the cc structure
    P <- cbind(
      apply(true_data$Y[, true_data$D[1, ] == 0], 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(true_data$X), sum = 1)$x}),
      apply(true_data$Y[, true_data$D[1, ] == 1], 2, function(y) {lsei::pnnls(a = X2, b = y, k = ncol(X) - ncol(true_data$X), sum = 1)$x})
    )
  } else {
    ## for the rest
    P <- apply(true_data$X %*% true_data$P_star + true_data$E, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(true_data$X), sum = 1)$x})
  }
  if (intercept) {
    return(P[-1, ])
  } else {
    return(P)
  }
}

## compare methods
# X <- model.matrix(~ 1 + true_data$X)
# sol_cvxr <- cvxr_NNLS_ext(y = true_data$Y[, 1], X = X, k_no = 0, alg = "pnnls", solver = "ECOS")
# sol_lsei <- lsei::pnnls(a = X, b = true_data$Y[, 1], k = ncol(X) - ncol(true_data$X), sum = 1)
# 
# sol_lsei$x  - sol_cvxr
# 
# sol_cvxr <- cvxr_NNLS_ext(y = true_data$Y[, 1], X = X, k_no = 0, alg = "nnls", solver = "ECOS")
# sol_lsei <- lsei::pnnls(a = X, b = true_data$Y[, 1], k = ncol(X) - ncol(true_data$X), sum = NULL)
# sol_lsei$x  - sol_cvxr

