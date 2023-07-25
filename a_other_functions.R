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
