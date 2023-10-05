## update dSVA extension for this setting
## what if we assume we know q = 1
## allow an intercept

# Y <- true_data$Y
# Theta <- true_data$X

dsva_for_sim <- function(Y, Theta, n_comp = 0, alg = c("nnls", "pnnls"), intercept = TRUE, test = FALSE, solver = c("lsei", "cvxr"), ...) {
  n <- ncol(Y)
  m <- nrow(Y)
  
  ## add an intercept if chosen
  if (intercept) {
    X <- model.matrix(~1 + Theta)
  } else {
    X <- Theta
  }
  
  # if no hidden factor is estimated, return the P_hat from the simplified model
  # if (q == 0) {
  #   message("Use PA to estimate the number of hidden factors (under development).")
  # }
  if (n_comp == 0) {
    message("Specifying number of latent factors = 0.")
    if (solver == "lsei") {
      if (alg == "pnnls") {
        P_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(Theta), sum = 1)$x})
        if (intercept) P_hat <- P_hat[-1, ]
        
      } else if (alg == "nnls") {
        P_tilde <- apply(Y, 2, function(y) {lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(Theta))$x})
        if (intercept) P_tilde <- P_tilde[-1, ]
        P_hat <- apply(P_tilde, 2, function(x) x/sum(x))
      }
      
    } else if (solver == "cvxr") {
      P_hat <- apply(Y, 2, function(y) cvxr_NNLS_ext(y = y, X = X, k_no = ncol(X) - ncol(Theta), alg = alg))
      
    }
    return(P_hat)
  }
  
  ## step 1: obtain the canonical model residual (using linear regression here instead of PNNLS)
  # svd_X <- svd(X)
  # U_x <- svd_X$u
  # Sigma_x <- svd_X$d
  # V_x <- svd_X$v
  # B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y  
  B_star_hat <- solve(t(X) %*% X) %*% t(X) %*% Y  
  R <- Y - X %*% B_star_hat
  
  ## step 2: svd on the residual space
  svd_R <- svd(R)
  if (n_comp == 1) {
    U_q <- as.matrix(svd_R$u[, 1]) 
    Psi_hat <- rbind(svd_R$d[1] %*% t(svd_R$v)[1, ])
    
  } else {
    U_q <- svd_R$u[, seq(n_comp)]
    Psi_hat <- diag(svd_R$d[seq(n_comp)]) %*% t(svd_R$v)[seq(n_comp), ]
  }
  
  ## step 3: estimate the surrogate variable
  # M_jn <- matrix(1/n, n, n)
  # D_jn <- diag(1, n) - M_jn
  # C2 <- Psi_hat %*% D_jn %*% t(Psi_hat)
  Psi_hat1 <- Psi_hat - apply(Psi_hat, 1, mean)
  C <- Psi_hat1 %*% t(Psi_hat1) 
  if (n_comp == 1) {
    if (1 / C[1, 1] < 1e-15) message("q_hat = 1 and the inverse of Psi_hat %*% D_jn %*% t(Psi_hat) is smaller than 1e-15.")
    # Gamma_hat2 <- U_q + (X %*% B_star_hat %*% D_jn %*% t(Psi_hat)) / C[1, 1]
    Gamma_hat <- U_q + (X %*% B_star_hat %*% t(Psi_hat1)) / C[1, 1]
    
  } else {
    if (1 / kappa(C) < 1e-15) message("q_hat > 1 and Psi_hat %*% D_jn %*% t(Psi_hat) is close to singular, i.e. its reverse condition number < 1e-15.")
    svd_C <- svd(C)  # use svd to invert C to mitigate invertibility issues
    C_inv <- svd_C$u %*% diag(1/svd_C$d) %*% t(svd_C$v)
    # Gamma_hat <- U_q + X %*% B_star_hat %*% D_jn %*% t(Psi_hat) %*% C_inv
    Gamma_hat <- U_q + X %*% B_star_hat %*% t(Psi_hat1) %*% C_inv
  }

  ## step 4: fitting the model again with the surrogate variable (using NNLS or PNNLS)
  X_new <- cbind(Gamma_hat, X)
  
  if (solver == "lsei") {
    if (alg == "pnnls") {
      B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = ncol(X_new) - ncol(Theta), sum = 1)$x})
      P_hat <- B_hat[-seq(ncol(X_new) - ncol(Theta)), ]
      
    } else if (alg == "nnls") {
      B_hat <- apply(Y, 2, function(y) {lsei::pnnls(a = X_new, b = y, k = ncol(X_new) - ncol(Theta), sum = NULL)$x})
      P_tilde <- B_hat[-seq(ncol(X_new) - ncol(Theta)), ]
      P_hat <- apply(P_tilde, 2, function(x) x/sum(x))
    }
    
  } else if (solver == "cvxr") {
    P_hat <- apply(Y, 2, function(y) cvxr_NNLS_ext(y = y, X = X_new, k_no = ncol(X_new) - ncol(Theta), alg = alg, ...))
  }
  
  ## return the P_hat
  if (!test) {
    return(P_hat)
  } else {
    Y_lat_hat <- Y - Theta %*% P_hat
    return(list(P_hat = P_hat,
                Y_lat_hat = Y_lat_hat))
  }
}



## a first pass method to find the number of components
estimate_n_comp_cutoff <- function(R = NULL, Y = NULL, Theta = NULL, intercept = TRUE) {
  
  if (is.null(R)) {
    
    if (is.null(Y) | is.null(Theta)) stop("Y and Theta both need to be specified if R is not.")
    
    ## add an intercept 
    if (intercept) {
      X <- model.matrix(~1 + Theta)
    } else {
      X <- Theta
    }
    
    ## perform first pass regression
    Bhat <- solve(t(X) %*% X) %*% t(X) %*% Y  
    R <- Y - X %*% Bhat
  } 
  
  ## use the gap of log singular values
  d_R <- log(prcomp(t(R), scale = TRUE, center = TRUE)$sdev)
  d <- diff(-d_R, lag = 1)
  d <- d[-length(d)]
  n_comp <- which(d == max(d))
  n_comp
}


estimate_n_comp_be <- function(R = NULL, Y = NULL, Theta = NULL, intercept = TRUE, B = 19, seed = NULL, threshold = .05) {
  
  if (is.null(R)) {
    
    if (is.null(Y) | is.null(Theta)) stop("Y and Theta both need to be specified if R is not.")
    
    ## add an intercept 
    if (intercept) {
      X <- model.matrix(~1 + Theta)
    } else {
      X <- Theta
    }
    
    ## perform first pass regression
    Bhat <- solve(t(X) %*% X) %*% t(X) %*% Y  
    R <- Y - X %*% Bhat
  } 
  
  ## use the PA method by Buja & Eyuboglu (1992)
  eigen0 <- prcomp(t(R), scale = TRUE, center = TRUE)$sdev # scaling for a fair comparison 
  nge <- rep(0, length(eigen0)) # TODO: make a B x n matrix for paralleling
  
  if (!is.null(seed)) set.seed(seed)
  for (b in 1:B) {
    
    ## start from the original residuals
    R_tmp <- R
    
    ## rotate genes within columns 2 thru n
    for (i in 2:ncol(R)) R_tmp[, i] <- sample(R_tmp[, i], nrow(R))
    
    ## calculate eigenvalues
    eigen_tmp <- prcomp(t(R_tmp), scale = TRUE, center = TRUE)$sdev
    
    ## compare eigen_tmp to eigen0
    nge <- nge + as.numeric(eigen_tmp > eigen0)
  }
  
  ## calculate cumulative p_values
  pvals <- (nge + 1)/(B + 1)
  
  ## find the significant cutoff
  inds <- which(p_vals[-ncol(R)] < threshold)
  n_comp <- inds[length(inds)]
  
  n_comp
}

## an upper-level function to include both methods
estimate_n_comp <- function(R = NULL, Y = NULL, Theta = NULL, method = "be", B = 19, seed = NULL, intercept = TRUE, ...) {
  
  if (method == "be") {
    n_comp <- estimate_n_comp_be(R = R, Y = Y, Theta = Theta, seed = seed, intercept = intercept, ...)
  } else if (method == "cutoff") {
    n_comp <- estimate_n_comp_cutoff(R = R, Y = Y, Theta = Theta, intercept = TRUE)
  } else {
    stop("Methods has to be either 'be' or 'cutoff'.")
  }
  
  n_comp
}
  
# estimate_n_comp(R = R, method = "be")
# estimate_n_comp(R = R, method = "cutoff")
# R = NULL
# Y = true_data$Y
# Theta = true_data$X
# intercept = TRUE

## use CVXR to re-build a lsei solver for NNLS and pNNLS
# k_no: the number of components not NN-restricted 
# input one column of Y (y) at a time
cvxr_NNLS_ext <- function(y, X, k_no, alg = c("nnls", "pnnls"), solver = "ECOS", parallel = FALSE) {
  
  ## define the variables
  b <- CVXR::Variable(rows = ncol(X))
  obj <- CVXR::Minimize(0.5 * sum((y - X %*% b)^2))
  prob <- CVXR::Problem(objective = obj)
  
  ## define constraints
  if (alg == "pnnls") {
    ## can solve using ECOS
    E <- diag(c(rep(0, k_no), rep(1, ncol(X) - k_no))) 
    F <- matrix(c(rep(0, k_no), rep(1, ncol(X) - k_no)), nrow = 1)  
    cons <- list(-E %*% b <= 0, F %*% b == 1)
    
  } else if (alg == "nnls") {
    ## can solve using ECOS or OSQP
    E <- diag(c(rep(0, k_no), rep(1, ncol(X) - k_no))) 
    cons <- list(-E %*% b <= 0)
  }
  
  CVXR::constraints(prob) <- cons
  sol <- CVXR::solve(prob, solver = solver, parallel = parallel)
  p_hat <- sol$getValue(b)[-seq(k_no)]
  # print(sol$getValue(b))
  
  ## normalize to satisfy sum-to-one if NNLS
  if (alg == "nnls") p_hat <- p_hat/sum(p_hat)
  
  p_hat
}


## code for tests
# Y <- true_data$Y
# Theta <- true_data$X
# intercept <- TRUE

# estimate_n_comp(true_data$Y, true_data$X)
