## This file checks lsei vs. cvxr
true_data <- readRDS("data/example.rds")

## write the function for testing
cvxr_NNLS_ext2 <- function(y, X, k_no, alg = c("nnls", "pnnls"), solver = "ECOS", parallel = FALSE) {
  
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
  sol
}

## let's test them!
y = Y[, 1]
X = model.matrix(~1 + true_data$X)
obj_lsei <- lsei::pnnls(a = X, b = y, k = ncol(X) - ncol(true_data$X), sum = 1)
obj_cvxr <- cvxr_NNLS_ext2(y = y, X = X, k_no = 1, alg = "pnnls", solver = "ECOS")


