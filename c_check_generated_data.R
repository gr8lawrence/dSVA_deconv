n <- 20
m <- 1000
K <- 5
q <- 2
err <- TRUE
p_sig <- 0.8
lambda <- 5
gamma <- 3

true_data <- dSVA_model_sim_intercept(m, n, K, q, p_sig, lambda, gamma, err = err)
Y <- true_data$Y
Theta <- true_data$X
X <- model.matrix(~1 + Theta)
var(c(true_data$Y_lat))/var(c(true_data$Y - true_data$Y_lat)) # the variance of Y_lat cannot dominate!!

## solve a LS and get the residual
svd_X <- svd(X)
U_x <- svd_X$u
Sigma_x <- svd_X$d
V_x <- svd_X$v
B_star_hat <- V_x %*% diag(1/Sigma_x^2) %*% t(V_x) %*% t(X) %*% Y   
R <- Y - X %*% B_star_hat

## check the singular values of R
svd_R <- svd(R)
svd_R$u

# all.equal(R, svd_R$u %*% diag(svd_R$d) %*% t(svd_R$v)) 
dSVA <- function (Y, X, ncomp = 0) 
{
  X1 <- model.matrix(~1 + X)
  n <- dim(X1)[1]
  if (ncomp == 0) {
    mod <- X1
    id <- num.sv(t(Y), mod, B = 20, method = "be", seed = NULL)
    ncomp <- id
  }
  q = ncomp
  Bhat.Org = solve(t(X1) %*% X1) %*% (t(X1) %*% Y)
  Ehat = Y - X1 %*% Bhat.Org
  current.svd = svd(Ehat)
  u = current.svd$u
  d = diag(current.svd$d)
  v = current.svd$v
  if (q > 1) {
    Delta = d[1:q, 1:q] %*% (t(v)[1:q, ])
  }
  else {
    Delta = rbind(d[1:q, 1:q] * (t(v)[1:q, ]))
  }
  Delta1 <- t(Delta - apply(Delta, 1, mean))
  Bhat.Org1 <- Bhat.Org - apply(Bhat.Org, 1, mean)
  rho_est1 <- (solve(t(Delta1) %*% Delta1) %*% t(Delta1) %*% 
                 t(Bhat.Org1))
  Bhat.Adj <- Bhat.Org - t(rho_est1) %*% Delta
  MZ1 <- rho_est1
  Z.est = u[, 1:q] + X1 %*% t(MZ1)
  re <- Get_Pvalue(Y, X1, Z.est)
  re$Z = Z.est
  re$ncomp = ncomp
  return(re)
}