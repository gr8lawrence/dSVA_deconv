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

## use RUV to obtain a cleaner solution
## TODO: check how to re-assemble the Y after removing batch effects
# ctl <- rep(TRUE, m)
# m1 <- ceiling(m/4)
# m2 <- floor(m/4)
# ctl[c(seq(m1), seq(m - m2 + 1, m))] <- FALSE
# fit <- ruv::RUV4(Y = log10(t(Y)), X = matrix(1, n, 1), Z = NULL, eta = NULL, ctl = ctl, k = q, include.intercept = FALSE)
# Y_rm <- t(exp(log10(t(Y)) - fit$W %*% fit$alpha))


## use LIMMA - log-space removal and then exponential
## examples
# y <- matrix(rnorm(10*9),10,9)
# y[,1:3] <- y[,1:3] + 5
# batch <- c("A","A","A","B","B","B","C","C","C")
# y2 <- removeBatchEffect(y, batch)
# par(mfrow=c(1,2))
# boxplot(as.data.frame(y),main="Original")
# boxplot(as.data.frame(y2),main="Batch corrected")

## Our example
log_Y <- log(Y)
batch <- true_data$D[1, ]
log_Y2 <- limma::removeBatchEffect(x = log_Y, batch = batch)
Y2 <- exp(log_Y2)

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
## Below is Seunggeung's code for dSVA
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