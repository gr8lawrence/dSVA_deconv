A = matrix(0.1, 10, 40)
A[, 21:40] = c(rep(0.05, 5), rep(0.15, 5))
B = cbind(matrix(0, 10, 20), A[, 21:40])
cov(t(A), t(B))


C = A
C[, 21:40] = c(rep(0, 5), rep(0.2, 5))
D = cbind(matrix(0, 10, 20), C[, 21:40])
cov(t(C), t(D))


E = matrix(0.1, 10, 40)
E[, 39:40] = c(rep(0.025, 5), rep(0.175, 5))
FM = cbind(matrix(0, 10, 38), E[, 39:40])
cov(t(E), t(FM))

# str_to_title("ALGORITHMS FOR CELL-TYPE DECONVOLUTION UNDER CONSTRAINTS AND
# LATENT STRUCTURAL UNWANTED VARIATION")
