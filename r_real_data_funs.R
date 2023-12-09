## functions for performing the real data analysis

## function for running the whole dSVA procedure
run_dSVA <- function(bulk_mat, sig_mat, intercept = TRUE, alg = "nnls", solver = "lsei", method = "be", ...) {
  ## estimate q
  q_hat <- estimate_n_comp(Y = bulk_mat, Theta = sig_mat, intercept = TRUE, method = "be", ...)
  
  ## estimate P
  P <- dsva_for_sim(Y = bulk_mat, 
                    Theta = sig_mat, 
                    n_comp = q_hat, 
                    intercept = TRUE, 
                    alg = alg, 
                    solver = solver) # using pnnls the results get interesting
  rownames(P) <- colnames(sig_mat)
  colnames(P) <- colnames(bulk_mat)
  return(P)
}

## function for getting the long form of estimated P
get_long_form <- function(P_hat) {
  P_df <- as_tibble(t(P_hat), rownames = "subject_id")
  P_long <- P_df %>% 
    pivot_longer(-subject_id, 
                 names_to = "Cell_type",
                 values_to = "Proportion") %>% 
    left_join(bulk_cols, by = "subject_id")
  P_long
}