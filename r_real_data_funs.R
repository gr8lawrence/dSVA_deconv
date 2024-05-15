## functions for performing the real data analysis

## function for running the whole dSVA procedure
run_dSVA <- function(bulk_mat, sig_mat, intercept = TRUE, exclude = NULL, alg = "nnls", solver = "lsei", method = "be", q = NULL, ...) {
  
  if (!is.null(q)) {
    q_hat = q
  } else {
    ## estimate q
    if (method == "none") {
      q_hat <- 0
    } else { 
      if (!is.null(exclude)) {
        q_hat <- estimate_n_comp(Y = bulk_mat, Theta = sig_mat[, -exclude], intercept = intercept, method = method, ...)
      } else {
        q_hat <- estimate_n_comp(Y = bulk_mat, Theta = sig_mat, intercept = intercept, method = method, ...)
      }
    }
  }
  
  print(q_hat)
  ## estimate P
  # P <- dsva_for_sim(Y = bulk_mat, 
  #                   Theta = sig_mat, 
  #                   n_comp = q_hat, 
  #                   intercept = intercept, 
  #                   alg = alg, 
  #                   solver = solver) # using pnnls the results get interesting
  # rownames(P) <- colnames(sig_mat)
  # colnames(P) <- colnames(bulk_mat)
  # q_hat = 1
  # q_hat = 0
  P_ls <- dsva_deconv_v2(Y = bulk_mat, 
                         Theta = sig_mat, 
                         n_comp = q_hat, 
                         intercept = intercept, 
                         exclude = exclude,
                         alg = alg, 
                         solver = solver) # using pnnls the results get interesting
  if (!is.null(exclude)) rownames(P_ls$P_hat) <- colnames(sig_mat[, -exclude]) else rownames(P_ls$P_hat) <- colnames(sig_mat)
  colnames(P_ls$P_hat) <- colnames(bulk_mat)
  return(P_ls)
}

## function for getting the long form of estimated P
get_long_form <- function(P_hat, bulk_cols = bulk_cols) {
  P_df <- as_tibble(t(P_hat), rownames = "subject_id")
  P_long <- P_df %>% 
    pivot_longer(-1, 
                 names_to = "Cell_type",
                 values_to = "Proportion") %>% 
    left_join(bulk_cols, by = "subject_id")
  P_long
}

get_long_form2 <- function(P_hat, bulk_cols = bulk_cols) {
  P_df <- as_tibble(t(P_hat), rownames = "geo_accession")
  P_long <- P_df %>% 
    pivot_longer(-1, 
                 names_to = "Cell_type",
                 values_to = "Proportion") %>% 
    left_join(bulk_cols, by = "geo_accession")
  P_long
}
