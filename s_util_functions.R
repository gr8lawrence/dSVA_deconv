## TODO: write a function that plots heatmap -> for visualizing hidden structures
## function 1: reducing the matrix dimension so plotting becomes easier
# example: genmatred <- redim_matrix(genmat, target_height = 600, target_width = 50)
redim_matrix <- function(
    mat,
    target_height = 100,
    target_width = 100,
    # summary_func = function(x) mean(x, na.rm = TRUE),
    summary = c("mean", "max"),
    output_type = 0.0, #vapply style
    n_core = 1 # parallel processing
) {
  ## choose which summary function to use
  if (summary == "mean") {
    summary_func <- function(x) mean(x, na.rm = TRUE)
  } else if (summary == "max") {
    summary_func <- function(x) max(x, na.rm = TRUE)
  }
  ## take log on values if needed to 
  if(target_height > nrow(mat) | target_width > ncol(mat)) {
    stop("Input matrix must be bigger than target width and height.")
  }
  seq_height <- round(seq(1, nrow(mat), length.out = target_height + 1))
  seq_width  <- round(seq(1, ncol(mat), length.out = target_width  + 1))
  # complicated way to write a double for loop
  do.call(rbind, parallel::mclapply(seq_len(target_height), function(i) { # i is row
    vapply(seq_len(target_width), function(j) { # j is column
      summary_func(
        mat[
          seq(seq_height[i], seq_height[i + 1]),
          seq(seq_width[j] , seq_width[j + 1] )
        ]
      )
    }, output_type)
  }, mc.cores = n_core))
}

## function 2: plot the matrix heatmap
plot_heatmap <- function(mat, plot_title = "Latent Variable Structure", take_log = TRUE) {
  if (take_log) mat <- log(mat + 1)
  image(
    t(mat),
    axes = FALSE,
    col = colorRampPalette(c("white", "darkorange", "black"))(30),
    breaks = c(seq(0, 3, length.out = 30), 100),
    main = plot_title
  )
  box()
}

## experiment
# dev.off()
# Y_lat_redim <- redim_matrix(true_data$Y_lat, summary = "mean", target_height = 500, target_width = 10)
# plot_heatmap(true_data$Y_lat, take_log = TRUE)
# plot_heatmap(redim_matrix(true_data$X, target_height = 500, target_width = 5, summary = "max"), take_log = TRUE)
# 
# 
# library(ComplexHeatmap)
# Heatmap(true_data$Y_lat, cluster_rows = FALSE, cluster_columns = FALSE)
# Heatmap(true_data$Y, column_title = "Bulk expression", name = "Value",
#         cluster_rows = FALSE, cluster_columns = FALSE)
