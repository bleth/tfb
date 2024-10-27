#' Kernel Matrix for TFB with KRLS
#'
#' This function creates a kernel matrix with size n_rows by n_cols using an inner product. for each row of `X_rows`, its inner product with (similarity to) every row of `X_cols` is calculated.
#' @param X_rows covariate matrix defining rows of K
#' @param X_cols covariate matrix defining cols of K
#' @param kernelfn inner product to use
#'
#' @returns A numeric matrix, the kernel/gram matrix produced by `X_rows` and `X_cols`.
#'
#' @keywords tfb
#' @examples
#' X <- as.matrix(iris[,2:4])
#' d <- iris[, 5] == "setosa"
#' set.seed(1221)
#' # fold indices for iris data from tfb_fold
#' i <- sample(rep(c(TRUE,FALSE),75))
#' tfb:::tfb_gram_matrix(X[i, ], X[!i & d, ], tfb:::tfb_kernelfn)


### Kernel Matrix Function
tfb_gram_matrix <- function(
  X_rows, # input data for rows of K
  X_cols, # input data for cols of K
  kernelfn  # function that takes two rows of X and returns a value
){

  n <- nrow(X_rows) # number of rows in X
  p <- nrow(X_cols) # number of columns in X
  M <- matrix(nrow = n, ncol = p) # basic matrix to alter

  for (i in 1:n) {
    for (j in 1:p) {
      product = tfb_kernelfn(X_rows[i, ], X_cols[j, ])
      M[i, j] = product
    }
  }

  return(M)

}
