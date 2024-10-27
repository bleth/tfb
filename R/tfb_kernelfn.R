#' Default Kernel Function
#'
#' This is an inner product function, taking two rows of X and returning a value. The default is the gaussian kernel with Hainmueller and Hazlett's suggested bandwidth parameter of `ncol(X)`. However, if a non-default `kernel` parameter is used when fitting KRLS, the user will have to supply the equivalent kernel function in order for `tfb` to return correct estimates. This is because KRLS standardizes its matrices for simpler computation, so in order to cross-fit KRLS `tfb` needs to be able to re-standardize the kernel matrix used in optimization.
#' @param X_i ith row of X
#' @param X_j jth row of X
#'
#' @returns A double, the guassian kernel of `X_i` and `X_j`.
#'
#' @keywords tfb
#' @examples
#' X <- as.matrix(iris[,2:4])
#' tfb:::tfb_kernelfn(X[1,],X[2,])


### Gaussian Kernel Function
tfb_kernelfn <- function(
    X_i,
    X_j
){
  bandwidth <- length(X_i)
  sqr_diff = sum((X_i - X_j)^2)
  return(exp(-sqr_diff/bandwidth))
}
