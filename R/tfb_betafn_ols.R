#' OLS Beta Function
#'
#' This is a helper function that returns a bootstrapped beta for OLS.
#' @param X covariate matrix
#' @param y outcome vector
#' @param params additional arguments
#'
#' @returns A numeric vector of coefficients.
#'
#' @importFrom stats lm
#' @keywords tfb
#' @examples
#' set.seed(1221)
#' # fold indices for iris data from tfb_fold
#' i <- sample(rep(c(TRUE,FALSE),75))
#' d <- iris[, 5] == "setosa"
#' X <- as.matrix(iris[!i & !d, 2:4])
#' y <- iris[!i & !d, 1]
#' boot <- sample(1:nrow(X), replace=TRUE)
#' tfb:::tfb_betafn_ols(X[boot,],y[boot],list())

tfb_betafn_ols <- function(X,y,params) {
  model <- do.call(lm,c(list(formula = y ~ X),params))
  return(unname(model$coefficients[-1]))
}
