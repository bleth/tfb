#' Elasticnet Beta Function
#'
#' This is a helper function that returns a bootstrapped beta for elasticnet
#' @param X covariate matrix
#' @param y outcome vector
#' @param params additional arguments
#'
#' @returns A numeric vector of coefficients.
#'
#' @import glmnet
#' @keywords tfb
#' @examples
#' set.seed(1221)
#' i <- sample(rep(c(T,F),75)) # fold indices for iris data from tfb_fold
#' d <- iris[, 5] == "setosa"
#' X <- as.matrix(iris[!i & !d, 2:4])
#' y <- iris[!i & !d, 1]
#' boot <- sample(1:nrow(X), replace=T)
#' tfb:::tfb_betafn_elasticnet(X[boot,],y[boot],list())

tfb_betafn_elasticnet <- function(X,y,params) {
  model <- do.call(cv.glmnet, c(list(x = X, y = y),params))
  return(model$glmnet.fit$beta[, model$index[2]])
}
