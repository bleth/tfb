#' BART Beta Function
#'
#' returns a bootstrapped beta for BART
#' @param X covariate matrix
#' @param y outcome vector
#' @param params additional arguments
#' @keywords tfb
#' @import BART
#' @examples
#' tfb_betafn_bart()

tfb_betafn_bart <- function(X,y,params) {
  model <- do.call(wbart,c(list(x.train = X, y.train = y),params))
  return(colMeans(model$yhat.test))
}
