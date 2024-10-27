#' Generalized TFB WDIM
#'
#' This function returns a wdim and dim within a generalized fold.
#' @param i split indices
#' @param w weights
#' @param d treatment vector
#' @param y outcome vector
#' @param estimand estimand being estimated
#'
#' @returns A list:
#' * `wdim`: The weighted difference in means within the fold.
#' * `dim`: The unweighted difference in means within the fold.
#'
#' @keywords tfb
#' @examples
#' set.seed(1221)
#' i <- sample(rep(c(TRUE,FALSE),75))
#' # weights for iris data from tfb_balance_att
#' w <- c(rep(1,32),c(1.79,1.69,1.92e8,2.05e2,1.57,1.13e2,2.25,1.76,1.46,2.07,
#' 1.46,1.43,2.69e1,5.2,5.66,4.56,2.83,1.7,3.41,2.38e8,2.56,1.02,9.48e-1,1.2,
#' 1.52,1.06,6.69e-1,1.39,1.13,6.63e-1,1.7,1.04,8.14e-1,1.64,8.49e-1,7.09e-1,
#' 8.97e-1,9.19e-1,1.25,1.06,1.69,1.43,1.37)*1e-7)
#' d <- iris[, 5] == "setosa"
#' y <- iris[, 1]
#' tfb:::tfb_estimate_fold(i,w,d,y,"att")


### TFB Summary Statistic Function
tfb_estimate_fold <- function(

  i, # change function name
  w,
  d,
  y,
  estimand

){

  if (estimand == "atc") {
    d <- 1 - d
  }

  wdim <- (mean(w[d[i == 1] == 1] * y[(i == 1) & (d == 1)]) - mean(w[d[i == 1] == 0] * y[(i == 1) & (d == 0)]))
  dim <- mean(y[(i == 1) & (d == 1)]) - mean(y[(i == 1) & (d == 0)])

  return(list(wdim,dim))

}
