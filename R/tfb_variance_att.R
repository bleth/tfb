#' Generalized WDIM Variance for the ATT
#'
#' This function returns an estimate of variance for the TFB estimate of the ATT within a fold.
#' @param yhat_c predicted values produced by f_0
#' @param e_c residuals of f_0
#' @param i fold indices
#' @param w weights
#' @param d treatment vector
#' @param wdim final tfb estimate
#' @param y response vector
#'
#' @returns A double, the estimated variance of the WDIM-ATT within a fold.
#'
#' @keywords tfb
#' @examples
#' X <- as.matrix(iris[, 2:4])
#' y <- iris[, 1]
#' d <- iris[, 5] == "setosa"
#' set.seed(1221)
#' i <- sample(rep(c(TRUE,FALSE),75))
#' # coefficients for iris data from tfb_target_ols
#' yhat_c <- (cbind(1, X) %*% c(1.8,0.565,0.812,-0.686))[i]
#' e_c <- y[i & !d] - yhat_c[!d[i]]
#' # weights for iris data from tfb_balance_att
#' w <- c(rep(1,32),c(1.79,1.69,1.92e8,2.05e2,1.57,1.13e2,2.25,1.76,1.46,
#' 2.07,1.46,1.43,2.69e1,5.2,5.66,4.56,2.83,1.7,3.41,2.38e8,2.56,1.02,
#' 9.48e-1,1.2,1.52,1.06,6.69e-1,1.39,1.13,6.63e-1,1.7,1.04,8.14e-1,1.64,
#' 8.49e-1,7.09e-1,8.97e-1,9.19e-1,1.25,1.06,1.69,1.43,1.37)*1e-7)
#' wdim <- -0.0315 # average wdim for iris data from tfb_estimate_fold
#' tfb:::tfb_variance_att(yhat_c,e_c,i,w,d,wdim,y)

tfb_variance_att <- function(

  yhat_c,
  e_c,
  i,
  w,
  d,
  wdim,
  y

){

  v_hat <- (1 / (sum(d[i == 1])^2)) * sum((y[(i == 1) & (d == 1)] - yhat_c[d[i == 1] == 1] - wdim)^2) + (1 / (sum(1 - d[i == 1])^2)) * sum((w[d[i == 1] == 0] * e_c)^2)

  return(v_hat)

}
