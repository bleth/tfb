#' Generalized WDIM Variance for the ATE
#'
#' This function returns an estimate of variance for the TFB estimate of the ATE within a fold.
#' @param yhat_c predicted values produced by f_0
#' @param yhat_t predicted values produced by f_1
#' @param e_c residuals of f_0
#' @param e_t residuals of f_1
#' @param i fold indices
#' @param w weights
#' @param d treatment vector
#' @param wdim final tfb estimate
#'
#' @returns A double, the estimated variance of the WDIM-ATE within a fold.
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
#' yhat_t <- (cbind(1, X) %*% c(2.28,0.532,0.632,-0.251))[i]
#' e_c <- y[i & !d] - yhat_c[!d[i]]
#' e_t <- y[i & d] - yhat_t[d[i]]
#' # weights for iris data from tfb_balance_ate
#' w <- c(9.66e-1,8.81e-1,1.27,1.45,9.29e-1,1.18,1.76,1.97,9.33e-1,1.01,
#' 4.87,1.29,1.26,8.71,2.04,6.57e-1,3.97,3.2e8,2.3,1.15,1.69,1.58,
#' 1.62,1.23,7.58e-1,1.21,1.45,1,2.96,1.05,1.76,1.1,1.09e1,4.8e6,
#' 1.94e1,4.83,1.04e1,2.71e8,3.11e1,1.27e1,3.73,2.98e1,7.19,6.8,
#' 3.77e2,8.27,1.07e1,2.05e1,1.54e8,1.9e1,1.09e1,1.17e1,1.04e1,
#' 3.43,4.99,3.71,4.56,4.68,1.73,2.64,3.87,2.3,3.91,5.32,5.15,6.42,
#' 4.88,1.5e1,3.53,2.74,5,5.55,3.62,4.51,6.16)*1e-7
#' wdim <- -0.806 # average wdim for iris data from tfb_estimate_fold
#' tfb:::tfb_variance_ate(yhat_c,yhat_t,e_c,e_t,i,w,d,wdim)

tfb_variance_ate <- function(

  yhat_c,
  yhat_t,
  e_c,
  e_t,
  i,
  w,
  d,
  wdim

){

  v_hat <- (1 / (length(d[i == 1])^2)) * sum((yhat_t - yhat_c - wdim)^2) + (1 / (sum(1-d[i == 1])^2)) * sum((w[d[i == 1] == 0] * e_c)^2) + (1 / (sum(d[i == 1])^2)) * sum((w[d[i == 1] == 1] * e_t)^2)

  return(v_hat)

}
