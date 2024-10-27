#' Sample-Specific Summary Statistics for TFB
#'
#' This function produces summary statistics for a given sample.
#' @param d treatment vector
#' @param confidence confidence level for the CI produced by tfb
#' @param wdim tfb estimate within the sample
#' @param weights list of vectors of weights
#' @param indices list of fold indices
#' @param v_hats vector of variance estimates
#'
#' @returns A list:
#' * `estimate`: the final tfb estimate for the sample
#' * `variance`: the final variance of the tfb estimate for the sample
#' * `weights`: the weights for the sample
#' * `confidence_interval`: the confidence interval for the estimate for the sample
#' * `p_value`: the p-value for the estimate
#'
#' @importFrom stats qnorm pnorm
#' @keywords tfb
#' @examples
#' d <- iris[, 5] == "setosa"
#' # average of tfb estimates for iris data from tfb_estimate_fold
#' wdim <- -0.0315
#' # weights for iris data from tfb_balance_att
#' weights <- list(
#'   c(rep(1,32),c(1.79,1.69,1.92e8,2.05e2,1.57,1.13e2,2.25,1.76,1.46,2.07,
#'   1.46,1.43,2.69e1,5.2,5.66,4.56,2.83,1.7,3.41,2.38e8,2.56,1.02,9.48e-1,
#'   1.2,1.52,1.06,6.69e-1,1.39,1.13,6.63e-1,1.7,1.04,8.14e-1,1.64,8.49e-1,
#'   7.09e-1,8.97e-1,9.19e-1,1.25,1.06,1.69,1.43,1.37)*1e-7),
#'   c(rep(1,18),c(2.15,2.44,1.92,3.34,2.18,2.18,1.04e1,2.94,2.81,2.54,2.38,
#'   2.14,5.25,4.67,2.03,2.6,2.35,1.69,2.38,2.51,2.13,2.28,3.71,2.36,2.72,
#'   2.72,2.6,5.7e9,3.2,1.71,1.32,1.43,1.36,1.08,2.25,1.15,1.3,1.8,1.58,1.76,
#'   1.76,1.1,1.9,1.98,1.43,1.21,1.44,1.7,2.02,1.58,1.48,1.8,1.71,1.34,1.45,
#'   1.75,1.75)*1e-8)
#' )
#' set.seed(1221)
#' i <- sample(rep(c(TRUE,FALSE),75))
#' indices <- list(i,!i) # indices for iris datafrom tfb_fold
#' v_hats <- c(0.021,0.0903) # variances for iris data from tfb_variance_att
#' tfb:::tfb_summary_sample(d,0.95,wdim,weights,indices,v_hats)

tfb_summary_sample <- function(

  d,
  confidence,
  wdim,
  weights,
  indices,
  v_hats

){

  v_hat <- 1/(length(v_hats)^2) * sum(v_hats)

  w <- rep(1,length(d))
  for (i in seq_along(indices)) {
    w[indices[[i]] == 1] <- weights[[i]]
  }

  ci_lower <- wdim - qnorm((1 + confidence) / 2) * sqrt(v_hat)
  ci_upper <- wdim + qnorm((1 + confidence) / 2) * sqrt(v_hat)

  p_val = min(pnorm(wdim, 0, sqrt(v_hat)), 1 - pnorm(wdim, 0, sqrt(v_hat)))

  out <- list(
      estimate = wdim,
      variance = v_hat,
      weights = w,
      confidence_interval = list(lower = ci_lower, upper = ci_upper),
      p_value = p_val
    )

  return(out)

}
