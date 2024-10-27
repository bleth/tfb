#' Final Summary Statistics for TFB
#'
#' This function produces the final summary statistics for TFB across all samples.
#' @param d treatment vector
#' @param y outcome vector
#' @param estimand causal inference estimand being estimated.
#' @param confidence confidence level for the CI produced by tfb
#' @param wdims_samples vector of wdims, one for each sample
#' @param v_hats_samples vector of variances, one for each sample
#'
#' @returns A list:
#' * `estimate`: the final tfb estimate
#' * `variance`: the final variance of the tfb estimate
#' * `confidence_interval`: the confidence interval for the estimate
#' * `p_value`: the p-value for the estimate
#' * `dim`: the naive difference in means for the data
#'
#' @import stats
#' @keywords tfb
#' @examples
#' d <- iris[, 5] == "setosa"
#' y <- iris[, 1]
#' wdims_samples <- -0.0315 # wdim for the iris data from tfb_summary_sample
#' v_hats_samples <- 0.0278 # variance for the estimate for the iris data from tfb_summary_sample
#' tfb:::tfb_summary(d,y,"att",0.95,wdims_samples,v_hats_samples)

tfb_summary <- function(

  d,
  y,
  estimand,
  confidence,
  wdims_samples,
  v_hats_samples

){

  if (estimand == "atc") {
    d <- 1 - d
  }

  dim <- mean(y[d==1]) - mean(y[d==0])

  wdim <- stats::median(wdims_samples)

  v_hat <- stats::median(v_hats_samples + (wdims_samples - wdim)^2/length(d))

  ci_lower <- wdim - stats::qnorm((1 + confidence) / 2) * sqrt(v_hat)
  ci_upper <- wdim + stats::qnorm((1 + confidence) / 2) * sqrt(v_hat)

  p_val = min(stats::pnorm(wdim, 0, sqrt(v_hat)), 1 - stats::pnorm(wdim, 0, sqrt(v_hat)))

  out <- list(
    final = list(
      estimate = wdim,
      variance = v_hat,
      confidence_interval = list(lower = ci_lower, upper = ci_upper),
      p_value = p_val,
      dim = dim
    )
  )

  return(out)

}
