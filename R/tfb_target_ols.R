#' TFB OLS Initial Regression
#'
#' This function fits an OLS model on a specific fold of the data and returns related information.
#' @param X covariate matrix
#' @param d treatment vector
#' @param y outcome vector
#' @param i fold indices
#' @param bstrap_cov whether to bootstrap covariance
#' @param bstrap_reps how many times to bootstrap covariance
#' @param reg_d whether to include the treatment as a variable for regression
#' @param treatment target treatment value
#' @param quiet whether to suppress console output
#' @param params additional named parameters to the fit
#'
#' @returns
#' A list:
#' * `beta`: The regression coefficients of the model, excluding the intercept.
#' * `V`: The covariance matrix.
#' * `e`: The residuals of the model within the fold.
#' * `sigma2`: The estimate of the model's standard error.
#' * `yhat`: The predicted outcomes within the fold.
#' * `X_tf`: The covariates as they were used to fit the model.
#'
#' @importFrom stats lm vcov
#' @keywords tfb
#' @examples
#' X <- as.matrix(iris[, 2:4])
#' d <- iris[, 5] == "setosa"
#' y <- iris[, 1]
#' set.seed(1221)
#' # fold indices for iris data from tfb_fold
#' i <- sample(rep(c(TRUE,FALSE),75))
#' tfb:::tfb_target_ols(X,d,y,i,FALSE,NULL,FALSE,"c",TRUE,list())

tfb_target_ols <- function(

  X,
  d,
  y,
  i,
  bstrap_cov,
  bstrap_reps,
  reg_d,
  treatment,
  quiet,
  params

) {

  i_in <- i

  if (all(i)) {
    i_out <- i
  } else {
    i_out <- !i
  }

  treatment <- treatment == "t"

  # for a given fold, train model on other folds
  if (reg_d) {
    model <- do.call(lm, c(list(formula = y[i_out] ~ d[i_out] + X[i_out, ]),params))
    coeffs <- unname(c(model$coefficients[1] + model$coefficients[2] * treatment, model$coefficients[-(1:2)]))
  } else {
    model <- do.call(lm, c(list(formula = y[i_out & (d == treatment)] ~ X[i_out & (d == treatment), ]),params))
    coeffs <- unname(model$coefficients)
  }

  beta <- coeffs[-1]
  sigma2 <- mean((y[i_out & (d == treatment)] - (cbind(1, X) %*% coeffs)[i_out & (d == treatment)])^2)
  yhat <- (cbind(1, X) %*% coeffs)[i_in]
  e <- y[i_in & (d == treatment)] - yhat[d[i_in] == treatment]

  V <- if (bstrap_cov) {
    if (reg_d) {
      tfb_bootstrapped_covariance(cbind(d,X)[i_out, ], y[i_out], tfb_betafn_ols, bstrap_reps, params)[-1,-1]
    } else {
      tfb_bootstrapped_covariance(X[i_out & (d == treatment), ], y[i_out & (d == treatment)], tfb_betafn_ols, bstrap_reps, params)
    }
  } else {
    if (reg_d) {
      vcov(model)[-(1:2),-(1:2)]
    } else {
      vcov(model)[-1, -1]
    }
  }

  X_tf <- X[i_in, ]

  return(list(beta,V,e,sigma2,yhat,X_tf))

}
