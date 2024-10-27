#' TFB KRLS Initial Regression
#'
#' This function fits a KRLS model on a specific fold of the data and returns related information.
#' @param X covariate matrix
#' @param y outcome vector
#' @param d treatment vector
#' @param i split indices
#' @param bstrap_cov whether to bootstrap covariance
#' @param bstrap_reps how many times to bootstrap covariance
#' @param reg_d whether to include the treatment as a variable for regression
#' @param treatment target treatment value
#' @param quiet whether to suppress console output
#' @param kernelfn the inner product to restandardize by
#' @param params additional arguments for the fit
#'
#' @returns
#' A list:
#' * `beta`: The regression coefficients of the model, excluding the intercept.
#' * `V`: The covariance matrix.
#' * `e`: The residuals of the model within the fold.
#' * `sigma2`: The estimate of the model's standard error.
#' * `yhat`: The predicted outcomes within the fold.
#' * `X_tf`: The covariates as they were used to fit the model. As in, the kernel matrix.
#'
#' @import KRLS
#' @import stats
#' @keywords tfb
#' @examples
#' X <- as.matrix(iris[, 2:4])
#' d <- iris[, 5] == "setosa"
#' y <- iris[, 1]
#' set.seed(1221)
#' i <- sample(rep(c(T,F),75)) # fold indices for iris data from tfb_fold
#' tfb:::tfb_target_krls(X,d,y,i,F,NULL,F,"c",T,tfb:::tfb_kernelfn,list())

tfb_target_krls <- function(

  X,
  d,
  y,
  i,
  bstrap_cov,
  bstrap_reps,
  reg_d,
  treatment,
  quiet,
  kernelfn,
  params

) {

  i_in <- i

  if (all(i)) {
    i_out <- i
  } else {
    i_out <- !i
  }

  treatment <- treatment == "t"

  # if (quiet) {
  #   krlsfn <- krls_quiet
  #   betafn <- tfb_betafn_krls_quiet
  #   params <- c(print.level = 0, params)
  # } else {
  #   krlsfn <- KRLS::krls
  #   betafn <- tfb_betafn_krls
  # }

  ## REMOVE ONCE A SOLUTION TO THE KRLS CONSOLE OUTPUT IS DECIDED ON ##

  if (quiet) {params <- c(print.level = 0, params)}
  #krlsfn <- KRLS::krls
  betafn <- tfb_betafn_krls

  #####################################################################

  if (reg_d) {
    bandwidth <- ncol(X) + 1
    model <- do.call(krlsfn, c(list(X = cbind(d,X)[i_out, ], y = y[i_out]),params))
    beta <- model$coeffs * stats::sd(y[i_out])
    fitted <- KRLS::predict.krls(model, cbind(treatment,X))$fit
  } else {
    bandwidth <- ncol(X)
    model <- do.call(krlsfn, c(list(X = X[i_out & (d == treatment), ], y = y[i_out & (d == treatment)]),params))
    beta <- model$coeffs * stats::sd(y[i_out & (d == treatment)])
    fitted <- KRLS::predict.krls(model, X)$fit
  }

  sigma2 <- mean((y[i_out & (d == treatment)] - fitted[i_out & (d == treatment)])^2)
  yhat <- fitted[i_in]
  e <- y[i_in & (d == treatment)] - fitted[i_in & (d == treatment)]

  if (reg_d) {
    if (bstrap_cov) {
      V <- tfb_bootstrapped_covariance(cbind(d = treatment,X)[i_out, ],y[i_out],betafn,bstrap_reps,params)
    } else {
      V <- model$vcov.c * stats::var(y[i_out])
    }
  } else {
    if (bstrap_cov) {
      V <- tfb_bootstrapped_covariance(X[i_out & (d == treatment), ],y[i_out & (d == treatment)],betafn,bstrap_reps,params)
    } else {
      V <- model$vcov.c * stats::var(y[i_out & (d == treatment)])
    }
  }

  if (reg_d) {
    center <- attr(scale(cbind(d,X)[i_out, ]), which = "scaled:center")
    scale  <- attr(scale(cbind(d,X)[i_out, ]), which = "scaled:scale")
    X_tf <- tfb_gram_matrix(scale(cbind(d = treatment,X), center, scale)[i_in, ],scale(cbind(d,X), center, scale)[i_out, ], kernelfn) # kernel matrix K
  } else {
    center <- attr(scale(X[i_out & (d == treatment), ]), which = "scaled:center")
    scale  <- attr(scale(X[i_out & (d == treatment), ]), which = "scaled:scale")
    X_sc <- scale(X, center, scale)
    X_tf <- tfb_gram_matrix(X_sc[i_in, ], X_sc[i_out & (d == treatment), ], kernelfn)
  }

  return(list(beta,V,e,sigma2,yhat,X_tf))

}
