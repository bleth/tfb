#' TFB Elasticnet Initial Regression
#'
#' This function fits an elasticnet model on a specific fold of the data and returns related information.
#' @param X covariate matrix
#' @param d treatment vector
#' @param y outcome vector
#' @param i split indices
#' @param bstrap_cov whether to bootstrap covariance
#' @param bstrap_reps how many times to bootstrap covariance
#' @param reg_d whether to include the treatment as a variable for regression
#' @param treatment target treatment value
#' @param quiet whether to suppress console output
#' @param params additional arguments to be passed to the fit
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
#' @import glmnet
#' @keywords tfb
#' @examples
#' X <- as.matrix(iris[, 2:4])
#' d <- iris[, 5] == "setosa"
#' y <- iris[, 1]
#' set.seed(1221)
#' i <- sample(rep(c(T,F),75)) # fold indices for iris data from tfb_fold
#' tfb:::tfb_target_elasticnet(X,d,y,i,T,400,F,"c",T,list())

tfb_target_elasticnet <- function(

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
    model <- do.call(glmnet::cv.glmnet, c(list(x = cbind(d,X)[i_out, ], y = y[i_out]),params))
    beta <- model$glmnet.fit$beta[, model$index[2]][-1]
    fitted <- predict(model, cbind(d = treatment,X), s = model$lambda.1se)
  } else {
    model <- do.call(glmnet::cv.glmnet, c(list(x = X[i_out & (d == treatment), ], y = y[i_out & (d == treatment)]),params))
    beta <- model$glmnet.fit$beta[, model$index[2]]
    fitted <- predict(model, X, s = model$lambda.1se)
  }

  sigma2 <- mean((y[i_out & (d == treatment)] - fitted[i_out & (d == treatment)])^2)
  yhat <- fitted[i_in]
  e <- y[i_in & (d == treatment)] - fitted[i_in & (d == treatment)]

  V <- if (reg_d) {
    tfb_bootstrapped_covariance(cbind(d,X)[i_out & (d == treatment), ], y[i_out & (d == treatment)], tfb_betafn_elasticnet, bstrap_reps, params)[-1,-1]
  } else {
    tfb_bootstrapped_covariance(X[i_out & (d == treatment), ], y[i_out & (d == treatment)], tfb_betafn_elasticnet, bstrap_reps, params)
  }

  X_tf <- X[i_in, ]

  return(list(beta,V,e,sigma2,yhat,X_tf))

}
