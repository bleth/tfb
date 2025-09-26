#' TFB BART Initial Regression
#'
#' This function fits a BART model on a specific fold of the data and returns related information.
#' @param X covariate matrix
#' @param y outcome vector
#' @param d treatment vector
#' @param i split indices
#' @param bstrap_cov whether to bootstrap covariance
#' @param bstrap_reps how many times to bootstrap covariance
#' @param reg_d whether to include the treatment as a variable for regression
#' @param treatment target treatment value
#' @param quiet whether to suppress console output
#' @param params additional arguments for the fit
#'
#' @returns
#' A list:
#' * `beta`: The regression coefficients of the model, excluding the intercept. Here, the vector of BART predictions.
#' * `V`: The covariance matrix of the posterior distribution.
#' * `e`: The residuals of the model within the fold.
#' * `sigma2`: The estimate of the model's standard error.
#' * `yhat`: The predicted outcomes within the fold.
#' * `X_tf`: The covariates as they were used to fit the model. Here, the identity matrix
#'
#' @import BART
#' @importFrom stats var sd
#' @keywords tfb
#' @examples
#' X <- as.matrix(iris[, 2:4])
#' d <- iris[, 5] == "setosa"
#' y <- iris[, 1]
#' set.seed(1221)
#' # fold indices for iris data from tfb_fold
#' i <- sample(rep(c(TRUE,FALSE),75))
#' tfb:::tfb_target_bart(X,d,y,i,FALSE,NULL,FALSE,"c",TRUE,list())

tfb_target_bart <- function(

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

  ## MAY NEED TO PUT IN SOME CODE HERE TO GET RID OF BART OUTPUT ##

  bartfn <- BART::gbart

  #####################################################################

  if (reg_d) {
    model <- do.call(bartfn, c(list(x.train = cbind(d,X)[i_out, ], y.train = y[i_out], x.test = cbind(treatment,X)),params))
    beta <- unname(model$yhat.test.mean)[i_in]
    fitted <- unname(model$yhat.test.mean)
  } else {
    model <- do.call(bartfn, c(list(x.train = X[i_out & (d == treatment), ], y.train = y[i_out & (d == treatment)], x.test = X),params))
    beta <- unname(model$yhat.test.mean)[i_in]
    fitted <- unname(model$yhat.test.mean)
  }

  sigma2 <- mean((y[i_out & (d == treatment)] - fitted[i_out & (d == treatment)])^2)
  yhat <- fitted[i_in]
  e <- y[i_in & (d == treatment)] - fitted[i_in & (d == treatment)]
  V <- var(model$yhat.test[, i_in])
  X_tf <- diag(rep(1, sum(i_in)))

  return(list(beta,V,e,sigma2,yhat,X_tf))

}
