#' TFB BART Initial Regression
#'
#' Returns initial regression information for BART in a specific fold of the data
#' @param X covariate matrix
#' @param d treatment vector
#' @param y outcome vector
#' @param i split indices
#' @param bstrap_cov whether to bootstrap covariance
#' @param bstrap_reps how many times to bootstrap covariance
#' @param reg_d whether to include the treatment as a variable for regression
#' @param treatment target treatment value
#' @param params additional named parameters to the fit
#' @import BART
#' @keywords tfb
#' @examples
#' tfb_target_bart()

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

  # for a given fold, train model on other folds
  if (reg_d) {
    model <- do.call(wbart, c(list(x.train = cbind(d,X)[i_out, ], y.train = y[i_out], x.test = cbind(treatment,X[i_in, ])),params))
    sigma2 <- mean((y[i_out] - colMeans(model$yhat.train))^2)
  } else {
    model <- do.call(wbart, c(list(x.train = X[i_out & (d == treatment), ], y.train = y[i_out & (d == treatment)], x.test = X[i_in, ]),params))
    sigma2 <- mean((y[i_out & (d == treatment)] - colMeans(model$yhat.train))^2)
  }

  beta <- colMeans(model$yhat.test)
  yhat <- beta
  e <- y[i_in & (d == treatment)] - yhat[d[i_in] == treatment]

  if (reg_d) {
    V <- tfb_bootstrapped_covariance(cbind(d,X)[i_out, ], y[i_out], tfb_betafn_bart, bstrap_reps, c(x.test = cbind(treatment,X[i_in, ]), params))
  } else {
    V <- tfb_bootstrapped_covariance(X[i_out & (d == treatment), ], y[i_out & (d == treatment)], tfb_betafn_bart, bstrap_reps, c(x.test = X[i_in, ], params))
  }

  X_tf <- diag(1,nrow = nrow(X[i_in, ]))

  return(list(beta,V,e,sigma2,yhat,X_tf))

}
