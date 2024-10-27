#' Bootstrapping The Covariance Matrix
#'
#' This function bootstraps a given fit's covariance matrix by re-sampling with replacement and re-fitting the model repeatedly.
#' @param X covariates
#' @param y responses
#' @param betafn function taking X,y that returns beta
#' @param bstrap_reps number of repetitions
#' @param params additional arguments for the fit
#'
#' @returns A numeric matrix, the bootstrapped covariance matrix of the model.
#'
#' @import stats
#' @keywords tfb
#' @examples
#' set.seed(1221)
#' # fold indices for iris data from tfb_fold
#' i <- sample(rep(c(TRUE,FALSE),75))
#' d <- iris[, 5] == "setosa"
#' X <- as.matrix(iris[!i & !d, 2:4])
#' y <- iris[!i & !d, 1]
#' tfb:::tfb_bootstrapped_covariance(X,y,tfb:::tfb_betafn_elasticnet,100,list())


### KRLS Initial Regression Function
tfb_bootstrapped_covariance <- function(

  X,              # covariate matrix
  y,              # response vector
  betafn,         # function taking X,y that returns beta
  bstrap_reps,    # number of repetitions
  params          # additional arguments

){

  init <- F

  for (i in 1:bstrap_reps) {
    boot <- sample(1:nrow(X), replace=T)
    if (init) {
      betas <- rbind(betas,matrix(do.call(betafn,list(X=X[boot,],y=y[boot],params=params)),nrow = 1))
    } else {
      betas <- matrix(do.call(betafn,list(X=X[boot,],y=y[boot],params=params)),nrow = 1)
      init <- T
    }
  }

  return(stats::cov(betas))
}
