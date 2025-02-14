% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tfb_target_ols.R
\name{tfb_target_ols}
\alias{tfb_target_ols}
\title{TFB OLS Initial Regression}
\usage{
tfb_target_ols(
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
)
}
\arguments{
\item{X}{covariate matrix}

\item{d}{treatment vector}

\item{y}{outcome vector}

\item{i}{fold indices}

\item{bstrap_cov}{whether to bootstrap covariance}

\item{bstrap_reps}{how many times to bootstrap covariance}

\item{reg_d}{whether to include the treatment as a variable for regression}

\item{treatment}{target treatment value}

\item{quiet}{whether to suppress console output}

\item{params}{additional named parameters to the fit}
}
\value{
A list:
\itemize{
\item \code{beta}: The regression coefficients of the model, excluding the intercept.
\item \code{V}: The robust covariance matrix.
\item \code{e}: The residuals of the model within the fold.
\item \code{sigma2}: The estimate of the model's standard error.
\item \code{yhat}: The predicted outcomes within the fold.
\item \code{X_tf}: The covariates as they were used to fit the model.
}
}
\description{
This function fits an OLS model on a specific fold of the data and returns related information.
}
\examples{
X <- as.matrix(iris[, 2:4])
d <- iris[, 5] == "setosa"
y <- iris[, 1]
set.seed(1221)
# fold indices for iris data from tfb_fold
i <- sample(rep(c(TRUE,FALSE),75))
tfb:::tfb_target_ols(X,d,y,i,FALSE,NULL,FALSE,"c",TRUE,list())
}
\keyword{tfb}
