% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tfb_betafn_elasticnet.R
\name{tfb_betafn_elasticnet}
\alias{tfb_betafn_elasticnet}
\title{Elasticnet Beta Function}
\usage{
tfb_betafn_elasticnet(X, y, params)
}
\arguments{
\item{X}{covariate matrix}

\item{y}{outcome vector}

\item{params}{additional arguments}
}
\value{
A numeric vector of coefficients.
}
\description{
This is a helper function that returns a bootstrapped beta for elasticnet
}
\examples{
set.seed(1221)
# fold indices for iris data from tfb_fold
i <- sample(rep(c(TRUE,FALSE),75))
d <- iris[, 5] == "setosa"
X <- as.matrix(iris[!i & !d, 2:4])
y <- iris[!i & !d, 1]
boot <- sample(1:nrow(X), replace=TRUE)
tfb:::tfb_betafn_elasticnet(X[boot,],y[boot],list())
}
\keyword{tfb}
