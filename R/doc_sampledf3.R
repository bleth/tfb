#' TFB Sample Dataset 3
#'
#' This is a sample dataset for use with `tfb`.
#' The dataset is theoretical, so the covariates, treatments, and outcomes do not have any special meaning.
#' More data of the same distribution can be generated using the dgp found in data-raw/.
#' The treatment assignment mechanism is known for this dataset, which makes it useful for comparing the efficiency of causal inference estimators.
#'
#' @format ## `tfb_sampledf3`
#' There is no treatment effect in this dgp.
#' A numeric matrix with 1,000 rows and 6 columns:
#' \describe{
#'   \item{x1,x2}{covariates}
#'   \item{d}{treatment}
#'   \item{y}{observed outcomes}
#'   \item{y0,y1}{potential outcomes}
#' }
"tfb_sampledf3"
