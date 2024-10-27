#' Sample Splitting
#'
#' This function splits a covariate matrix into folds for TFB.
#' @param d treatment status vector
#' @param folds number of folds to split the data across
#'
#' @returns
#' A list of logical vectors. The ith vector specifies which observations belong to the ith fold.
#'
#' @keywords tfb
#' @examples
#' d <- iris[, 5] == "setosa"
#' set.seed(1221)
#' tfb:::tfb_fold(d,2)

### Sample splitting
tfb_fold <- function(

  d,
  folds

){
  n <- length(d)

  sampling <- sample(rep(seq_len(folds),(n %/% folds)+1)[1:n])

  out <- list()
  for (fold in seq_len(folds)) {
    out <- c(out,list(sampling == fold))
  }

  return(out)
}
