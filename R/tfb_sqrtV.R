#' TFB Covariance Matrix Square Root
#'
#' This function returns the square root of the covariance matrix or an approximation.
#' @param V covariance matrix
#' @param ev_approx whether to approximate the eigenpairs
#' @param ev_n how many eigenvalues to approximate
#' @param quiet whether to suppress console output
#' @param sample sample of the current iteration
#' @param fold fold of the current iteration
#' @param treatment target treatment value of the current iteration
#'
#' @returns A numeric matrix, the square root of `V`.
#'
#' @import RSpectra
#' @keywords tfb
#' @examples
#' # covariance matrix for iris data from tfb_target_ols
#' V <- matrix(
#'   c(0.0319,-0.00368,-0.00433,
#'   -0.00368,0.0122,-0.0175,
#'   -0.00433,-0.0175,0.0436),nrow=3
#' )
#' tfb:::tfb_sqrtV(V,FALSE,NULL,TRUE,1,1,"c")

tfb_sqrtV <- function(

  V,
  ev_approx,
  ev_n,
  quiet,
  sample,
  fold,
  treatment

) {

  if (!quiet) {
    if (treatment == "c") {
      group <- "control"
    } else {
      group <- "treatment"
    }
    iteration <- paste0("Sample ", sample, ", fold ", fold, ", estimation of the CEF of the ", group, " potential outcomes: ")
  }

  if (ev_approx) {
    if (is.null(ev_n)) {ev_n <- min(5,nrow(V))}
    if (ev_n >= nrow(V)) {
      ev_n <- nrow(V)
      if (!quiet) {message(paste0(iteration, "`ev_n` exceeds size of the covariance matrix, performing full decomposition."))}
    }
    decomposition <- RSpectra::eigs_sym(V, ev_n)
  } else {
    decomposition <- eigen(V, T)
  }

  evalues <- decomposition$values
  evectors <- decomposition$vectors
  sqrtV <- evectors %*% diag(sqrt(evalues)) %*% t(evectors)

  if (ev_approx & !quiet) {message(paste0(iteration, "estimated eigenvalues account for ", 100 * sum(evalues) / sum(diag(V)), " percent of the total variance."))}

  return(sqrtV)

}
