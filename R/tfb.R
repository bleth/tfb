#' Targeted Function Balancing
#'
#' @description `tfb` performs Targeted Function Balancing, a weighting method for causal inference. `tfb` weights observations in your dataset by minimizing targeted function imbalance from an initial regression, and returns a causal estimate as well as summary statistics. `tfb` is intended for use with data from observational studies where inference about the effect of a treatment is desired even though the treatment assignment mechanism is unknown or dependent on the covariates.
#' @param X The covariate dataset, containing all observed confounders of the relationship between the treatment and the outcome. `X` can contain categorical covariates if a dataframe is supplied. A categorical covariate with `n` levels will be converted to `n-1` indicator variables for `tfb`'s optimization. A numeric matrix or a dataframe. No default.
#' @param d The dichotomous treatment vector, which takes the value of `0` or `FALSE` when an observation is in the control group, and the value of `1` or `TRUE` when an observation is in the treatment group. If a dichotomous character vector is supplied, the first value in the vector will be interpreted as indicating treatment. A numeric, logical, or character vector. No default.
#' @param y The outcome vector, containing the measured response for each observation. A numeric vector. No default.
#' @param fit The hypothesized functional form of the relationship between the outcome and the covariates. This decides what regression model will be fit on the data and targeted by the targeted function balancing. Can be set to one or two of the following:
#' * `"ols"`: ordinary least squares ( `lm()` ) will be used.
#' * `"krls"`: kernel regularized least squares ( `krls()` ) will be used.
#' * `"elasticnet"`: an elastic net penalized regression ( `glmnet()` ) will be used.
#'
#' `fit` should only include two values when estimating the ATE and the treatment and control potential outcomes are hypothesized to be produced by different families of functions on the covariates. A character or character vector. No default.
#' @param estimand The causal inference estimand to be estimated via targeted function balancing. Can be set to one of the following:
#' * `"att"`: the average treatment effect on the treated (ATT) will be estimated.
#' * `"atc"`: the average treatment effect on the controlled (ATC) will be estimated.
#' * `"ate"`: the average treatment effect (ATE) will be estimated.
#'
#'  A character. No default.
#' @param quiet If set to `TRUE`, suppresses console output. A logical. Default is `TRUE`.
#' @param folds The number of folds for cross-fitting. The data will be partitioned uniformly at random into this many folds. It is suggested to use at least two folds, as using only one makes cross-fitting impossible and may introduce bias. An integer. Default is `2`.
#' @param samples The number of times to sample the data; the number of times the data is split into folds. An integer. Default is `1`.
#' @param bstrap_cov Whether to bootsrap the covariance. If `FALSE`, `tfb` will estimate the covariance matrix via closed formula where applicable. `bstrap_cov` is ignored when there is no closed form for the covariance matrix (elasticnet). Setting this to `TRUE` may significantly extend calculation times. A logical. Default is `FALSE`.
#' @param ... Additional named arguments. See *Details* for more information.
#'
#' @returns
#' A nested list:
#' * `final`: Summary statistics and estimates for the dataset, including:
#' * `sample_i`: Summary statistics for the ith sample, including:
#'
#' Containing generally:
#' * `estimate`: The TFB estimate/weighted difference in means.
#' * `variance`: The variance of the estimate.
#' * `confidence_interval`: The confidence interval for the estimate.
#' * `p_value`: The p-value for the estimate under the null hypothesis that there is no treatment effect.
#' * `dim`: The unweighted difference in means.
#'
#' Containing by sample:
#' * `weights`: TFB's weights.
#'
#' Containing by fold:
#' * `indices`: The fold indices.
#' * `betas`: The covariate coefficients for the targeted function.
#' * `covariances`: The covariance matrices for the targeted function.
#' * `errors`: The residuals within the fold.
#' * `standard_errors`: The variance estimates of error terms within the regression.
#' * `dims`: The unweighted differences in means.
#' * `estimates`: The TFB estimates/weighted differences in means.
#' * `variances`: The variances of the estimates.
#'
#' @details
#' The dynamic dots also allow for a number of additional, more technical parameters:
#' * `bstrap_reps`: How many times to bootstrap the covariance matrix. This will have a significant effect on computation time, especially when targeting KRLS, as the reparametrization of X to its kernel matrix K greatly increases the number of covariates, provided n >> p. An integer. Default is `400.`
#' * `confidence`: The choice of confidence level for the reported confidence intervals. A double. Default is `0.95`.
#' * `chi_q`: The choice of quantile for the chi-squared distribution, which effectively controls the spread of the hypothesized distribution of \eqn{\hat\beta_d} about the true coefficient vector \eqn{\beta_d} in \eqn{\mathbb{R}^n}, supposing the CEF has been correctly specified. `chi_q` is used to specify the squidward constant. A double. Default is `0.95`.
#' * `ev_approx`: Whether to abbreviate calculation of covariance matrices using `RSpectra`, for when computation of the covariance matrix would be particularly taxing (large p). A logical. Default is `False.`
#' * `ev_n`: How many eigenpairs `RSpectra` should calculate. The user is recommended to analyze a complete sample eigenmatrix first, to look for a common breakpoint such that the calculated eigenvalues still account for most all of the total explained variance. If the number is too big with KRLS a full decomposition will be performed instead. An integer. The default is typically `5`.
#' * `kernelfn`: A way to specify the inner product kernel function alternative to the KRLS default (the gaussian kernel). `tfb` has to restandardize all of the KRLS kernel matrices in order for calculations to turn out correctly, but this requires the use of the correct inner product. In short, if the user passes KRLS an alternative `whichkernel` through the dynamic dots, they must also pass `tfb` an inner product matching the alternative kernel, taking two rows of the covariate matrix as input and returning a positive semidefinite similarity score as output. A function. The default is the gaussian kernel, which the user may wish to view as a template: `tfb:::tfb_kernelfn`.
#' * `reg_d`: Whether to incorporate the treatment vector into the covariates, allowing it to be used for prediction in the targeting stage of `tfb`. \eqn{\hat f_0} and \eqn{\hat f_1} are then respectively calculated by forcing \eqn{D} to be a constant term of zero for \eqn{\hat f_0} and one for \eqn{\hat f_1}. This also implies that `tfb` will train \eqn{\hat f_0} and \eqn{\hat f_1} on both treatment and control units. A logical. The default is `False`.
#' * `X_c`,`X_t`: A way to modify training data for \eqn{\hat f_0} and \eqn{\hat f_1} _separately_ of the rest of the dataset. For example, you could fit an \eqn{\hat f_0 = \beta_0+\beta_1x_1+\beta_2x_2^2} by letting `X_c = cbind(X[, 1],X[, 2]^2)` while leaving \eqn{\hat f_1 = \beta_0 + \beta_1x_1+\beta_2x_2}. Numeric matrices. Default is `X`.
#'
#'
#' @import fastDummies
#' @import stats
#' @keywords tfb
#' @export
#' @examples
#'
#' X <- iris[, 2:4]
#' d <- iris[, 5] == "setosa"
#' y <- iris[, 1]
#' # estimating the effect of `Species == "setosa"` on `Sepal.Length` given `Sepal.Width`, `Petal.Length`, and `Petal.Width`
#' tfb(X, d, y, fit = "ols", estimand = "att")$final$estimate
#'

tfb <- function(
    X,              # covariate matrix
    d,              # treatment vector
    y,              # response vector
    fit,            # c("ols", "krls", "elasticnet")
    estimand,       # c("atc", "att", "ate")
    quiet = T,      # suppress console output
    folds = 2,      # number of folds for cross-fitting
    samples = 1,    # number of times to sample and fit the data
    bstrap_cov = F, # whether to bootsrap the covariance
    ...             # additional named arguments
){

  # additional tfb arguments
  args_extended <- c(
    "bstrap_reps",
    "confidence",
    "chi_q",
    "ev_approx",
    "ev_n",
    "kernelfn",
    "reg_d",
    "X_c",
    "X_t"
  )

  # all ellipsis arguments
  args_ellipsis <- names(list(...))

  # arguments for the initial regression
  args_params <- args_ellipsis[!(args_ellipsis %in% args_extended)]

  # assigning extended arguments for local reference
  for (arg in args_extended) {
    if (arg %in% args_ellipsis) {
      assign(arg,list(...)[args_ellipsis %in% arg][[1]])
    }
  }

  # function that checks for integer input
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){abs(x - round(x)) < tol}

  # check basic argument quality
  if (!(is.matrix(X) | is.data.frame(X))) {stop("`X` must be a matrix or data frame.\n")}
  if (sum(sapply(1:ncol(X), function(i) class(X[,i])) %in% c("factor", "character")) > 0) {
    message("Categorical covariates are being expanded to one or more binary covariates for their respective levels.\n")
    x <- fastDummies::dummy_cols(X, remove_first_dummy = T)
  }
  X <- as.matrix(X)
  if (nrow(X) < 1 | ncol(X) < 1) {stop("`X` must contain data.\n")}
  if (!is.vector(d) | !(is.numeric(d) | is.logical(d) | is.character(d))) {stop("`d` must be a numeric, logical, or character vector.\n")}
  if (length(unique(d)) < 2) {stop("`d` cannot be a degenerate random variable.\n")}
  if (length(unique(d)) > 2) {stop("`d` cannot take more than 2 values, as `tfb` has not been generalized to non-dichotomous treatments at this time.\n")}
  if (is.character(d)) {d <- d == unique(d)[1]}
  if (length(d) != nrow(X)) {stop("`d` must be the same length as `X`.\n")}
  if (!is.vector(y) | !(is.numeric(y) | is.logical(y))) {stop("`y` must be a numeric or logical vector.\n")}
  if (length(y) != nrow(X)) {stop("`y` must be the same length as `X`.\n")}
  if (!(estimand %in% c("att","atc","ate"))) {stop("`estimand` must be one of `\"att\",\"atc\",\"ate\".\n")}
  if (estimand == "atc") {d <- 1-d} # inverting treatment for atc
  if (any(!(fit %in% c("ols","krls","elasticnet","bart")))) {stop("`fit` must be one of `\"ols\",\"krls\",\"elasticnet\"`.\n")}
  if (length(fit) > 1 & estimand != "ate") {warning("Multiple fits is only intended for `estimand=\"ate\". Using only the first fit specified.\n")}
  if (length(fit) > 2 & estimand == "ate") {warning("`estimand=\"ate\" uses at most two different fits. Using only the first two fits specified.\n")}
  if (length(fit) == 1 & estimand == "ate") {fit <- rep(fit, 2)}
  if (!is.logical(quiet)) {stop("`quiet` must be a logical.\n")}
  if (!is.wholenumber(folds)) {stop("The number of `folds` must be an integer.\n")}
  if (folds < 1 | folds > nrow(X)) {stop("The number of `folds` must be at least `1` and at most `n`.\n")}
  if (!is.wholenumber(samples)) {stop("The number of `samples` must be an integer.\n")}
  if (samples < 1) {stop("The number of `samples` must be at least `1`.\n")}
  if (!is.logical(bstrap_cov)) {stop("`bstrap_cov` must be a logical.\n")}
  if (all(fit %in% "elasticnet") & bstrap_cov == F) {
    if (quiet == F) {warning("Ignoring user input of `bstrap_cov = F`. The covariance must be bootstrapped when `fit = \"elasticnet\"`.\n")}
  }

  # check extended argument quality
  if (!exists("bstrap_reps")) {
    bstrap_reps <- 400
  } else {
    if (!is.wholenumber(bstrap_reps)) {stop("`bstrap_reps` must be an integer.\n")}
    if (bstrap_reps < 1) {stop("`bstrap_reps` must be at least `1`.\n")}
  }

  if (!exists("confidence")) {
    confidence <- 0.95
  } else {
    if (!is.double(confidence)) {stop("`confidence` must be a proportion.\n")}
    if (confidence <= 0 | confidence >= 1) {stop("`confidence` must be a proportion.\n")}
  }

  if (!exists("chi_q")) {
    chi_q <- 0.95
  } else {
    if (!is.double(chi_q)) {stop("`chi_q` must be a proportion.\n")}
    if (chi_q <= 0 | chi_q >= 1) {stop("`chi_q` must be a proportion.\n")}
  }

  if (!exists("ev_approx")) {
    ev_approx <- F
  } else {
    if (!is.logical(ev_approx)) {stop("`ev_approx` must be a logical.\n")}
    if ("krls" %in% fit & ev_approx) {warning("Use of `ev_approx` with KRLS is not recommended.\n")}
  }

  if (!exists("ev_n")) {
    ev_n <- NULL
  } else {
    if (!is.wholenumber(ev_n)) {stop("`ev_n` must be an integer.\n")}
    if (!("krls" %in% fit) & (ev_n < 1 | ev_n > ncol(X))) {stop(paste0("For `fit` of ", fit, " `ev_n` must be at least 1 and at most ", ncol(X), "\n."))}
    if ("krls" %in% fit & (ev_n < 1)) {stop(paste0("For `fit` of ", fit, " `ev_n` must be at least 1.\n"))}
    if ("krls" %in% fit & !is.null(ev_n)) {warning("The set value of `ev_n` may equal or exceed the maximum dimensions of the eigenmatrix. In the case that this happens, `tfb` will perform a complete eigenvalue decomposition.\n")}
  }

  if (!exists("kernelfn")) {
    kernelfn <- tfb_kernelfn
  } else {
    if (!is.function(kernelfn)) {stop("`kernelfn` must be a function.\n")}
    if (length(formals(kernelfn)) != 2) {stop("`kernelfn` must be an inner product, taking two rows and returning one value.\n")}
  }

  if (!exists("reg_d")) {
    reg_d <- F
  } else {
    if (!is.logical(reg_d)) {stop("`reg_d` must be a logical.\n")}
  }

  if (!exists("X_c")) {
    X_c <- X
  } else {
    if (!is.numeric(X_c) | !is.matrix(X_c)) {stop("`X_c` must be a numeric matrix.\n")}
    if (nrow(X_c) != nrow(X)) {stop("`X_c` must have as many rows as `X`.\n")}
  }

  if (!exists("X_t")) {
    X_t <- X
  } else {
    if (!is.numeric(X_t) | !is.matrix(X_t)) {stop("`X_t` must be a numeric matrix.\n")}
    if (nrow(X_t) != nrow(X)) {stop("`X_t` must have as many rows as `X`.\n")}
  }

  # creating iterator for specific estimand
  treatments <- c("c")
  if (estimand == "ate") {
    treatments <- c(treatments,"t")
    estimand_type <- "ate"
  } else {
    estimand_type <- "att"
  }

  # initializing empty output variables
  output <- list()
  wdims_samples <- v_hats_samples <- c()

  # iterating by sample
  for (sample in seq_len(samples)) {

    # split data into folds
    out <- tfb_fold(d, folds)
    args <- paste0("i_",seq_len(folds))
    for (i in seq_along(out)) {assign(args[i],out[[i]])}

    # iterating by fold
    for (fold in seq_len(folds)) {

      # ate fits twice, att/atc fit once
      for (treatment in treatments) {
        # fitting approximations of the CEF
        args <- c(paste0("X_",treatment),"d","y",paste0("i_",fold),"bstrap_cov","bstrap_reps","reg_d","treatment","quiet")
        if (fit[(treatment == "t") + 1] == "krls") {args <- c(args,"kernelfn")}
        out <- do.call(paste0("tfb_target_",fit[(treatment == "t") + 1]),c(lapply(args,as.symbol),list(params = list(...)[args_ellipsis %in% args_params])))
        args <- paste0(c("beta_","V_","e_","sigma2_","yhat_","X_"),treatment,fold)
        for (i in seq_along(out)) {assign(args[i],out[[i]])}

        # decomposition/square root of the covariance matrix
        args <- c(paste0("V_",treatment,fold), "ev_approx", "ev_n", "quiet", "sample", "fold", "treatment")
        out <- do.call("tfb_sqrtV",lapply(args,as.symbol))
        assign(paste0("sqrtV_",treatment,fold), out)
      }

      # optimization/weighting to minimize objective function
      args <- unlist(lapply(c("X_","beta_","sqrtV_" ,"sigma2_"),paste0,treatments))
      args <- c(paste0(c(args,"i_"),fold), "d", "chi_q", "quiet")
      out <- do.call(paste0("tfb_balance_",estimand_type),lapply(args,as.symbol))
      assign(paste0("w_",fold),out)

      # calculating wdim and dim within the fold
      args <- c("i_","w_")
      args <- c(paste0(args,fold), "d", "y", "estimand")
      out <- do.call(tfb_estimate_fold, lapply(args,as.symbol))
      args <- paste0(c("wdim_","dim_"),fold)
      for (i in seq_along(out)) {assign(args[i],out[[i]])}
    }

    args <- ls()
    wdims <- args[startsWith(args,"wdim_")]
    wdim <- mean(do.call(c,lapply(wdims,as.symbol)))
    wdims_samples <- c(wdims_samples,wdim)

    # calculating variance by fold
    for (fold in seq_len(folds)) {
      args <- unlist(lapply(c("yhat_","e_"),paste0,treatments))
      args <- c(paste0(c(args,"i_","w_"), fold),"d","wdim")
      if (estimand_type == "att") {args <- c(args,"y")}
      out <- do.call(paste0("tfb_variance_", estimand_type), lapply(args,as.symbol))
      assign(paste0("v_hat_",fold),out)
    }

    # calculating summary stats for the sample
    args <- ls()
    v_hats <- args[startsWith(args,"v_hat_")]
    out <- tfb_summary_sample(
      d,confidence,wdim,
      lapply(lapply(args[startsWith(args,"w_")], as.symbol),eval,envir=environment()),
      lapply(lapply(args[startsWith(args,"i_")], as.symbol),eval,envir=environment()),
      do.call(c,lapply(v_hats,as.symbol))
    )
    v_hats_samples <- c(v_hats_samples,out$variance)

    # binding output for the sample
    output <- c(output, list(. = list(
      final = out,
      indices = stats::setNames(lapply(lapply(args[startsWith(args,"i_")], as.symbol),eval,envir=environment()),args[startsWith(args,"i_")]),
      betas = stats::setNames(lapply(lapply(args[startsWith(args,"beta_")], as.symbol),eval,envir=environment()),args[startsWith(args,"beta_")]),
      covariances = stats::setNames(lapply(lapply(args[startsWith(args,"V_")], as.symbol),eval,envir=environment()),args[startsWith(args,"V_")]),
      errors = stats::setNames(lapply(lapply(args[startsWith(args,"e_")], as.symbol),eval,envir=environment()),args[startsWith(args,"e_")]),
      standard_errors = stats::setNames(lapply(lapply(args[startsWith(args,"sigma2_")], as.symbol),eval,envir=environment()),args[startsWith(args,"sigma2_")]),
      dims = stats::setNames(lapply(lapply(args[startsWith(args,"dim_")], as.symbol),eval,envir=environment()),args[startsWith(args,"dim_")]),
      estimates = stats::setNames(lapply(lapply(args[startsWith(args,"wdim_")], as.symbol),eval,envir=environment()),args[startsWith(args,"wdim_")]),
      variances = stats::setNames(lapply(lapply(args[startsWith(args,"v_hat_")], as.symbol),eval,envir=environment()),args[startsWith(args,"v_hat_")])
    )))

  }

  output <- stats::setNames(output,paste0("sample_",seq_len(samples)))

  # binding final summary stats/output
  out <- tfb_summary(d,y,estimand,confidence,wdims_samples,v_hats_samples)
  output <- c(out,output)

  return(output)

}
