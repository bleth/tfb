#' Targeted Function Balancing
#'
#' @description `tfb` performs Targeted Function Balancing ([Wainstein and Bai, 2025](https://arxiv.org/abs/2203.12179)), a weighting method for causal inference. `tfb` weights observations in your dataset by minimizing targeted function imbalance from an initial regression, and returns a causal estimate as well as summary statistics. `tfb` is intended for use with data from observational studies where inference about the effect of a treatment is desired even though the treatment assignment mechanism is unknown or dependent on the covariates.
#' @param X The covariate dataset, containing all observed confounders of the relationship between the treatment and the outcome. `X` can contain categorical covariates if a dataframe is supplied. A categorical covariate with `n` levels will be converted to `n-1` indicator variables for `tfb`'s optimization. A numeric matrix or a dataframe. No default.
#' @param d The dichotomous treatment vector, which takes the value of `0` or `FALSE` when an observation is in the control group, and the value of `1` or `TRUE` when an observation is in the treatment group. If a dichotomous character vector is supplied, the first value in the vector will be interpreted as indicating treatment. A numeric, logical, or character vector. No default.
#' @param y The outcome vector, containing the measured response for each observation. A numeric vector. No default.
#' @param fit The hypothesized functional form of the relationship between the outcome and the covariates. This decides what regression model will be fit on the data and targeted by the targeted function balancing. Can be set to one or two of the following:
#' * `"ols"`: ordinary least squares ( `lm()` ) will be used.
#' * `"krls"`: kernel regularized least squares ( `krls()` from [`KRLS`](https://cran.r-project.org/web/packages/KRLS/index.html)) will be used.
#' * `"elasticnet"`: an elastic net penalized regression ( `glmnet()` from [`glmnet`](https://cran.r-project.org/web/packages/glmnet/index.html)) will be used.
#' * `"bart"`: a bayesian additive regression tree ( `gbart()` from [`BART`](https://cran.r-project.org/web/packages/BART/index.html)) will be used.
#' * `"custom"`: a custom fit. Requires user to manually specify `beta_c` (or `fhat_c`) and `Vhat_c` for the control group, and/or `beta_t` (or `fhat_t`) and `Vhat_t` for the treated group. Does not allow sample-splitting (i.e., `folds=1` is required.)
#' When estimating the ATT or the ATC, `fit` should only be one value. When estimating the ATE, `fit` may include two different values to allow for the control (first value) and treatment (second value) potential outcomes to be modeled by different families of functions on the covariates. If `fit` is only one value when estimating the ATE, then `tfb` will use the same family of functions to model the treatment and control potential outcomes. A character or character vector. No default.
#' @param estimand The causal inference estimand to be estimated via targeted function balancing. Can be set to one of the following:
#' * `"att"`: the average treatment effect on the treated (ATT) will be estimated.
#' * `"atc"`: the average treatment effect on the controlled (ATC) will be estimated.
#' * `"ate"`: the average treatment effect (ATE) will be estimated.
#'  A character. No default.
#' @param quiet If set to `TRUE`, suppresses console output. A logical. Default is `TRUE`.
#' @param folds The number of folds for sample-splitting, if desired. The data will be partitioned uniformly at random into this many folds. An integer. Default is `1`, which means no sample-splitting.
#' @param samples The number of times to sample the data; the number of times the data is split into folds. An integer. Default is `1`. Must be `1` if `folds=1`.
#' @param bstrap_cov Whether to bootstrap the covariance. If `FALSE`, `tfb` will estimate the covariance matrix via closed formula where applicable. `bstrap_cov` is ignored when there is no closed form for the covariance matrix (e.g., `fit="elasticnet"`). Setting this to `TRUE` may significantly extend calculation times. A logical. Default is `FALSE`.
#' @param ... Additional named arguments. See *Details* for more information.
#'
#' @returns
#'
#' When sample-splitting is not employed (i.e., `folds=1`), the output is a nested list that contains:
#' * `estimate`: The TFB estimate/weighted difference in means.
#' * `variance`: The variance of the estimate.
#' * `weights`: TFB's weights
#' * `confidence_interval`: The confidence interval for the estimate.
#' * `p_value`: The p-value for the estimate under the null hypothesis that there is no treatment effect.
#' * `betas`: The covariate coefficients for the targeted function(s).
#' * `Xs`: The covariate matrices for the targeted function(s).
#' * `covariances`: The covariance matrices for the targeted function(s).
#' * `yhats`: The predicted values from the targeted function(s).
#'
#' When sample-splitting is employed (i.e., `folds>1`), the output is a nested list that contains:
#' * `final`: A list of summary statistics and estimates for the dataset.
#' * `sample_i`: A list of summary statistics for the ith sample-split iteration.
#'
#' The `final` list contains:
#' * `estimate`: The TFB estimate/weighted difference in means. This is the median of all of the estimates from each sample-split iteration.
#' * `variance`: The variance of the above `estimate`.
#' * `confidence_interval`: The confidence interval for the `estimate`.
#' * `p_value`: The p-value for the `estimate` under the null hypothesis that there is no treatment effect.
#'
#' Each `sample_i` list contains:
#' * `final`: A list of summary statistics and estimates for the sample-split iteration. The estimate here is the average of the estimates across the folds for the sample-split iteration.
#'
#' as well various results for each fold within the sample-split:
#' * `indices`: The fold indices.
#' * `betas`: The covariate coefficients for the targeted function(s).
#' * `Xs`: The covariate matrices for the targeted function(s).
#' * `covariances`: The covariance matrices for the targeted function(s).
#' * `yhats`: The predicted values from the targeted function(s).
#' * `estimates`: The TFB estimates/weighted differences in means(s).
#' * `variances`: The variances of the above `estimates`(s).
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
#' * `X_c`, `X_t`: A way to modify training data for \eqn{\hat f_0} and \eqn{\hat f_1} _separately_ of the rest of the dataset. For example, you could fit an \eqn{\hat f_0 = \beta_0+\beta_1x_1+\beta_2x_2^2} by letting `X_c = cbind(X[, 1],X[, 2]^2)` while leaving \eqn{\hat f_1 = \beta_0 + \beta_1x_1+\beta_2x_2}. Numeric matrices. Default is `X`.
#' * `beta_c`, `beta_t`: User-chosen \eqn{\hat \beta_0} and \eqn{\hat \beta_1} coefficient vectors for the control and treatment group, respectively. In this case, \eqn{\hat f_0 (X) = X_c \hat \beta_0} and \eqn{\hat f_1 (X) = X_t \hat \beta_1}. `beta_c`/`beta_t` (or `fhat_c`/`fhat_t`) must be specified for the control/treatment group, respectively, if the corresponding `fit = "custom"`. Can only be used if not sample splitting (i.e., `folds=1`). Numeric vectors. If the length of `beta_c`/`beta_t` is the same as the number of columns of `X_c`/`X_t` (after removing collinear and constant columns), then no intercept is assumed. If the length of `beta_c`/`beta_t` is the same as one plus the number of columns of `X_c`/`X_t` (after removing colliner and constant columns), then the first element of `beta_c`/`beta_t` is assumed to be the intercept. Default is `NULL`.
#' * `fhat_c`, `fhat_t`: User-chosen \eqn{\hat f_0} and \eqn{\hat f_1} predicted value vectors for the control and treatment group, respectively. If `beta_c`/`beta_t` is not specified, then `fhat_c`/`fhat_t` must be specified for the control/treatment group, respectively, if the corresponding `fit = "custom"`. Numeric vectors. Can only be used if not sample splitting (i.e., `folds=1`). Default is `NULL`.
#' * `Vhat_c`, `Vhat_t`: User-chosen variance matrices for `beta_c` and `beta_t` (or `fhat_c` and `fhat_t`), respectively. `Vhat_c`/`Vhat_t` must be specified for the control/treatment group, respectively, if the corresponding `fit = "custom"`. Numeric matrices. Can only be used if not sample splitting (i.e., `folds=1`). Default is `NULL`.
#'
#'
#' @import fastDummies
#' @importFrom stats setNames
#' @keywords tfb
#' @export
#' @examples
#'
#' X <- iris[, 2:4]
#' d <- iris[, 5] == "setosa"
#' y <- iris[, 1]
#' # estimating the effect of `Species == "setosa"` on `Sepal.Length`
#' # given `Sepal.Width`, `Petal.Length`, and `Petal.Width`
#' tfb(X, d, y, fit = "ols", estimand = "att")$estimate
#'

tfb <- function(
    X,              # covariate matrix
    d,              # treatment vector
    y,              # response vector
    fit,            # c("ols", "krls", "elasticnet", "bart")
    estimand,       # c("atc", "att", "ate")
    quiet = T,      # suppress console output
    folds = 1,      # number of folds for sample-splitting
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
    "X_t",
    "beta_c",
    "beta_t",
    "fhat_c",
    "fhat_t",
    "Vhat_c",
    "Vhat_t"
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
    X <- fastDummies::dummy_cols(X, remove_first_dummy = T, remove_selected_columns = T)
  }
  X <- as.matrix(X)
  if (sum(apply(X, 2, var)==0)>0) {
    message("Constant columns removed from X.\n")
    to_remove <- which(apply(X, 2, var)==0)
    X <- as.matrix(X[, -to_remove])
  }
  if (nrow(X) < 1 | ncol(X) < 1) {stop("`X` must contain data.\n")}
  if (!is.vector(d) | !(is.numeric(d) | is.logical(d) | is.character(d))) {stop("`d` must be a numeric, logical, or character vector.\n")}
  if (length(unique(d)) < 2) {stop("`d` cannot be a degenerate random variable.\n")}
  if (length(unique(d)) > 2) {stop("`d` cannot take more than 2 values, as `tfb` has not been generalized to non-dichotomous treatments at this time.\n")}
  if (is.character(d)) {d <- d == unique(d)[1]}
  if (length(d) != nrow(X)) {stop("`d` must be the same length as `X`.\n")}
  if (!is.vector(y) | !(is.numeric(y) | is.logical(y))) {stop("`y` must be a numeric or logical vector.\n")}
  if (sum(is.na(lm(y ~ X)$coefficients[-1]))>0) {
    message("Collinear columns removed from X.\n")
    to_remove <- which(is.na(lm(y ~ X)$coefficients[-1]))
    X <- as.matrix(X[, -to_remove])
  }
  if (length(y) != nrow(X)) {stop("`y` must be the same length as `X`.\n")}
  if (!(estimand %in% c("att","atc","ate"))) {stop("`estimand` must be one of `\"att\",\"atc\",\"ate\".\n")}
  if (any(!(fit %in% c("ols","krls","elasticnet","bart", "custom")))) {stop("`fit` must be one of `\"ols\",\"krls\",\"elasticnet\",\"bart\",\"custom\"`.\n")}
  if (length(fit) > 1 & estimand != "ate") {warning("Multiple fits is only intended for `estimand=\"ate\". Using only the first fit specified.\n")}
  if (length(fit) > 2 & estimand == "ate") {warning("`estimand=\"ate\" uses at most two different fits. Using only the first two fits specified.\n")}
  if (length(fit) == 1 & estimand == "ate") {fit <- rep(fit, 2)}
  if (!is.logical(quiet)) {stop("`quiet` must be a logical.\n")}
  if (!is.wholenumber(folds)) {stop("The number of `folds` must be an integer.\n")}
  if (folds < 1 | folds > nrow(X)) {stop("The number of `folds` must be at least `1` and at most `n`.\n")}
  if (any(fit %in% c("custom")) & folds!=1) {stop("Cross-fitting is not allowed (i.e., `folds=1` is required) if one of the `fit = \"custom\"`.\n")}
  if (!is.wholenumber(samples)) {stop("The number of `samples` must be an integer.\n")}
  if (samples < 1) {stop("The number of `samples` must be at least `1`.\n")}
  if (folds==1 & samples>1) {stop("If not sample-splitting (i.e., `folds=1`) then `samples` must be `1`.\n")}
  if (!is.logical(bstrap_cov)) {stop("`bstrap_cov` must be a logical.\n")}
  if (all(fit %in% "elasticnet") & bstrap_cov == F) {
    if (quiet == F) {warning("Ignoring user input of `bstrap_cov = F`. The covariance must be bootstrapped when `fit = \"elasticnet\"`.\n")}
  }
  if (all(fit %in% "bart") & bstrap_cov == T) {
    if (quiet == F) {warning("Ignoring user input of `bstrap_cov = T`. The covariance must be the covariance of posterior distribution when `fit = \"bart\"`.\n")}
  }

  # check extended argument quality
  if (!exists("bstrap_reps", inherits = FALSE)) {
    bstrap_reps <- 400
  } else {
    if (!is.wholenumber(bstrap_reps)) {stop("`bstrap_reps` must be an integer.\n")}
    if (bstrap_reps < 1) {stop("`bstrap_reps` must be at least `1`.\n")}
  }

  if (!exists("confidence", inherits = FALSE)) {
    confidence <- 0.95
  } else {
    if (!is.double(confidence)) {stop("`confidence` must be a proportion.\n")}
    if (confidence <= 0 | confidence >= 1) {stop("`confidence` must be a proportion.\n")}
  }

  if (!exists("chi_q", inherits = FALSE)) {
    chi_q <- 0.95
  } else {
    if (!is.double(chi_q)) {stop("`chi_q` must be a proportion.\n")}
    if (chi_q <= 0 | chi_q >= 1) {stop("`chi_q` must be a proportion.\n")}
  }

  if (!exists("ev_approx", inherits = FALSE)) {
    ev_approx <- F
  } else {
    if (!is.logical(ev_approx)) {stop("`ev_approx` must be a logical.\n")}
    if ("krls" %in% fit & ev_approx) {warning("Use of `ev_approx` with KRLS is not recommended.\n")}
  }

  if (!exists("ev_n", inherits = FALSE)) {
    ev_n <- NULL
  } else {
    if (!is.wholenumber(ev_n)) {stop("`ev_n` must be an integer.\n")}
    if (!("krls" %in% fit) & (ev_n < 1 | ev_n > ncol(X))) {stop(paste0("For `fit` of ", fit, " `ev_n` must be at least 1 and at most ", ncol(X), "\n."))}
    if ("krls" %in% fit & (ev_n < 1)) {stop(paste0("For `fit` of ", fit, " `ev_n` must be at least 1.\n"))}
    if ("krls" %in% fit & !is.null(ev_n)) {warning("The set value of `ev_n` may equal or exceed the maximum dimensions of the eigenmatrix. In the case that this happens, `tfb` will perform a complete eigenvalue decomposition.\n")}
  }

  if (!exists("kernelfn", inherits = FALSE)) {
    kernelfn <- tfb_kernelfn
  } else {
    if (!is.function(kernelfn)) {stop("`kernelfn` must be a function.\n")}
    if (length(formals(kernelfn)) != 2) {stop("`kernelfn` must be an inner product, taking two rows and returning one value.\n")}
  }

  if (!exists("reg_d", inherits = FALSE)) {
    reg_d <- F
  } else {
    if (!is.logical(reg_d)) {stop("`reg_d` must be a logical.\n")}
  }

  if (!exists("X_c", inherits = FALSE)) {
    X_c <- as.matrix(X)
  } else {
    if (!is.numeric(X_c) | !is.matrix(X_c)) {stop("`X_c` must be a numeric matrix.\n")}
    if (nrow(X_c) != nrow(X)) {stop("`X_c` must have as many rows as `X`.\n")}
    if (sum(apply(X_c, 2, var)==0)>0) {
      message("Constant columns removed from `X_c`.\n")
      to_remove <- which(apply(X_c, 2, var)==0)
      X_c <- as.matrix(X_c[, -to_remove])
    }
    if (sum(is.na(lm(y ~ X_c)$coefficients[-1]))>0) {
      message("Collinear columns removed from `X_c`.\n")
      to_remove <- which(is.na(lm(y ~ X_c)$coefficients[-1]))
      X_c <- as.matrix(X_c[, -to_remove])
    }
  }

  if (!exists("X_t", inherits = FALSE)) {
    X_t <- as.matrix(X)
  } else {
    if (!is.numeric(X_t) | !is.matrix(X_t)) {stop("`X_t` must be a numeric matrix.\n")}
    if (nrow(X_t) != nrow(X)) {stop("`X_t` must have as many rows as `X`.\n")}
    if (sum(apply(X_t, 2, var)==0)>0) {
      message("Constant columns removed from `X_t`.\n")
      to_remove <- which(apply(X_t, 2, var)==0)
      X_t <- as.matrix(X_t[, -to_remove])
    }
    if (sum(is.na(lm(y ~ X_t)$coefficients[-1]))>0) {
      message("Collinear columns removed from `X_t`.\n")
      to_remove <- which(is.na(lm(y ~ X_t)$coefficients[-1]))
      X_t <- as.matrix(X_t[, -to_remove])
    }
  }

  if (!exists("beta_c", inherits = FALSE)) {
    beta_c <- NULL
    intercept_c <- NULL
  } else {
    if (!is.numeric(beta_c)) {stop("`beta_c` must be a numeric vector.\n")}
    if (length(beta_c) != ncol(X_c) & length(beta_c) != (ncol(X_c) + 1)){
      stop("`beta_c` must have as many components as columns of `X` (plus 1, if including intercept term).\n")
    }
  }
  if(!is.null(beta_c) & length(beta_c) == ncol(X_c)){
   intercept_c <- FALSE
   message("No intercept assumed in `beta_c`.")
  }
  if(!is.null(beta_c) & length(beta_c) == (ncol(X_c) + 1)){
    intercept_c <- TRUE
    message("First element of `beta_c` assumed to be intercept term.")
  }

  if (!exists("beta_t", inherits = FALSE)) {
    beta_t <- NULL
    intercept_t <- NULL
  } else {
    if (!is.numeric(beta_t)) {stop("`beta_t` must be a numeric vector.\n")}
    if (length(beta_t) != ncol(X_t) & length(beta_t) != (ncol(X_t) + 1)) {
      stop("`beta_t` must have as many components as columns of `X` (plus 1, if including intercept term).\n")
    }
  }
  if(!is.null(beta_t) & length(beta_t) == ncol(X_t)){
    intercept_t <- FALSE
    message("No intercept assumed in `beta_t`.")
  }
  if(!is.null(beta_t) & length(beta_t) == (ncol(X_t) + 1)){
    intercept_t <- TRUE
    message("First element of `beta_t` assumed to be intercept term.")
  }

  if (!exists("fhat_c", inherits = FALSE)) {
    fhat_c <- NULL
  } else {
    if (!is.numeric(fhat_c)) {stop("`fhat_c` must be a numeric vector.\n")}
    if (length(fhat_c) != nrow(X)) {stop("`fhat_c` must have as many components as rows of `X`.\n")}
  }
  if(!is.null(beta_c) & !is.null(fhat_c)){message("`beta_c` and `fhat_c` both specified. Defaulting to `fhat_c`.\n")}
  if(estimand=="att"){
    if(fit=="custom" & is.null(beta_c) & is.null(fhat_c)){stop("`beta_c` or `fhat_c` must be specified if `fit`=\"custom\" for the control group.\n")}
  }
  if(estimand=="ate"){
    if(fit[1]=="custom" & is.null(beta_c) & is.null(fhat_c)){stop("`beta_c` or `fhat_c` must be specified if `fit`=\"custom\" for the control group.\n")}
  }

  if (!exists("fhat_t", inherits = FALSE)) {
    fhat_t <- NULL
  } else {
    if (!is.numeric(fhat_t)) {stop("`fhat_t` must be a numeric vector.\n")}
    if (length(fhat_t) != nrow(X)) {stop("`fhat_t` must have as many components as rows of `X`.\n")}
  }
  if(!is.null(beta_t) & !is.null(fhat_t)){message("`beta_t` and `fhat_t` both specified. Defaulting to `fhat_t`.\n")}
  if(estimand=="atc"){
    if(fit=="custom" & is.null(beta_t) & is.null(fhat_t)){stop("`beta_t` or `fhat_t` must be specified if `fit`=\"custom\" for the treated group.\n")}
  }
  if(estimand=="ate"){
    if(fit[2]=="custom" & is.null(beta_t) & is.null(fhat_t)){stop("`beta_t` or `fhat_t` must be specified if `fit`=\"custom\" for the treated group.\n")}
  }

  if (!exists("Vhat_c", inherits = FALSE)) {
    Vhat_c <- NULL
  } else {
    if (!is.numeric(Vhat_c) | !is.matrix(Vhat_c)) {stop("`Vhat_c` must be a numeric matrix.\n")}
    if (nrow(Vhat_c) != ncol(Vhat_c)) {stop("`Vhat_c` must be a square matrix.\n")}
    if (!is.null(fhat_c) & ncol(Vhat_c) != length(fhat_c)) {stop("`Vhat_c` must have the same number of rows/columns as the length of `fhat_c`.\n")}
    if (is.null(fhat_c) & !is.null(beta_c) & ncol(Vhat_c) != length(beta_c)) {stop("`Vhat_c` must have the same number of rows/columns as the length of `beta_c`.\n")}
  }
  if (!is.null(fhat_c) & is.null(Vhat_c)){stop("`Vhat_c` must be specified if `fhat_c` is specified.\n")}
  if (is.null(fhat_c) & !is.null(beta_c) & is.null(Vhat_c)){stop("`Vhat_c` must be specified if `beta_c` is specified.\n")}
  if(estimand=="att"){
    if(fit=="custom" & is.null(Vhat_c)){stop("`Vhat_c` must be specified if `fit`=\"custom\" for the control group.\n")}
  }
  if(estimand=="ate"){
    if(fit[1]=="custom" & is.null(Vhat_c)){stop("`Vhat_c` must be specified if `fit`=\"custom\" for the control group.\n")}
  }

  if (!exists("Vhat_t", inherits = FALSE)) {
    Vhat_t <- NULL
  } else {
    if (!is.numeric(Vhat_t) | !is.matrix(Vhat_t)) {stop("`Vhat_t` must be a numeric matrix.\n")}
    if (nrow(Vhat_t) != ncol(Vhat_t)) {stop("`Vhat_t` must be a square matrix.\n")}
    if (!is.null(fhat_t) & ncol(Vhat_t) != length(fhat_t)) {stop("`Vhat_t` must have the same number of rows/columns as the length of `fhat_t`.\n")}
    if (is.null(fhat_t) & !is.null(beta_t) & ncol(Vhat_t) != length(beta_t)) {stop("`Vhat_t` must have the same number of rows/columns as the length of `beta_t`.\n")}
  }
  if (!is.null(fhat_t) & is.null(Vhat_t)){stop("`Vhat_t` must be specified if `fhat_t` is specified.\n")}
  if (is.null(fhat_t) & !is.null(beta_t) & is.null(Vhat_t)){stop("`Vhat_t` must be specified if `beta_t` is specified.\n")}
  if(estimand=="atc"){
    if(fit=="custom" & is.null(Vhat_t)){stop("`Vhat_t` must be specified if `fit`=\"custom\" for the treated group.\n")}
  }
  if(estimand=="ate"){
    if(fit[1]=="custom" & is.null(Vhat_t)){stop("`Vhat_t` must be specified if `fit`=\"custom\" for the treated group.\n")}
  }

  # Reversing things for ATC -- willl need to take this into account for naming later.
  if(estimand == "atc"){
    d <- 1-d # inverting treatment
    beta_c <- beta_t # inverting betas
    beta_t <- NULL
    fhat_c <- fhat_t # inverting fhats
    fhat_t <- NULL
    Vhat_c <- Vhat_t # inverting Vhats
    Vhat_t <- NULL
    intercept_c <- intercept_t # inverting Vhats
    intercept_t <- NULL
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
        if(fit[(treatment == "t") + 1] != "custom"){
          args <- c(paste0("X_",treatment),"d","y",paste0("i_",fold),"bstrap_cov","bstrap_reps","reg_d","treatment","quiet")
          if (fit[(treatment == "t") + 1] == "krls") {args <- c(args,"kernelfn")}
          out <- do.call(paste0("tfb_target_",fit[(treatment == "t") + 1]),c(lapply(args,as.symbol),list(params = list(...)[args_ellipsis %in% args_params])))
          args <- paste0(c("beta_","V_","e_","sigma2_","yhat_","X_"),treatment,fold)
          for (i in seq_along(out)) {assign(args[i],out[[i]])}
        }
        if(fit[(treatment == "t") + 1] == "custom"){
          # Get indices for the fold
          temp_i_in <- get(paste0("i_",fold))
          if (all(temp_i_in)) {
            temp_i_out <- temp_i_in
          } else {
            temp_i_out <- !temp_i_in
          }

          # First, grab the beta and fhat vectors
          temp_beta <- get(paste0("beta_", treatment))
          temp_fhat <- get(paste0("fhat_", treatment))

          # If coefficient is supplied instead of fhat vector
          if(!is.null(temp_beta) & is.null(temp_fhat)){
            # Get beta
            temp_intercept <- get(paste0("intercept_", treatment))
            assign(
              paste0("beta_",treatment,fold),
              temp_beta[(1 + 1*temp_intercept):length(temp_beta)]
            )

            # Get V
            temp_V <- get(paste0("Vhat_", treatment))
            if(temp_intercept==TRUE) temp_V <- temp_V[2:nrow(temp_V), 2:ncol(temp_V)]
            assign(
              paste0("V_",treatment,fold),
              temp_V
            )

            # Get X
            temp_X <- get(paste0("X_", treatment))
            assign(
              paste0("X_",treatment,fold),
              temp_X[temp_i_in, ]
            )

            # Get yhat
            temp <- temp_X
            if(temp_intercept==TRUE) temp <- cbind(1, temp_X)
            temp_yhat <- as.vector(temp %*% temp_beta)
            assign(
              paste0("yhat_",treatment,fold),
              temp_yhat[temp_i_in]
            )
          }

          # If fhat vector is supplied
          if(!is.null(temp_fhat)){
            # Get beta
            assign(
              paste0("beta_",treatment,fold),
              temp_fhat[temp_i_in]
            )

            # Get V
            temp_V <- get(paste0("Vhat_", treatment))
            assign(
              paste0("V_",treatment,fold),
              temp_V[temp_i_in, temp_i_in]
            )

            # Get X
            temp_X <- diag(length(temp_fhat))
            assign(
              paste0("X_",treatment,fold),
              temp_X[temp_i_in, temp_i_in]
            )

            # Get yhat
            temp_yhat <- temp_fhat
            assign(
              paste0("yhat_",treatment,fold),
              temp_yhat[temp_i_in]
            )
          }

          # Get e
          if(treatment=="c") temp_e <- y[d==0 & temp_i_in] - temp_yhat[d==0 & temp_i_in]
          if(treatment=="t") temp_e <- y[d==1 & temp_i_in] - temp_yhat[d==1 & temp_i_in]
          assign(
            paste0("e_",treatment,fold),
            temp_e
          )

          # sigma2
          if(treatment=="c") temp <- y[d==0 & temp_i_out] - temp_yhat[d==0 & temp_i_out]
          if(treatment=="t") temp <- y[d==1 & temp_i_out] - temp_yhat[d==1 & temp_i_out]
          temp_sigma2 <- mean(temp^2)
          assign(
            paste0("sigma2_",treatment,fold),
            temp_sigma2
          )
        }

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
      if (estimand_type == "att") {args <- c(args,"y","estimand")}
      out <- do.call(paste0("tfb_variance_", estimand_type), lapply(args,as.symbol))
      assign(paste0("v_hat_",fold),out)
    }

    # calculating summary stats for the sample
    args <- ls()
    to_remove <- c("beta_c", "beta_t", "fhat_c", "fhat_t", "Vhat_c", "Vhat_t", "X_c", "X_t")
    to_remove <- which(args %in% to_remove)
    args <- args[-to_remove]
    args_to_loop <- args
    if(estimand=="atc"){ # need to correct names for ATC output
      for(nm in args_to_loop){
      for(fold in seq_len(folds)){
        if(endsWith(nm, paste0("_c", fold))){
          temp <- get(nm)
          rm(list=c(nm))
          new_nm <- gsub(paste0("_c", fold), paste0("_t", fold), nm)
          assign(new_nm, temp)
          if(new_nm %in% args) args <- args[-match(new_nm, args)]
          args[args==nm] <- new_nm
        }
      }
      }
    }

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
      indices = setNames(lapply(lapply(args[startsWith(args,"i_")], as.symbol),eval,envir=environment()),args[startsWith(args,"i_")]),
      betas = setNames(lapply(lapply(args[startsWith(args,"beta_")], as.symbol),eval,envir=environment()),args[startsWith(args,"beta_")]),
      Xs = setNames(lapply(lapply(args[startsWith(args,"X_")], as.symbol),eval,envir=environment()),args[startsWith(args,"X_")]),
      covariances = setNames(lapply(lapply(args[startsWith(args,"V_")], as.symbol),eval,envir=environment()),args[startsWith(args,"V_")]),
      yhats = setNames(lapply(lapply(args[startsWith(args,"yhat_")], as.symbol),eval,envir=environment()),args[startsWith(args,"yhat_")]),
      # errors = setNames(lapply(lapply(args[startsWith(args,"e_")], as.symbol),eval,envir=environment()),args[startsWith(args,"e_")]), # removing these -- not needed
      # standard_errors = setNames(lapply(lapply(args[startsWith(args,"sigma2_")], as.symbol),eval,envir=environment()),args[startsWith(args,"sigma2_")]), # removing these -- not needed
      # dims = setNames(lapply(lapply(args[startsWith(args,"dim_")], as.symbol),eval,envir=environment()),args[startsWith(args,"dim_")]), # removing these -- not needed
      estimates = setNames(lapply(lapply(args[startsWith(args,"wdim_")], as.symbol),eval,envir=environment()),args[startsWith(args,"wdim_")]),
      variances = setNames(lapply(lapply(args[startsWith(args,"v_hat_")], as.symbol),eval,envir=environment()),args[startsWith(args,"v_hat_")])
    )))

  }

  output <- setNames(output,paste0("sample_",seq_len(samples)))

  # binding final summary stats/output
  out <- tfb_summary(d,y,estimand,confidence,wdims_samples,v_hats_samples)
  output <- c(out,output)

  # Simplify output when not sample splitting
  if(folds==1 & samples==1){
    temp <- output$sample_1$final
    temp$betas <- output$sample_1$betas
    names(temp$betas) <- gsub("_c1", "_c", names(temp$betas))
    names(temp$betas) <- gsub("_t1", "_t", names(temp$betas))
    temp$Xs <- output$sample_1$Xs
    names(temp$Xs) <- gsub("_c1", "_c", names(temp$Xs))
    names(temp$Xs) <- gsub("_t1", "_t", names(temp$Xs))
    temp$covariances <- output$sample_1$covariances
    names(temp$covariances) <- gsub("_c1", "_c", names(temp$covariances))
    names(temp$covariances) <- gsub("_t1", "_t", names(temp$covariances))
    temp$yhats <- output$sample_1$yhats
    names(temp$yhats) <- gsub("_c1", "_c", names(temp$yhats))
    names(temp$yhats) <- gsub("_t1", "_t", names(temp$yhats))
    output <- temp
  }

  return(output)

}
