#' Targeted Function Balancing Optimization with CVXR -- ATE
#'
#' This function solves the linear optimization problem presented by TFB with CVXR, given information from an initial regression, with the ATE as the estimand.
#' @param X_c covariate matrix used to fit f_0
#' @param X_t covariate matrix used to fit f_1
#' @param beta_c coefficient vector for the initial regression on the control set
#' @param beta_t coefficient vector for the initial regression on the test set
#' @param sqrtV_c square root of variance matrix for the initial regression on the control set
#' @param sqrtV_t square root of variance matrix for the initial regression on the test set
#' @param sigma2_c estimated residual variance for the initial regression on the control set
#' @param sigma2_t estimated residual variance for the initial regression on the test set
#' @param i fold indices
#' @param d treatment status vector
#' @param chi_q probability threshold
#' @param quiet whether to suppress console output
#' @param solver solver from CVXR
#' @param rtol tolerance level for CVXR solver
#'
#' @returns A numeric vector of weights.
#'
#' @import CVXR
#' @importFrom stats qchisq
#' @importFrom methods as
#' @keywords tfb
#' @examples
#' set.seed(1221)
#' i <- sample(rep(c(TRUE,FALSE),75))
#' # only observations in the split are used to fit the initial regression
#' X <- as.matrix(iris[, 2:4])[i,]
#' d <- iris[, 5] == "setosa"
#' # coefficients for iris data from tfb_target_ols
#' beta_c <- c(0.565,0.812,-0.686)
#' beta_t <- c(0.532,0.632,-0.251)
#' # root of covariance matrix for iris data from tfb_sqrtV
#' sqrtV_c <- matrix(
#'   c(0.177,-0.017,-0.0143,
#'   -0.017,0.0904,-0.0614,
#'   -0.0143,-0.0614,0.199),nrow=3
#' )
#' sqrtV_t <- matrix(
#'   c(0.158,-0.0511,0.00316,
#'   -0.0511,0.423,-0.126,
#'   0.00316,-0.126,0.614),nrow=3
#' )
#' # standard error for iris data from tfb_target_ols
#' sigma2_c <- 0.12
#' sigma2_t <- 0.0616
#' tfb:::tfb_balance_rmosek_ate(
#'   X,X,beta_c,beta_t,sqrtV_c,sqrtV_t,sigma2_c,sigma2_t,i,d,0.95,TRUE, "CLARABEL", 1e-6
#' )


### Optimization Function
tfb_balance_cvxr_ate <- function(
    X_c,
    X_t,
    beta_c,
    beta_t,
    sqrtV_c,
    sqrtV_t,
    sigma2_c,
    sigma2_t,
    i,
    d,
    chi_q,
    quiet,
    solver,
    rtol

){

  if (quiet) {
    verb <- FALSE
  } else {
    verb <- TRUE
  }
  d <- d[i == 1]

                   # Define values we need
  n_c <- sum(d==0) # control sample size
  n_t <- sum(d==1) # treated sample size
  p_c <- ncol(X_c) # number of covariates in X_c
  p_t <- ncol(X_t) # number of covariates in X_t

  # Define linear function
  c_CVXR <- c(
    rep(0, n_c),               # w_c
    rep(0, n_t),               # w_t
    rep(0, p_c),               # v_rawc
    rep(0, p_t),               # v_rawt
    rep(0, p_c),               # v_tfc
    rep(0, p_t),               # v_tft
    rep(0, 2),                 # t_1c
    rep(0, 2),                 # t_1t
    0,                         # t_2c
    0,                         # t_2t
    0,                         # t_3
    c(1, 0),                   # u_1
    c(sigma2_c / (n_c^2), 0),  # u_2c
    c(sigma2_t / (n_t^2), 0),  # u_2t
    rep(0, 2),                 # tilde u_2c
    rep(0, 2)                  # tilde u_2t
  )                            # final problem: $c^t\theta=u_1^{(1)}+\frac{\hat\sigma_c^2}{n_c^2}u_{2c}^{(1)}+\frac{\hat\sigma_t^2}{n_t^2}u_{2t}^{(1)}$

  # Make CVXR Variables
  dim <- length(c_CVXR) # dimension of theta
  parvec <- Variable(dim)

  # Define Linear constraint matrices
  A.sum_to_n_c <- t(as.matrix(c( # Define Linear constraint matrices A
    rep(1, n_c),                 # w_c
    rep(0, dim - n_c)            # everything else
  )))                            # w_c sum to n_c: $\sum_{i=1}^{n_c}w_c^{(i)}=n_c$
  bc.sum_to_n_c <- matrix(n_c, nrow=2) # $n_c\leq A_{n_c}\cdot\theta\leq n_c$

  A.sum_to_n_t <- t(as.matrix(c(
    rep(0, n_c),                 # w_c
    rep(1, n_t),                 # w_t
    rep(0, dim - n_c - n_t)      # everything else
  )))                            # w_t sum to n_t: $\sum_{i=1}^{n_t}w_t^{(i)}=n_t$
  bc.sum_to_n_t <- matrix(n_t, nrow=2) # $n_t\leq A_{n_t}\cdot\theta\leq n_t$

  A.raw_c <- cbind(
    (1 / n_c) * t(X_c[d==0, ]),                         # w_c
    matrix(0, nrow = p_c, ncol = n_t),                  # w_t
    as.matrix(diag(rep(1, p_c+1))[1:p_c, 1:p_c]),       # v_rawc
    matrix(0, nrow = p_c, ncol = dim - n_c - n_t - p_c) # everything else
  )                                                     # raw imbalance: $V_{1c}=\frac{1}{n}\sum_{i=1}^nX_i-\frac{1}{n_c}\sum_{i=1}^{n_c}w_c^{(i)}X_i$
  bc.raw_c <- matrix(colMeans(X_c), nrow=2, ncol=p_c, byrow=T) # $\frac{1}{n}\sum_{i=1}^nX_i\leq A_{\text{raw}_c}\cdot\theta\leq\frac{1}{n}\sum_{i=1}^nX_i$

  A.raw_t <- cbind(
    matrix(0, nrow = p_t, ncol = n_c),                        # w_c
    (1 / n_t) * t(X_t[d==1, ]),                               # w_t
    matrix(0, nrow = p_t, ncol = p_c),                        # v_rawc
    as.matrix(diag(rep(1, p_t+1))[1:p_t, 1:p_t]),             # v_rawt
    matrix(0, nrow = p_t, ncol = dim - n_c - n_t - p_c - p_t) # everything else
  )                                                           # raw imbalance: $V_{1t}=\frac{1}{n}\sum_{i=1}^nX_i-\frac{1}{n_t}\sum_{i=1}^{n_t}w_t^{(i)}X_i$
  bc.raw_t <- matrix(colMeans(X_t), nrow=2, ncol=p_t, byrow=T) # $\frac{1}{n}\sum_{i=1}^nX_i\leq A_{\text{raw}_t}\cdot\theta\leq\frac{1}{n}\sum_{i=1}^nX_i$

  A.tf_c <- cbind(
    (1 / n_c) * sqrtV_c %*% t(X_c[d==0, ]),             # w_c
    matrix(0, nrow = p_c, ncol = n_t),                  # w_t
    as.matrix(diag(rep(0, p_c+1))[1:p_c, 1:p_c]),       # v_rawc
    matrix(0, nrow = p_c, ncol = p_t),                  # v_rawt
    as.matrix(diag(rep(1, p_c+1))[1:p_c, 1:p_c]),       # v_tfc
                                                        # everything else
    matrix(0, nrow = p_c, ncol = dim - n_c - n_t - 2 * p_c - p_t)
  )                                                     # transformed imbalance: $V_{2c}=\hat V_{\beta_0}^{\frac{1}{2}}left(\frac{1}{n}\sum_{i=1}^nX_i-\frac{1}{n_c}\sum_{i=1}^{n_c}w_c^{(i)}X_i\right)$
  bc.tf_c <- matrix(colMeans(X_c %*% sqrtV_c), nrow=2, ncol=p_c, byrow=T)

  A.tf_t <- cbind(
    matrix(0, nrow = p_t, ncol = n_c),               # w_c
    (1 / n_t) * sqrtV_t %*% t(X_t[d==1, ]),          # w_t
    matrix(0, nrow = p_t, ncol = p_c),               # v_rawc
    as.matrix(diag(rep(0, p_t+1))[1:p_t, 1:p_t]),    # v_rawt
    matrix(0, nrow = p_t, ncol = p_c),               # v_tfc
    as.matrix(diag(rep(1, p_t+1))[1:p_t, 1:p_t]),    # v_tft
                                                     # everything else
    matrix(0, nrow=p_t, ncol= dim - n_c - n_t - 2 * p_c - 2 * p_t)
  )                                                  # transformed imbalance: $V_{2t}=\hat V_{\beta_1}^{\frac{1}{2}}left(\frac{1}{n}\sum_{i=1}^nX_i-\frac{1}{n_t}\sum_{i=1}^{n_c}w_c^{(i)}X_i\right)$
  bc.tf_t <- matrix(colMeans(X_t %*% sqrtV_t), nrow=2, ncol=p_t, byrow=T)

  A.mag_c <- c(
    rep(0, n_c),                        # w_c
    rep(0, n_t),                        # w_t
    beta_c,                             # v_1c
    rep(0, p_t),                        # v_1t
    rep(0, p_c),                        # v_2c
    rep(0, p_t),                        # v_2t
    c(-1, 1),                           # t_1c
                                        # everything else
    rep(0, dim - n_c - n_t - 2 * p_c - 2 * p_t - 2)
  )                                     # magnified imbalance: $\beta_0^TV_{1c}=t_{1c}^{(1)}-t_{1c}^{(2)}$
  A.mag_c <- t(as.matrix(A.mag_c))

  A.mag_t <- c(
    rep(0, n_c),                        # w_c
    rep(0, n_t),                        # w_t
    rep(0, p_c),                        # v_1c
    beta_t,                             # v_1t
    rep(0, p_c),                        # v_2c
    rep(0, p_t),                        # v_2t
    rep(0, 2),                          # t_1c
    c(-1, 1),                           # t_1t
                                        # everything else
    rep(0, dim - n_c - n_t - 2 * p_c - 2 * p_t - 4)
  )                                     # magnified imbalance: $\beta_1^TV_{1t}=t_{1t}^{(1)}-t_{1t}^{(2)}$
  A.mag_t <- t(as.matrix(A.mag_t))

  bc.mag_c <- matrix(0, nrow=2)
  bc.mag_t <- matrix(0, nrow=2)

  A.bias <- c(
    rep(0, n_c),                        # w_c
    rep(0, n_t),                        # w_t
    rep(0, p_c),                        # v_1c
    rep(0, p_t),                        # v_1t
    rep(0, p_c),                        # v_2c
    rep(0, p_t),                        # v_2t
    rep(1, 2),                          # t_1c
    rep(1, 2),                          # t_1t
    sqrt(qchisq(p = chi_q, df=p_c)),    # t_2c
    sqrt(qchisq(p = chi_q, df=p_c)),    # t_2t
    -1,                                 # t_3
                                        # everything else
    rep(0, dim - n_c - n_t - 2 * p_c - 2 * p_t - 7)
  )                                     # bias: $t_3=(t_{1c}^{(1)}+t_{1c}^{(2)})+(t_{1t}^{(1)}+t_{1t}^{(2)})+\sqrt{Q_q\Chi_p^2}(t_{2c}+t_{2t})$
  A.bias <- t(as.matrix(A.bias))
  bc.bias <- matrix(0, nrow=2)

  A.u2s <- cbind(
    matrix(0, nrow=3, ncol=n_c),        # w_c
    matrix(0, nrow=3, ncol=n_t),        # w_t
    matrix(0, nrow=3, ncol=p_c),        # v_1c
    matrix(0, nrow=3, ncol=p_t),        # v_1t
    matrix(0, nrow=3, ncol=p_c),        # v_2c
    matrix(0, nrow=3, ncol=p_t),        # v_2t
    matrix(0, nrow=3, ncol=2),          # t_1c
    matrix(0, nrow=3, ncol=2),          # t_1t
    matrix(0, nrow=3, ncol=1),          # t_2c
    matrix(0, nrow=3, ncol=1),          # t_2t
    matrix(0, nrow=3, ncol=1),          # t_3
    matrix(0, nrow=3, ncol=1),          # u_1^(1)
    matrix(c(1, 0, 0), nrow=3, ncol=1), # u_1^(2)
    matrix(0, nrow=3, ncol=1),          # u_2c^(1)
    matrix(c(0, 1, 0), nrow=3, ncol=1), # u_2c^(2)
    matrix(0, nrow=3, ncol=1),          # u_2t^(1)
    matrix(c(0, 0, 1), nrow=3, ncol=1),  # u_2t^(2)
    matrix(0, nrow=3, ncol=4)          # tilde u2s
  )                                     # u2s: $u_1^{(2)}=\frac{1}{2}$ $u_{2c}^{(2)}=\frac{1}{2}$ $u_{2t}^{(2)}=\frac{1}{2}$
  bc.u2s <- matrix(0.5, nrow=2, ncol=3)

  A.tildeu2s <- cbind(
    matrix(0, nrow=4, ncol=n_c),        # w_c
    matrix(0, nrow=4, ncol=n_t),        # w_t
    matrix(0, nrow=4, ncol=p_c),        # v_1c
    matrix(0, nrow=4, ncol=p_t),        # v_1t
    matrix(0, nrow=4, ncol=p_c),        # v_2c
    matrix(0, nrow=4, ncol=p_t),        # v_2t
    matrix(0, nrow=4, ncol=2),          # t_1c
    matrix(0, nrow=4, ncol=2),          # t_1t
    matrix(0, nrow=4, ncol=1),          # t_2c
    matrix(0, nrow=4, ncol=1),          # t_2t
    matrix(0, nrow=4, ncol=1),          # t_3
    matrix(0, nrow=4, ncol=1),          # u_1^(1)
    matrix(0, nrow=4, ncol=1),          # u_1^(2)
    matrix(c(1/sqrt(2), 1/sqrt(2), 0, 0), nrow=4, ncol=1),          # u_2c^(1)
    matrix(c(-1/sqrt(2), 1/sqrt(2), 0, 0), nrow=4, ncol=1),         # u_2c^(2)
    matrix(c(0, 0, 1/sqrt(2), 1/sqrt(2)), nrow=4, ncol=1),          # u_2t^(1)
    matrix(c(0, 0, -1/sqrt(2), 1/sqrt(2)), nrow=4, ncol=1),         # u_2t^(2)
    matrix(c(-1, 0, 0, 0), nrow=4, ncol=1),  # for tilde u_2c^(1) -- this is the minus one
    matrix(c(0, -1, 0, 0), nrow=4, ncol=1),  # for tilde u_2c^(2) -- this is the plus one
    matrix(c(0, 0, -1, 0), nrow=4, ncol=1),  # for tilde u_2t^(1) -- this is the minus one
    matrix(c(0, 0, 0, -1), nrow=4, ncol=1)  # for tilde u_2t^(2) -- this is the plus one
  )
  bc.tildeu2s <- matrix(0, nrow=2, ncol=4)

  A_CVXR <- rbind(     # linear constraints
    A.sum_to_n_c,
    A.sum_to_n_t,
    A.raw_c,
    A.raw_t,
    A.tf_c,
    A.tf_t,
    A.mag_c,
    A.mag_t,
    A.bias,
    A.u2s,
    A.tildeu2s
  )

  bc_CVXR <- as.matrix(cbind(     # linear constraint bounds
    bc.sum_to_n_c,
    bc.sum_to_n_t,
    bc.raw_c,
    bc.raw_t,
    bc.tf_c,
    bc.tf_t,
    bc.mag_c,
    bc.mag_t,
    bc.bias,
    bc.u2s,
    bc.tildeu2s
  )[1, ])

  # Define CVXR equality constraints
  equality_constraints_CVXR <- A_CVXR %*% parvec == bc_CVXR

  # Define upper and lower bounds for parameters
  lower <- c(       # parameter lower bounds
    rep(0, n_c),    # w_c
    rep(0, n_t),    # w_t
    rep(-Inf, p_c), # v_1c
    rep(-Inf, p_t), # v_1t
    rep(-Inf, p_c), # v_2c
    rep(-Inf, p_t), # v_2t
    rep(0, 2),      # t_1c
    rep(0, 2),      # t_1t
    0,              # t_2c
    0,              # t_2t
    0,              # t_3
    rep(0, 2),      # u_1
    rep(0, 2),      # u_2c
    rep(0, 2),       # u_2t
    c(-0.5 / sqrt(2), 0.5 / sqrt(2), -0.5 / sqrt(2), 0.5 / sqrt(2)) # tilde u2s
  )
  lower_CVXR <- matrix(lower)

  upper <- c(      # parameter upper bounds
    rep(n_c, n_c), # w_c
    rep(n_t, n_t), # w_t
    rep(Inf, p_c), # v_1c
    rep(Inf, p_t), # v_1t
    rep(Inf, p_c), # v_2c
    rep(Inf, p_t), # v_2t
    rep(Inf, 2),   # t_1c
    rep(Inf, 2),   # t_1t
    Inf,           # t_2c
    Inf,           # t_2t
    Inf,           # t_3
    rep(Inf, 2),   # u_1
    rep(Inf, 2),   # u_2c
    rep(Inf, 2),   # u_2t
    rep(Inf, 4)    # tilde u2s --- (n_c+2p+5):(n_c+2p+10)
  )
  upper_CVXR <- matrix(upper)

  # Define bounds on parameters
  lower_finite <- as.vector(!is.infinite(lower_CVXR))
  lowerbound_constraints_CVXR <- diag(1*lower_finite)[lower_finite, ] %*% parvec >= matrix(lower_CVXR[lower_finite, ])

  upper_finite <- as.vector(!is.infinite(upper_CVXR))
  upperbound_constraints_CVXR <- diag(1*upper_finite)[upper_finite, ] %*% parvec <= matrix(upper_CVXR[upper_finite, ])

  # Quadratic Cone constraint
  isolate_v2c <- rep(0, dim)
  isolate_v2c[(n_c + n_t + p_c + p_t + 1):(n_c + n_t + 2 * p_c + p_t)] <- 1
  isolate_v2c <- diag(isolate_v2c)
  find_t2c <- rep(0, dim)
  find_t2c[n_c + n_t + 2 * p_c + 2 * p_t + 5] <- 1
  quadcone1_constraints_CVXR <- p_norm(isolate_v2c %*% parvec, 2) <= t(find_t2c) %*% parvec

  isolate_v2t <- rep(0, dim)
  isolate_v2t[(n_c + n_t + 2 * p_c + p_t + 1):(n_c + n_t + 2 * p_c + 2 * p_t)] <- 1
  isolate_v2t <- diag(isolate_v2t)
  find_t2t <- rep(0, dim)
  find_t2t[n_c + n_t + 2 * p_c + 2 * p_t + 6] <- 1
  quadcone2_constraints_CVXR <- p_norm(isolate_v2t %*% parvec, 2) <= t(find_t2t) %*% parvec

  # Rotated Cone constraint
  isolate_t3 <- rep(0, dim)
  isolate_t3[n_c + n_t + 2 * p_c + 2 * p_t + 7] <- 1
  isolate_t3 <- diag(isolate_t3)
  find_u11 <- rep(0, dim)
  find_u11[n_c + n_t + 2 * p_c + 2 * p_t + 8] <- 1
  rotquadcone1_constraints_CVXR <- sum_squares(isolate_t3 %*% parvec) <= t(find_u11) %*% parvec

  isolate_w_tildeu2c1 <- rep(0, dim)
  isolate_w_tildeu2c1[1:n_c] <- 1
  isolate_w_tildeu2c1[dim-3] <- 1
  isolate_w_tildeu2c1 <- diag(isolate_w_tildeu2c1)
  find_tildeu2c2 <- rep(0, dim)
  find_tildeu2c2[dim-2] <- 1
  rotquadcone2_constraints_CVXR <- p_norm(isolate_w_tildeu2c1 %*% parvec, 2) <= t(find_tildeu2c2) %*% parvec # Note this is a rotated cone constraint that has been recoded as a quadratic cone constraint.

  isolate_w_tildeu2t1 <- rep(0, dim)
  isolate_w_tildeu2t1[(n_c + 1):(n_c + n_t)] <- 1
  isolate_w_tildeu2t1[dim-1] <- 1
  isolate_w_tildeu2t1 <- diag(isolate_w_tildeu2t1)
  find_tildeu2t2 <- rep(0, dim)
  find_tildeu2t2[dim] <- 1
  rotquadcone3_constraints_CVXR <- p_norm(isolate_w_tildeu2t1 %*% parvec, 2) <= t(find_tildeu2t2) %*% parvec # Note this is a rotated cone constraint that has been recoded as a quadratic cone constraint.

  # Get CVXR solution
  prob <- Problem(
    Minimize(t(c_CVXR) %*% parvec),
    constraints = c(
      equality_constraints_CVXR,
      lowerbound_constraints_CVXR,
      upperbound_constraints_CVXR,
      quadcone1_constraints_CVXR,
      quadcone2_constraints_CVXR,
      rotquadcone1_constraints_CVXR,
      rotquadcone2_constraints_CVXR,
      rotquadcone3_constraints_CVXR
    )
  )
  result <- psolve(prob, verbose=verb, feastol=rtol, reltol=rtol, abstol=rtol, solver=solver)

  # Pull out the parameter values
  parameters <- value(parvec)
  w <- rep(0, n_c + n_t)
  w[d == 0] <- parameters[1:n_c]
  w[d == 1] <- parameters[(n_c+1):(n_c+n_t)]

  return(w)
}
