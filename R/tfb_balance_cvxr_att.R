#' Targeted Function Balancing Optimization with CVXR -- ATT
#'
#' This function solves the linear optimization problem presented by TFB with CVXR, given information from an initial regression, with the ATT as the estimand.
#' @param X covariate matrix
#' @param beta coefficient vector
#' @param sqrtV square root of variance matrix
#' @param sigma2 estimated variance
#' @param i fold indices
#' @param d treatment status vector
#' @param chi_q probability threshold
#' @param quiet whether to suppress console output
#' @param solver solver from CVXR
#' @param rtol tolerance level for CVXR solver
#' @param maxit maximum number of iterations for CVXR solver
#'
#' @returns A numeric vector of weights.
#'
#' @import CVXR
#' @import ECOSolveR
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
#' beta <- c(0.565,0.812,-0.686)
#' # root of covariance matrix for iris data from tfb_sqrtV
#' sqrtV <- matrix(
#'   c(0.177,-0.017,-0.0143,
#'   -0.017,0.0904,-0.0614,
#'   -0.0143,-0.0614,0.199),nrow=3
#' )
#' # standard error for iris data from tfb_target_ols
#' sigma2 <- 0.12
#' tfb:::tfb_balance_cvxr_att(X,beta,sqrtV,sigma2,i,d,0.95,TRUE, "CLARABEL", 1e-6, 1e3)


### Optimization Function
tfb_balance_cvxr_att <- function(

    X,
    beta,
    sqrtV,
    sigma2,
    i,
    d,
    chi_q,
    quiet,
    solver,
    rtol,
    maxit

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
    X_c <- as.matrix(X[d==0, ]) # control X
    X_t <- as.matrix(X[d==1, ]) # treated X
    p <- ncol(X)     # number of covariates

    # Define linear function
    c_CVXR <- c(
      rep(0, n_c),           # w_c
      rep(0, p),             # v_1c
      rep(0, p),             # v_2c
      rep(0, 2),             # t_1c
      0,                     # t_2c
      0,                     # t_3
      c(1, 0),               # u_1
      c(sigma2 / (n_c^2), 0), # u_2c
      rep(0, 2)               # tilde u_2c
    )

    # Make CVXR Variables
    length_parvec <- n_c + p + p + 2 + 1 + 1 + 2 + 2 + 2
    parvec <- Variable(length_parvec)

    # Define Linear constraint matrices
    A.sum_to_n_c <- c(
      rep(1, n_c),              # for weights
      rep(0, length(c_CVXR) - n_c) # everything else
    )
    A.sum_to_n_c <- t(as.matrix(A.sum_to_n_c)) # Just turning it into a matrix of the right dimension
    bc.sum_to_n_c <- matrix(n_c, nrow=2)         # upper and lower-bound on the constraints
    ### Here we'll want       LB <= A.sum_to_n_c * theta <= UB

    A.balance1 <- cbind(
      (1 / n_c) * t(X_c),                               # for weights
      as.matrix(diag(rep(1, p+1))[1:p, 1:p]),           # for v1
      as.matrix(diag(rep(0, p+1))[1:p, 1:p]),           # for v2
      matrix(0, nrow=p, ncol=4),                        # for t's
      matrix(0, nrow=p, ncol=6)                         # for u's
    )
    bc.balance1 <- matrix(colMeans(X_t), nrow=2, ncol=p, byrow=T)

    A.balance2 <- cbind(
      (1 / n_c) * sqrtV %*% t(X_c),                        # for weights
      as.matrix(diag(rep(0, p+1))[1:p, 1:p]),              # for v1
      as.matrix(diag(rep(1, p+1))[1:p, 1:p]),              # for v2
      matrix(0, nrow=p, ncol=4),                           # for t's
      matrix(0, nrow=p, ncol=6)                            # for u's
    )
    bc.balance2 <- matrix(colMeans(X_t %*% sqrtV), nrow=2, ncol=p, byrow=T)

    A.mag_imbal <- c(
      rep(0, n_c), # for w
      beta,       # for v1
      rep(0, p),  # for v2
      c(-1, 1),   # for t_1's
      rep(0, 2),  # for other t's
      rep(0, 6)   # for u's
    )
    A.mag_imbal <- t(as.matrix(A.mag_imbal))
    bc.mag_imbal <- matrix(0, nrow=2)

    A.bias <- c(
      rep(0, n_c),                # for w
      rep(0, p),                 # for v1,
      rep(0, p),                 # for v2,
      1,                         # t_1^(1)
      1,                         # t_1^(2)
      sqrt(qchisq(p = chi_q, df=p)), # t_2
      -1,                        # t_3
      rep(0, 6)                  # u's
    )
    A.bias <- t(as.matrix(A.bias))
    bc.bias <- matrix(0, nrow=2)

    A.u2s <- cbind(
      matrix(0, nrow=2, ncol=n_c),      # for w
      matrix(0, nrow=2, ncol=p),       # for v1
      matrix(0, nrow=2, ncol=p),       # for v2
      matrix(0, nrow=2, ncol=4),       # for t's
      matrix(0, nrow=2, ncol=1),       # for u_1^(1)
      matrix(c(1, 0), nrow=2, ncol=1), # for u_1^(2)
      matrix(0, nrow=2, ncol=1),       # for u_2^(1)
      matrix(c(0, 1), nrow=2, ncol=1),  # for u_2^(2)
      matrix(0, nrow=2, ncol=2)       # for tilde us
    )
    bc.u2s <- matrix(0.5, nrow=2, ncol=2)

    A.tildeu2s <- cbind(
      matrix(0, nrow=2, ncol=n_c),      # for w
      matrix(0, nrow=2, ncol=p),       # for v1
      matrix(0, nrow=2, ncol=p),       # for v2
      matrix(0, nrow=2, ncol=4),       # for t's
      matrix(0, nrow=2, ncol=2),       # for u_1s
      matrix(c(1/sqrt(2), 1/sqrt(2)), nrow=2, ncol=1),  # for u_2^(1)
      matrix(c(-1/sqrt(2), 1/sqrt(2)), nrow=2, ncol=1),  # for u_2^(2)
      matrix(c(-1, 0), nrow=2, ncol=1),  # for tilde u_2^(1) -- this is the minus one
      matrix(c(0, -1), nrow=2, ncol=1)  # for tilde u_2^(2) -- this is the plus one
    )
    bc.tildeu2s <- matrix(0, nrow=2, ncol=2)

    A_CVXR <- rbind(
      A.sum_to_n_c,
      A.balance1,
      A.balance2,
      A.mag_imbal,
      A.bias,
      A.u2s,
      A.tildeu2s
    )

    bc_CVXR <- as.matrix(cbind(
      bc.sum_to_n_c,
      bc.balance1,
      bc.balance2,
      bc.mag_imbal,
      bc.bias,
      bc.u2s,
      bc.tildeu2s
    )[1, ])

    # Define CVXR equality constraints
    equality_constraints_CVXR <- A_CVXR %*% parvec == bc_CVXR

    # Define upper and lower bounds for parameters
    lower <- c(
      rep(0, n_c),   # weights --- 1:n_c
      rep(-Inf, p), #      v1 --- (n_c+1):(n_c+p)
      rep(-Inf, p), #      v2 --- (n_c+p+1):(n_c+2p)
      rep(0, 4),    #     t's --- (n_c+2p+1):(n_c+2p+4)
      rep(0, 4),     #      regular us --- (n_c+2p+5):(n_c+2p+8)
      c(-0.5 / sqrt(2), 0.5 / sqrt(2))     #      tildes u2s --- (n_c+2p+9):(n_c+2p+10)
    )
    lower_CVXR <- matrix(lower)

    upper <- c(
      rep(n_c, n_c),  # weights --- 1:n_c
      rep(Inf, p), #      v1 --- (n_c+1):(n_c+p)
      rep(Inf, p), #      v2 --- (n_c+p+1):(n_c+2p)
      rep(Inf, 4), #     t's --- (n_c+2p+1):(n_c+2p+4)
      rep(Inf, 6)  #      us --- (n_c+2p+5):(n_c+2p+10)
    )
    upper_CVXR <- matrix(upper)

    # Define bounds on parameters
    lower_finite <- as.vector(!is.infinite(lower_CVXR))
    lowerbound_constraints_CVXR <- diag(1*lower_finite)[lower_finite, ] %*% parvec >= matrix(lower_CVXR[lower_finite, ])

    upper_finite <- as.vector(!is.infinite(upper_CVXR))
    upperbound_constraints_CVXR <- diag(1*upper_finite)[upper_finite, ] %*% parvec <= matrix(upper_CVXR[upper_finite, ])

    # Quadratic Cone constraint
    isolate_v2 <- diag(c(
      rep(0, n_c),  # weights --- 1:n_c
      rep(0, p), #      v1 --- (n_c+1):(n_c+p)
      rep(1, p), #      v2 --- (n_c+p+1):(n_c+2p)
      rep(0, 4), #     t's --- (n_c+2p+1):(n_c+2p+4)
      rep(0, 6)  #      us --- (n_c+2p+5):(n_c+2p+10)
    ))
    find_t2 <- c(
      rep(0, n_c),  # weights --- 1:n_c
      rep(0, p), #      v1 --- (n_c+1):(n_c+p)
      rep(0, p), #      v2 --- (n_c+p+1):(n_c+2p)
      c(0, 0, 1, 0), #     t's --- (n_c+2p+1):(n_c+2p+4)
      rep(0, 6)  #      us --- (n_c+2p+5):(n_c+2p+10)
    )
    quadcone_constraints_CVXR <- p_norm(isolate_v2 %*% parvec, 2) <= t(find_t2) %*% parvec

    # Rotated Cone constraint
    isolate_t3 <- diag(c(
      rep(0, n_c),  # weights --- 1:n_c
      rep(0, p), #      v1 --- (n_c+1):(n_c+p)
      rep(0, p), #      v2 --- (n_c+p+1):(n_c+2p)
      c(0, 0, 0, 1), #     t's --- (n_c+2p+1):(n_c+2p+4)
      rep(0, 6)  #      us --- (n_c+2p+5):(n_c+2p+8)
    ))
    find_u11 <- c(
      rep(0, n_c),   #     weights --- 1:n_c
      rep(0, p),     #     v1 --- (n_c+1):(n_c+p)
      rep(0, p),     #     v2 --- (n_c+p+1):(n_c+2p)
      rep(0, 4),     #     t's --- (n_c+2p+1):(n_c+2p+4)
      c(1, 0, 0, 0, 0, 0)  #     us --- (n_c+2p+5):(n_c+2p+8)
    )
    rotquadcone1_constraints_CVXR <- sum_squares(isolate_t3 %*% parvec) <= t(find_u11) %*% parvec

    isolate_w_tildeu21 <- diag(c(
      rep(1, n_c),  # weights --- 1:n_c
      rep(0, p), #      v1 --- (n_c+1):(n_c+p)
      rep(0, p), #      v2 --- (n_c+p+1):(n_c+2p)
      rep(0, 4), #     t's --- (n_c+2p+1):(n_c+2p+4)
      c(0, 0, 0, 0, 1, 0)  #      us --- (n_c+2p+5):(n_c+2p+8)
    ))
    find_tildeu22 <- c(
      rep(0, n_c),   #     weights --- 1:n_c
      rep(0, p),     #     v1 --- (n_c+1):(n_c+p)
      rep(0, p),     #     v2 --- (n_c+p+1):(n_c+2p)
      rep(0, 4),     #     t's --- (n_c+2p+1):(n_c+2p+4)
      c(0, 0, 0, 0, 0, 1)  #     us --- (n_c+2p+5):(n_c+2p+8)
    )
    rotquadcone2_constraints_CVXR <- p_norm(isolate_w_tildeu21 %*% parvec, 2) <= t(find_tildeu22) %*% parvec # Note this is a rotated cone constraint that has been recoded as a quadratic cone constraint.

    # Get CVXR solution
    prob <- Problem(
      Minimize(t(c_CVXR) %*% parvec),
      constraints = c(
        equality_constraints_CVXR,
        lowerbound_constraints_CVXR,
        upperbound_constraints_CVXR,
        quadcone_constraints_CVXR,
        rotquadcone1_constraints_CVXR,
        rotquadcone2_constraints_CVXR
      )
    )
    result <- psolve(prob, verbose=verb, feastol=rtol, reltol=rtol, abstol=rtol, solver=solver, num_iter=maxit)

    # Pull out the parameter values
    parameters <- value(parvec)
    w <- rep(1, n_c + n_t)
    w[d == 0] <- parameters[1:n_c]

    return(w)
}
