#' Targeted Function Balancing Optimization -- ATT
#'
#' This function solves the linear optimization problem presented by TFB, given information from an initial regression, with the ATT as the estimand.
#' @param X covariate matrix
#' @param beta coefficient vector
#' @param sqrtV square root of variance matrix
#' @param sigma2 estimated variance
#' @param i fold indices
#' @param d treatment status vector
#' @param chi_q probability threshold
#' @param quiet whether to suppress console output
#'
#' @returns A numeric vector of weights.
#'
#' @import Rmosek
#' @import SparseM
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
#' tfb:::tfb_balance_att(X,beta,sqrtV,sigma2,i,d,0.95,TRUE)


### Optimization Function
tfb_balance_att <- function(

    X,
    beta,
    sqrtV,
    sigma2,
    i,
    d,
    chi_q,
    quiet

  ){

    rtol <- 1e-16
    if (quiet) {
      verb <- 0
    } else {
      verb <- 10
    }
    d <- d[i == 1]

                     # Define values we need
    n_c <- sum(d==0) # control sample size
    n_t <- sum(d==1) # treated sample size
    X_c <- X[d==0, ] # control X
    X_t <- X[d==1, ] # treated X
    p <- ncol(X)     # number of covariates

    # Want minimization
    P <- list(sense = "min")

    # Define linear function
    P$c <- c(
      rep(0, n_c),           # w_c
      rep(0, p),             # v_1c
      rep(0, p),             # v_2c
      rep(0, 2),             # t_1c
      0,                     # t_2c
      0,                     # t_3
      c(1, 0),               # u_1
      c(sigma2 / (n_c^2), 0) # u_2c
    )

    # Define Linear constraint matrices
    A.sum_to_n_c <- c(
      rep(1, n_c),              # for weights
      rep(0, length(P$c) - n_c) # everything else
    )
    A.sum_to_n_c <- t(as.matrix(A.sum_to_n_c)) # Just turning it into a matrix of the right dimension
    bc.sum_to_n_c <- matrix(n_c, nrow=2)         # upper and lower-bound on the constraints
      ### Here we'll want       LB <= A.sum_to_n_c * theta <= UB

    A.balance1 <- cbind(
      (1 / n_c) * t(X_c),        # for weights
      diag(rep(1, p)),           # for v1
      diag(rep(0, p)),           # for v2
      matrix(0, nrow=p, ncol=4), # for t's
      matrix(0, nrow=p, ncol=4)  # for u's
    )
    bc.balance1 <- matrix(colMeans(X_t), nrow=2, ncol=p, byrow=T)

    A.balance2 <- cbind(
      (1 / n_c) * sqrtV %*% t(X_c), # for weights
      diag(rep(0, p)),              # for v1
      diag(rep(1, p)),              # for v2
      matrix(0, nrow=p, ncol=4),    # for t's
      matrix(0, nrow=p, ncol=4)     # for u's
    )
    bc.balance2 <- matrix(colMeans(X_t %*% sqrtV), nrow=2, ncol=p, byrow=T)

    A.mag_imbal <- c(
      rep(0, n_c), # for w
      beta,       # for v1
      rep(0, p),  # for v2
      c(-1, 1),   # for t_1's
      rep(0, 2),  # for other t's
      rep(0, 4)   # for u's
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
      rep(0, 4)                  # u's
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
      matrix(c(0, 1), nrow=2, ncol=1)  # for u_2^(2)
    )
    bc.u2s <- matrix(0.5, nrow=2, ncol=2)

    A <- rbind(
      A.sum_to_n_c,
      A.balance1,
      A.balance2,
      A.mag_imbal,
      A.bias,
      A.u2s
    )
    A <- SparseM::as.matrix.csr(A)
    P$A <- as(A, "CsparseMatrix")

    bc <- cbind(
      bc.sum_to_n_c,
      bc.balance1,
      bc.balance2,
      bc.mag_imbal,
      bc.bias,
      bc.u2s
    )
    P$bc <- bc

    # Define upper and lower bounds for parameters
    lower <- c(
      rep(0, n_c),   # weights --- 1:n_c
      rep(-Inf, p), #      v1 --- (n_c+1):(n_c+p)
      rep(-Inf, p), #      v2 --- (n_c+p+1):(n_c+2p)
      rep(0, 4),    #     t's --- (n_c+2p+1):(n_c+2p+4)
      rep(0, 4)     #      us --- (n_c+2p+5):(n_c+2p+8)
    )
    lower <- t(as.matrix(lower))

    upper <- c(
      rep(n_c, n_c),  # weights --- 1:n_c
      rep(Inf, p), #      v1 --- (n_c+1):(n_c+p)
      rep(Inf, p), #      v2 --- (n_c+p+1):(n_c+2p)
      rep(Inf, 4), #     t's --- (n_c+2p+1):(n_c+2p+4)
      rep(Inf, 4)  #      us --- (n_c+2p+5):(n_c+2p+8)
    )
    upper <- t(as.matrix(upper))

    P$bx <- rbind(
      lower, upper
    )

    # Define cones
    P$cones <- matrix(list(), 2, 3)
    rownames(P$cones) <- c("type", "sub")

    P$cones[1,1] <- "QUAD"
    P$cones[2,1] <- list(c(
      n_c+(2*p)+3, # t2
      (n_c+p+1):(n_c+(2*p)) # v2
    ))

    P$cones[1,2:3] <- "RQUAD"
    P$cones[2,2] <- list(c(
      n_c+(2*p)+5, n_c+(2*p)+6, # u_1^{(1)}, u_1^{(2)}
      n_c+(2*p)+4 # t3
    ))
    P$cones[2,3] <- list(c(
      n_c+(2*p)+7, n_c+(2*p)+8, # u_2^{(1)}, u_2^{(2)}
      1:n_c # w
    ))

    # define tolerance level
    P$dparam <- list(
      MSK_DPAR_ANA_SOL_INFEAS_TOL = rtol,
      MSK_DPAR_BASIS_REL_TOL_S = rtol
    )

    # Get solution
    solution <- Rmosek::mosek(P, opts = list(verbose = verb))

    if (!quiet) {message(solution$sol$itr$solsta)}

    parameters <- solution$sol$itr$xx # define parameters for the solution
    w <- rep(1, n_c + n_t)
    w[d == 0] <- parameters[1:n_c]

    return(w)
}
