
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tfb

<!-- badges: start -->
<!-- badges: end -->

## Overview

Targeted function balancing, or TFB, is a statistical method for causal
inference. `tfb` is the official package for the implementation of TFB
into the programming language R. TFB **Targets** a predictive
**Function** and **Balances** observations to determine its final
estimate. TFB thereby transforms the question of how best to weight
observational study data into how best to model the relationships in
that data.

## Installation

You can install the development version of `tfb` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bleth/tfb")
```

**PLEASE NOTE!** Currently, `tfb` solves its optimization problem with
`RMosek`, a package requiring additional installation steps. See
<https://docs.mosek.com/latest/rmosek/install-interface.html> for a
guide to installation. `RMosek` can be used for free with an educational
license.

## Use Case

### Causal Inference

`tfb` can help answer the question of what the *average effect* of a
binary variable is on an outcome. For example, we might want to know
what the effect of taking AP Calculus in high school had on college
students’ GPA. Or we might want to know whether the placement of a
product at eye level in a grocery store increases its sales; or any
other number of questions common to the practitioner of statistics. This
means that TFB answers questions of *inference*, not *prediction*.
Frequently inference is the subject of experimental studies: say,
medical researchers hire participants to determine the effectiveness of
their new medicine, through a hypothesis test or otherwise. This
methodology usually requires the *randomization* of treatment amongst
the participants, and may be referred to as *independence*. But what
happens when we aren’t able to fully randomize the treatment? Must we
always perform an experiment in order to make inference, or could we use
data from an observational study? *Causal inference* attempts to find
solutions to these problems, and TFB is a competitive method in the
field.

### Requirements

For a detailed treatment of the assumptions and conditions of `tfb`’s
best use case, please see the *Assumptions* and *Conditions* subsections
of the Background section of my paper, which can be found here:
<https://github.com/bleth/tfb/blob/main/papers/tfb_kbh.pdf>.

To properly use `tfb`, you’ll want a few things in your data:

1.  A matrix or data frame of covariates `X`. Covariates are also called
    variables or features. `X` is n by p in size and does not include
    the treatment.

2.  A dichotomous treatment vector `d`. This means `d` must take
    *exactly* two values. `d` is n by 1 in size.

3.  A numeric outcome vector `y`. `y` is n by 1 in size.

`tfb` will try to estimate the effect of `d` on `y`, given `X`.

You will also need to specify:

1.  The type of initial `fit` to use, representing what you think a
    reasonable representation of the relationship between `y` and `X`
    is. Currently, `tfb` supports ordinary least squares (`"ols"`),
    kernel-regularized least squares (`"krls"`), and elastic net
    regularization (`"elasticnet"`).

2.  The estimand, representing what subset of the data you want to
    estimate the average treatment effect for. Currently, `tfb` supports
    estimation of the average treatment effect on the treated (`"att"`),
    the average treatment effect on the controlled (`"atc"`), and the
    average treatment effect (`"ate"`).

## Simulation Study

### Data Exploration

TFB includes 4 sample datasets for the user. They have been created via
data generating process (DGP), the specifics of which are available at
<https://github.com/bleth/tfb/blob/main/data-raw/tfb_sampledat1.R>.
Looking at the first dataset:

``` r
df <- tfb::tfb_sampledat1
```

Our covariates here are `x1` and `x2`. Our treatment vector is `d`, and
our outcome vector is `y`. `y0` and `y1` are the *potential* outcomes,
which are not observable in practice.

``` r
X <- df[,1:2]
d <- df[,3]
y <- df[,4]
```

In this DGP, While `d` affects `y`, `x1` and `x2` also affect both `d`
and `y`. This means that the treatment effect of `d` is *confounded* by
the covariates. Fortunately, we are able to observe all of the
confounders in our data. We also know that `d` only takes two values,
and that the treatment values are *stable* (our observations are
independent and the treatment received is consistent between
observations). In a practical application of `tfb`, the user should
speculate on whether these assumptions hold, and exercise caution in
reporting estimates when they do not hold.

From our DGP, we also know that the relationship between the outcome and
the covariates is *linear*, meaning an initial `fit` of ordinary least
squares is appropriate. In practice, the efficacy of `tfb` depends on
how accurately the user can model this relationship. Empirically, we
might choose a predictive model through cross-validation:

``` r

library(tfb)
library(KRLS)
#> ## KRLS Package for Kernel-based Regularized Least Squares.
#> ## See Hainmueller and Hazlett (2014) for details.
set.seed(2062)


# splitting training and testing set
indices <- sample(1:nrow(df),200)
train <- df[indices, ]
test <- df[-indices, ]

X <- train[,1:2]
d <- train[,3]
y <- train[,4]


# linear cross-validated error
cve_ols <- function(X, y, folds) {
  
  data <- as.data.frame(cbind(y,X))
  colnames(data) <- c("y", paste0("x", 1:ncol(X)))
  formula <- as.formula(paste0(
    "y ~ ", paste(colnames(data)[-1], collapse = " + ")
  ))
  
  n <- nrow(data)
  sampling <- sample(rep(seq_len(folds),(n %/% folds)+1)[1:n])
  
  rmses <- c()
  for (i in seq_len(folds)) {
    model <- lm(formula, data = data[sampling != i, ])
    rmses <- c(rmses,sqrt(mean(
      (y[sampling == i] - predict(model, data[sampling == i, ]))^2
    )))
  }
  
  return(mean(rmses))
}


# nonlinear cross-validated error
cve_krls <- function(X, y, folds) {
  
  n <- nrow(X)
  sampling <- sample(rep(seq_len(folds),(n %/% folds)+1)[1:n])
  
  rmses <- c()
  for (i in seq_len(folds)) {
    model <- krls(
      X[sampling != i, ], y[sampling != i], print.level = 0
    )
    rmses <- c(rmses,sqrt(mean(
      (y[sampling == i] - predict(model, X[sampling == i, ])$fit)^2
    )))
  }
  
  return(mean(rmses))
}


# ols cves
cve_ols <- c(
  cve_ols(X, y, 5),
  cve_ols(X[d == 0, ], y[d == 0], 5),
  cve_ols(X[d == 1, ], y[d == 1], 5)
)


# krls cves
cve_control <- c(
  cve_krls(X, y, 5),
  cve_krls(X[d == 0, ], y[d == 0], 5),
  cve_krls(X[d == 1, ], y[d == 1], 5)
)

matrix(
  c(cve_ols,cve_control), ncol=2,
  dimnames = list(
    c("Full","Control","Treatment"),
    c("OLS CVE", "KRLS CVE")
  )
)
#>            OLS CVE KRLS CVE
#> Full      2.292737 2.427177
#> Control   2.183838 2.145914
#> Treatment 1.341460 1.932812
```

It seems like a simple linear model is a better choice here in terms of
predictive power compared to the highly flexible kernel regularized
least squares (KRLS). However, in practice we cannot be sure of the
relationship between the outcome and the covariates, and it may be best
to rely on field knowledge. Let’s further evaluate our linear fit:

``` r

# retrieving correlation and p-values
adj_rsq <- c(
  summary(lm(y ~ X))[9],
  summary(lm(y[d == 0] ~ X[d == 0, ]))[9],
  summary(lm(y[d == 1] ~ X[d == 1, ]))[9]
)

matrix(
  adj_rsq, ncol = 1, 
  dimnames = list(
    c("Full","Control","Treatment"),
    c("Adjusted R2")
  )
)
#>           Adjusted R2
#> Full      0.6565209  
#> Control   0.2751468  
#> Treatment 0.9170279
```

In truth, TFB targeting OLS (TFB-OLS) will be unbiased for any estimand
here, because the DGP is truly linear. But in practice, we may be
hesitant to estimate the ATE and the ATT, which train a model on the
control set. The predictive power in the treatment set seems to be good,
suggesting the choice of the ATC as an estimand. All that’s left to
check is the distributions of the covariates.

``` r

# treatment group ranges
x1_min <- min(X[d == 1, 1])
x2_min <- min(X[d == 1, 2])
x1_max <- max(X[d == 1, 1])
x2_max <- max(X[d == 1, 2])

# proportion of control observations uncovered in x1
mean((X[d == 0, 1] < x1_min) | (X[d == 0, 1] > x1_max))
#> [1] 0.08

# proportion of control observations uncovered in x2
mean((X[d == 0, 2] < x2_min) | (X[d == 0, 2] > x2_max))
#> [1] 0.03
```

For the ATC, we weight the treatment group to fit the control group, so
we ideally want the distribution of covariates in the control group to
be a *subset* of the distribution of covariates in the treatment group.
Disjoint or particularly disparate treatment and control groups can
introduce bias, but fortunately that isn’t a worry in this instance.
With our preliminary analysis concluded, we will use the TFB-OLS
estimator for the ATC:

``` r

X <- test[,1:2]
d <- test[,3]
y <- test[,4]

out <- tfb(X, d, y, fit = "ols", estimand = "atc")

matrix(out$final$estimate, dimnames = list(NULL, "Final Estimate"))
#>      Final Estimate
#> [1,]      -0.977017
```

In this sample, our final estimate is pleasantly close to the true ATC
of -1, because we have correctly specified our `fit`. In practice, we
don’t know the true treatment effect nor its empirical analog, but if
TFB is applied with care, the method provides an estimator that is less
naive than a difference in means (DIM).

### Analysis

What else can `tfb` tell us? Its output is a nested list, where we can
find summary statistics as well as raw information about the
optimization.

``` r

# the final estimate
# out$final$estimate

# the variance of the final estimate
# out$final$variance

# a confidence interval for the final estimate
# out$final$confidence_interval

# a p-value for the one-sided hypothesis test
# that the final estimate is greater or less than zero
# out$final$p_value

# the naive difference in means, for comparison
# out$final$dim

out$final
#> $estimate
#> [1] -0.977017
#> 
#> $variance
#> [1] 0.04743468
#> 
#> $confidence_interval
#> $confidence_interval$lower
#> [1] -1.403887
#> 
#> $confidence_interval$upper
#> [1] -0.5501466
#> 
#> 
#> $p_value
#> [1] 3.62954e-06
#> 
#> $dim
#> [1] 1.944528
```

From the summary statistics, we have our estimate and a measure of its
variance. The confidence interval and p-value both suggest that the true
ATC is less than zero. We can also see that the naive DIM is quite far
off from our estimate, and quite far off from the true ATC given our
knowledge of the DGP. The rest of the output is split up by sample. This
includes summary statistics for within the split, which will be
identical to those above in the case that only one sample is taken.
However, the lower level summary also includes TFB’s weights for that
sample:

``` r

w <- out$sample_1$final$weights

cbind(w, X, d, y)[1:10, ]
#>                  w         x1          x2 d          y
#>  [1,] 1.000000e+00 -0.3142317  0.23960407 0  0.2297696
#>  [2,] 7.052556e-05  2.4033751 -0.18627112 1  6.5574659
#>  [3,] 2.143326e+00 -0.7182640 -0.06589086 1 -1.9735341
#>  [4,] 1.000000e+00 -1.7606110 -0.08700839 0 -2.4294093
#>  [5,] 1.000000e+00 -1.1252812  0.37375539 0 -0.8127416
#>  [6,] 7.631146e-04  1.3102240  0.06165787 1  5.1013499
#>  [7,] 1.132537e+00  0.1521774  0.20096113 1  2.4857978
#>  [8,] 1.172403e+00  0.6533708 -0.23710188 1  2.8999608
#>  [9,] 1.000000e+00 -0.9473355  0.43233386 0  1.9303582
#> [10,] 1.000000e+00 -0.1760594 -0.24134621 0  0.4716542
```

TFB’s estimator is a weighted difference in means (WDIM), so naturally
the weights can be used to manually calculate its estimates. The weights
may also be investigated for trends in covariate importance. Estimating
the ATC with TFB involves weighting the treatment set, and we can
similarly see that all the control units have the default weight of 1.
The rest of the output within a sample is split up by fold.

``` r

# fold indices
# out$sample_1$indices

# model coefficients
out$sample_1$betas
#> $beta_c1
#> [1] 3.923368 3.053705
#> 
#> $beta_c2
#> [1] 3.848764 2.570067

# model covariance matrices
# out$sample_1$covariances

# model residuals
# out$sample_1$errors

# model standard error
# out$sample_1$standard_errors

# fold dims
# out$sample_1$dims

# fold wdims
out$sample_1$wdims
#> NULL

# fold wdim variances
# out$sample_1$variances
```

This can be used for detail work involving the predictive models fit as
a part of TFB. For example, we can see the coefficients of the linear
model fit on the treatment group in each fold, and the WDIMs within each
fold, which are averaged to determine the final estimate within the
sample.
