# sample linear dgp

set.seed(1111)

### Set sample size### Set sample size
n <- 1000

### Generate Zs (Observed)
z1 <- rnorm(n)
z2 <- rnorm(n)
z3 <- rnorm(n)
z4 <- rnorm(n)

### Generate Y (Observed)
error <- rnorm(n)
f0 <- 8*z1 + 4*z2 + 2*z3 + 1*z4
error_sigma <- 9.21
y <-  f0 + (error_sigma * error)

### Propensity Score (Unobserved)
mx <- 0.20*z1 + 0.20*z2 + 0.20*z3 + 0.20*z4
prob <- 1/(1 + exp(-mx))

### Treatment Assignment (Observed)
d <- 1*(runif(n, min=0, max=1)<prob)

### Generate A variables (Observed)
p_A <- 5
A <- d + matrix(rnorm(n*p_A), nrow=n)
colnames(A) <- paste("a", 1:p_A, sep="")

### Generate U variables (Observed)
p_U <- 10
U <- matrix(rnorm(n*p_U), nrow=n)
colnames(U) <- paste("u", 1:p_U, sep="")

### Combine into dataset
y0 <- y1 <- y
tfb_sampledf4 <- cbind(y, y0, y1, d, z1, z2, z3, z4, A, U)

# not run
#
# X <- df[,-(1:4)]
# d <- df[,4]
# y <- df[,1]
#
# out <- tfb(X, d, y, fit = "lm", estimand = "att")

usethis::use_data(tfb_sampledf4)
