# sample linear dgp

set.seed(211)

### Set sample size
n <- 1000

### Generate X's (Observed)
center_value <- 5
mean1 <- sample(c(0, 1), n, replace=T)
mean2 <- sample(c(0, 1), n, replace=T)

x1 <- center_value*mean1 + rnorm(n)
x2 <- center_value*mean2 + rnorm(n)

### Generate latent Z's (UNobserved)
z1 <- 1 / (sqrt((x1 - 0)^2 + (x2 - 0)^2) + 1)
z2 <- 1 / (sqrt((x1 - center_value)^2 + (x2 - 0)^2) + 1)
z3 <- 1 / (sqrt((x1 - 0)^2 + (x2 - center_value)^2) + 1)
z4 <- 1 / (sqrt((x1 - center_value)^2 + (x2 - center_value)^2) + 1)

### Create Propensity Score (UNobserved, based on unobserved Z's)
# Get which cluster they're in
center_index <- rep(0, n)
center_index[mean1==0 & mean2==0] <- 1
center_index[mean1==1 & mean2==0] <- 2
center_index[mean1==0 & mean2==1] <- 3
center_index[mean1==1 & mean2==1] <- 4

# Get true propensity score
temp <- 0.47
mx <-
  4*(center_index==1)*(z1-temp) +
  20*(center_index==2)*(z2-temp) +
  20*(center_index==3)*(z3-temp) +
  20*(center_index==4)*(z4-temp)
prob <- 1/(1 + exp(-mx))

### Generate D (Observed, based on unobserved propensity score)
d <- 1*(runif(n, min=0, max=1)<prob)

### Generate Y (Unobserved, based on unobserved Z's)
error <- rnorm(n)
f0 <- 10*z1 + 1*z2 + 1*z3 + 1*z4
error_sigma <- sqrt(2/3)*1.5
y <-  f0 + error_sigma*error
y0 <- y1 <- y

### Create data
tfb_sampledf3 <- cbind(x1,x2,d,y,y0,y1)

# not run
#
# X <- df[,1:2]
# d <- df[,3]
# y <- df[,4]
#
# out <- tfb(X, d, y, fit = "lm", estimand = "att")

usethis::use_data(tfb_sampledf3)
