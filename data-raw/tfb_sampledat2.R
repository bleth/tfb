# sample nonlinear dgp

makedf <- function(n) {
  x1 <- rchisq(n, df = 3) - 3 # E 0, V 6
  x2 <- rbeta(n, 0.5, 0.5)    # E 0.5, V 0.125
  prob <- 1/(1+exp(x1*x2+0.5))
  d <- runif(n) < prob
  e0 <- rnorm(n, 0, 0.25)
  e1 <- rnorm(n, 0, 0.5)
  y0 <- (x1 - 3) * (x2 - 2) + e0
  y1 <- (x1 - 3) * (x1 + 3) * (x2 - 2) * (x2 + 3) + e1
  y <- y0 * (1 - d) + y1 * d
  return(cbind(x1,x2,d,y,y0,y1))
}

set.seed(121)
tfb_sampledat2 <- makedf(1000)

# not run
#
# X <- df[,1:2]
# d <- df[,3]
# y <- df[,4]
#
# out <- tfb(X, d, y, fit = "krls", estimand = "att")

usethis::use_data(tfb_sampledat2)
