# sample linear dgp

makedf <- function(n = 200) {
  x1 <- rnorm(n)            # E 0, V 1
  x2 <- rnorm(n, sd = 0.25) # E 0, V 0.0625
  p <- exp(0.7*x1+x2)/(1+exp(0.7*x1+x2))
  d <- runif(n) < p
  e0 <- rnorm(n, 0, 2)
  e1 <- rnorm(n, 0, 1.5)
  y0 <- x1 + 2 * x2 + e0
  y1 <- x1 + 2 * x2 + 3 * x1 + x2 + e1
  y <- y0 * (1 - d) + y1 * d
  return(cbind(x1,x2,d,y,y0,y1))
}

set.seed(112)
tfb_sampledat1 <- makedf(1000)

# not run
#
# X <- df[,1:2]
# d <- df[,3]
# y <- df[,4]
#
# out <- tfb(X, d, y, fit = "lm", estimand = "att")

usethis::use_data(tfb_sampledat1)
