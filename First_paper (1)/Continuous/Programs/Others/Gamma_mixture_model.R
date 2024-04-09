set.seed(100)



logLikFun <- function(param, g, x) {
  p <- param[1]
  alpha_M <- param[2]
  beta_M <- param[3]
  alpha_U <- param[4]
  beta_U <- param[5]
  return( sum(g*log(p)+(1-g)*log(1-p))+sum(g*dgamma(x, shape = alpha_M, rate = beta_M, log = TRUE))
  + sum((1-g)*dgamma(x, shape = alpha_U, rate = beta_U, log = TRUE)) )
}

g <- c(rep(1,50), rep(0,4950))
x <- c(rgamma(50, shape = 1, rate = 0.5), rgamma(4950, shape = 9, rate = 0.2))
param <- c(1/100, 1,0.5,9,0.2)
logLikFun(param, g , x)

library(mixtools)
fit <- gammamixEM(x, verb = TRUE)
