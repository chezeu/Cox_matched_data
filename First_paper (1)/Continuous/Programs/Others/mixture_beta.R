
library(maxLik)
loglik <- function(param, x, g){
  p = param[1]
  alpha_M = param[2]
  beta_M = param[3]
  alpha_U = param[4]
  beta_U = param[5]
  
  llValue = g*dbeta(x, shape1 = alpha_M, shape2 = beta_M, log = TRUE) 
          + (1-g)*dbeta(x, shape1 = alpha_U, shape2 = beta_U, log = TRUE)
          + g*ln(p)+(1-g)*ln(1-p)
  return(sum(llValue))
}

set.seed(100)
x <- c(rbeta(200, shape1 = 1,shape2 = 5), rbeta(200, shape1 = 0.5, shape2 = 2))
g = rep(0.5,400)
maxLik(logLik = loglik )
