
# Monte-Carlo estimation of constant 
#censoring time that yields 40%, 30%, and 20%,  censoring

##covarites data
censor = NULL
for (i in 1:100000) {
  ##covarites data
  n=1000
  beta = c(0.5,-0.5, 1)
  
  X1 = rexp(n,0.5)
  X2 = rpois (n,2)
  # X3 = runif(n,0,1)
  X3 = rlnorm(n,0,1)
  X  = as.matrix(cbind(X1, X2, X3)) 
  
  U  = runif(n,0,1) # Time to event variable
  #Tt = sigma*(-log(U)/exp((X%*%beta)))^(1/alpha) # for weibull
  
  Tt = -log(U)/exp((X%*%beta)) #lambda=1
  
  # Constant right censored
  C =c(quantile(Tt, probs = 0.6), 
        quantile(Tt, probs = 0.7), # censor = 30%
        quantile(Tt, probs = 0.8))# censor=20%
  
  censor = rbind(censor, C)
  censor 
}
C = colMeans(censor)