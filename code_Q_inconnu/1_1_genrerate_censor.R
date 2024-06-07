
# Monte-Carlo estimation of constant censoring time that yields 30% censoring

##covarites data
censor = NULL
for (i in 1:100000) {
  ##covarites data
  n=1000
  beta = c(0.5,-0.5, 1)
  
  X1 = rnorm(n,0,1)
  X2 = rbinom(n,size = 1, prob = 0.7)
  X3 = rnorm(n,0,2)
  X  = as.matrix(cbind(X1, X2, X3)) 
  
  U  = runif(n,0,1) # Time to event variable
  Tt = -log(U)/exp((X%*%beta)) #lambda=1
  # Constant right censored
  c = quantile(Tt, probs = 0.7)
  
  censor = c(censor, c)
  censor 
}
c = mean(censor)

