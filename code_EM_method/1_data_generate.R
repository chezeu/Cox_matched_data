
########################################################

Generate_data <- function(m,n,beta){
  # set.seed(42)
  ##covarites data
  X1 = rnorm(m,0,1)
  X2 = rbinom(m,size = 1, prob = 0.8)
  X  = as.matrix(cbind(X1, X2))
  
  U  = runif(m,0,1) # Time to event variable
  Tt = -log(U)/exp((X%*%beta))#lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = 0.75)
  c  = 1.777907
  
  Time      = pmin(Tt,c) # le vrai temps
  delta     = as.numeric(Tt<=c)
  surv_data = data.frame(Time,delta,X)
  colnames(surv_data) = c('Time','delta',paste("X", 1:2, sep = ""))
  
  return(surv_data=surv_data )
}
