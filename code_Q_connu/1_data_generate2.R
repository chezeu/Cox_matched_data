
########################################################
# data generation of the bigger data base




Generate_data <- function(nB,nA,beta,Q){
  ##covarites data
  p=2
  X1 = rnorm(nB,0,1)
  X2 = rbinom(nB,size = 1, prob = 0.7)
  X  = as.matrix(cbind(X1, X2))
  
  L = sapply(1:nA, FUN = function(j){
    return( rmultinom(1,size=1, prob = Q[j,]))
  } )
  Z = t(L)%*%X 
  
  U  = runif(nA,0,1) # Time to event variable
  Tt = -log(U)/exp((Z%*%beta))#lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = 0.7)
  c  = 1.777907
  
  Time      = pmin(Tt,c) # le vrai temps
  delta     = as.numeric(Tt<=c)
  surv_data = data.frame(Time,delta,Z) # true data
  colnames(surv_data) = c('Time','delta',paste("X", 1:2, sep = ""))
  
  data_naive = data.frame(Time,delta,Z) 
  XB = as.matrix(X, ncol = p) 
  return(list(surv_data=surv_data, data_naive=data_naive,XB=XB) )
}


