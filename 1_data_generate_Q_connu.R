
########################################################
# data generation of the bigger data base

Generate_data <- function(nB,nA,beta,Q,censor){
  ##covarites data
  p=3
  X1 = rnorm(nB,0,1)
  X2 = rbinom(nB,size = 1, prob = 0.7)
  X3 = rnorm(nB,0,2)
  X  = as.matrix(cbind(X1, X2, X3))
  
  L_true = sapply(1:nA, FUN = function(j){
    return( rmultinom(1,size=1, prob = Q[j,]))
  } )
  X_true = t(L_true)%*%X
  
  U  = runif(nA,0,1) # Time to event variable
  Tt = -log(U)/exp((X_true%*%beta))#lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = 0.7)
  c= censor
  
  Time      = pmin(Tt,c) # le vrai temps
  delta     = as.numeric(Tt<=c)
  surv_data = data.frame(Time,delta,X_true) # ou bien Z
  colnames(surv_data) = c('Time','delta',paste("X", 1:p, sep = ""))
  
  L_naive =sapply(1:nA, FUN= function(i){  which.max(Q[i,])
    })  
  X_naive = X[L_naive,]
  data_naive = data.frame(Time,delta, X_naive) 
  
  XB = as.matrix(X, ncol = p) 
  return(list(surv_data=surv_data, data_naive=data_naive,XB=XB) )
}


