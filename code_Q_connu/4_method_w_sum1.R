

##################

#1. Finding the risk set R(t) given some time t
GetRiskSet <- function(time_of_interest, time_vector, event_vector) {
  
  return(which(((time_vector == time_of_interest & event_vector == 1) | (time_vector > time_of_interest))))
  
}
######################


################# estimating equation

equa_W_sum1 <- function(beta,Ts,event,XB, Q) {
  
  p = ncol(XB)
  X = as.matrix(XB, ncol = p) # covariates
  
  eXbeta = exp(X%*% beta)
  XeXbeta = X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  ##les sommes pour Z
  Z = Q%*%X
  
  ## les sommes q*exp(beta x) pour tilfe_f
  som1 =  Q%*%eXbeta
  
  ##  les sommes pour tilde_g
  som2 = Q%*%XeXbeta
  
  dat1 = cbind(Ts,Z)[which(event==1),] #new  data
  
  ## Estimating equation
  s = matrix(0,nrow=nrow(dat1), ncol = p)
  for (i in 1:nrow(dat1)) {
    ts = dat1[i, 1]
    Z1R = dat1[i,2:ncol(dat1)]
    
    risk = GetRiskSet(ts, Ts, event)
    nrisk = length(risk)
    if(nrisk==1){
      t2R = som1[risk] 
      t3R = matrix(som2, ncol = p)[risk,] 
      
    }else if(nrisk!=1){
      
      t2R = sum(som1[risk] )
      t3R = colSums(matrix(som2, ncol = p)[risk,] ) }
    
    s[i,] = (Z1R - (t3R/t2R))
    
  }
  
  s = colSums(s)
  
  return(s=s)
  
}
#####################
# ###########solve the equations 
coxph_w_sum1 <- function(Ts,event,XB, Q, maxiter = 20){
  
  f <- function(x){
    equa_W_sum1  (beta=x,Ts,event,XB, Q)
  }
  
  fit_manual = multiroot(f,start = rep(0,p), maxiter = maxiter)
  beta = fit_manual$root
  iterations = fit_manual$iter
  converge = as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list(coef = beta, converge = converge,iterations = iterations))
}

###################



