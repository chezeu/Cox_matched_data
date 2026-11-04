
source("2_risk_function.R")

################# estimating equation of the Q weighted average equation 

equa_Q_sum <- function(beta,Ts,event,XB, Q) {
  n=length(Ts)
  p = ncol(XB)
  X = as.matrix(XB, ncol = p) # covariates
  
  eXbeta = exp(X%*% beta)
  XeXbeta = X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  ##les sommes pour Z
  Z = Q%*%X
  ## the value q*exp(beta x) for tilfe_f
  som1 =  Q%*%eXbeta
  
  ##  the value q*x*exp(beta x) for tilde_g
  som2 = Q%*%XeXbeta
  
  ## Estimating equation
  s =0
  for (i in 1:n) {
    
    risk = GetRiskSet(Ts[i], Ts)
    nrisk = length(risk)
    if(nrisk==1){
      t2R = som1[risk] 
      t3R = som2[risk,] 
      
    }else if(nrisk!=1){
      
      t2R = sum(som1[risk] )
      t3R = colSums(som2[risk,] ) 
    }
    
    s = s + event[i]* (Z[i,] - (t3R/t2R))
    
  }
  s
  return(s=s)
  
}
#####################
# ###########solve the equations 
coxph_Q_sum <- function(Ts,event,XB, Q, maxiter = 20){
  
  f <- function(x){
    equa_Q_sum  (beta=x,Ts,event,XB, Q)
  }
  
  fit_manual = multiroot(f,start = rep(0,p), maxiter = maxiter)
  beta = fit_manual$root
  iterations = fit_manual$iter
  converge = as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  if (iterations  == maxiter) {
    cat("WARNING! NOT CONVERGENT FOR W1_NEWTON!", "\n")
    converge = FALSE
  }
  return(list(coef = beta, converge = converge,iterations = iterations))
}

###################
