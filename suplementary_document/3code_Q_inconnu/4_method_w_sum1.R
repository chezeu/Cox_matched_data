setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/3code_Q_inconnu")
source("2_EM_risk_function.R")

################# estimating equation of the weighted average equation 

equa_W_sum1 <- function(beta,Ts,event,XB, Q) {
  n=length(Ts)
  p = ncol(XB)
  X = as.matrix(XB, ncol = p) # covariates
  
  ##les average Z
  Z = Q%*%X
  ## the value q*exp(beta x) for tilfe_f
  eZbeta = exp(Z%*% beta)
  
  ##  the value q*x*exp(beta x) for tilde_g
  ZeZbeta = Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
  
  ## Estimating equation
  s =0
  for (i in 1:n) {
    
    risk = GetRiskSet(Ts[i], Ts)
    nrisk = length(risk)
    if(nrisk==1){
      t2R = eZbeta [risk] 
      t3R = ZeZbeta [risk,] 
      
    }else if(nrisk!=1){
      
      t2R = sum(eZbeta [risk] )
      t3R = colSums(ZeZbeta [risk,] ) 
    }
    
    s = s + event[i]* (Z[i,] - (t3R/t2R))
    
  }
  s
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
  if (iterations  == maxiter) {
    cat("WARNING! NOT CONVERGENT FOR W1_NEWTON!", "\n")
    converge = FALSE
  }
  return(list(coef = beta, converge = converge,iterations = iterations))
}

###################



