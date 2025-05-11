setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/1code_Q_connu_mixt_matrix")
source("2_risk_function.R")

################# estimating equation of the weighted average equation 

equa_W_sum1 <- function(beta,Ts,event,XB, Q) {
  n=length(Ts)
  p = ncol(XB)
  X = as.matrix(XB, ncol = p) # covariates
  
  ##les average Z
  Z1 = Q%*%X
  ## the value q*exp(beta x) for tilfe_f
  eZ1beta = exp(Z1%*% beta)
  
  ##  the value q*x*exp(beta x) for tilde_g
  Z1eZ1beta = Z1*matrix(rep(exp(Z1%*% beta),p), ncol = p)
  
  ## Estimating equation
  s =0
  for (i in 1:n) {
   
    risk = GetRiskSet(Ts[i], Ts)
    nrisk = length(risk)
    if(nrisk==1){
      t2R = eZ1beta [risk] 
      t3R = Z1eZ1beta [risk,] 
      
    }else if(nrisk!=1){
      
      t2R = sum(eZ1beta [risk] )
      t3R = colSums(Z1eZ1beta [risk,] ) 
      }
    
    s = s + event[i]* (Z1[i,] - (t3R/t2R))
   
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



