
#########  multinomiale Q
setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/suplementary document/3code_Q_inconnu")

source("2_EM_risk_function.R")

generate_naive_data <- function(naive_id,XB,surv_data){
  
  X_naive <- XB[naive_id,]
  data_naive = cbind(surv_data[,1:2],X_naive)
  data.frame(data_naive)                 
  return(data_naive=data_naive)
}

#################
# equation naive
equa_naive <- function(beta,Ts,event, Z) {
  
  p=ncol(Z)
  ezbeta = exp(Z%*% beta)
  zezbeta = Z*matrix(rep(exp(Z%*% beta),p), ncol = p) 
  
  n = length(Ts)
  
  s = 0
  for(i in 1:n){ 
  # at risk
   # Yi = as.numeric( (Ts == Ts[i] ) | (Ts > Ts[i]))# at risk
    risq = which(( Ts >= Ts[i]))
    if(length(risq)==1){
      
      num_naive = zezbeta[risq,]
      denum_naive = ezbeta[risq] 
      
    }else if (length(risq)!=1){
      
      num_naive = colSums( zezbeta[risq,])
      denum_naive = sum(ezbeta[risq])  }  
    
    s = s+ event[i]* ( Z[i,]- num_naive/denum_naive)
  }
  s
  return(H_naive = s) ##naive estimator 
}

##############

# solve the naive equation 

coxph_equa_naive <- function(Ts,event, Z, maxiter = 20){
  p=ncol(Z)
  
  f <- function(x){
    equa_naive (beta=x,Ts,event, Z)
  }
  
  fit_manual = multiroot(f,start = rep(0.2,p), maxiter = maxiter)
  beta = fit_manual$root
  iterations = fit_manual$iter
  converge = as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  if (iterations  == maxiter) {
    cat("WARNING! NOT CONVERGENT FOR NAIVE_NEWTON!", "\n")
    converge = FALSE
  }
  return(list(coef = beta, converge = converge, iterations = iterations))
}
#################
