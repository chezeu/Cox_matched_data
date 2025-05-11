

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/FS-method")
source("2_risk_function.R")
#################

equa_W_sum2 <- function(beta,Ts,event,XB, Q) {
  
  p = ncol(XB)
  n = length(Ts)
  X = as.matrix(XB, ncol = p) #Only covariates
  
  eXbeta = exp(X%*% beta)
  XeXbeta = X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  # the maximum values
  Z = matrix(0, nrow = n, ncol = p)
  som1 = vector()                  # denominator
  som2 = matrix(0, nrow = n, ncol = p) #numerator
  
  for (i in 1:n) {
    qi = Q[i,]
    i1 = which(qi == max(qi))
    l1 = max(qi)
    qi[i1] = 0
    i2 = which(qi == max(qi))
    if(length(i2)!= 1){i2 = i2[1]}
    qi[i1] = l1
    #normalize the probabilities
    p = qi[i1]+qi[i2]
    
    qix = ( (qi/p)*X)[c(i1,i2),]
    Z[i,] = colSums(qix)
    
    #### les sommes q*exp(beta x) pour tilfe_f   
    qixbeta = ((qi/p)*eXbeta)[c(i1,i2),]
    som1[i] = sum( qixbeta)
    
    ###  les sommes pour tilde_g   
    qiXeXbeta = ((qi/p)*XeXbeta)[c(i1,i2),]
    som2[i,] = colSums(qiXeXbeta)
  }  
  
  ## Estimating equation
  s = 0
  for (i in 1:n) {
 #at risk   
    risk = GetRiskSet(Ts[i], Ts)
    nrisk = length(risk)
    if(nrisk==1){
      t2R = som1[risk] 
      t3R = som2[risk,] 
      
    }else if(nrisk!=1){
      
      t2R = sum(som1[risk] )
      t3R = colSums(som2 [risk,] ) }
    
    s = s + event[i]*(Z[i,] - (t3R/t2R))
    
  }
  s 
  return(H_w_sum2 = s)  
}
# ########## 
# solve the equations 


coxph_w_sum2 <- function(Ts,event,XB, Q, maxiter = 20){
  p= ncol(XB)
  
  f <- function(x){
    equa_W_sum2  (beta=x,Ts,event,XB, Q)
  }
  
  fit_manual <- multiroot(f,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  iterations <- fit_manual$iter
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  if (iterations  == maxiter) {
    cat("WARNING! NOT CONVERGENT FOR W2_NEWTON!", "\n")
    converge = FALSE
  }
  return(list(coef = beta, converge = converge, iterations=iterations))
}
########