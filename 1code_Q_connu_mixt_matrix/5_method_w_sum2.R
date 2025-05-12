

#setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/1code_Q_connu_mixt_matrix")
source("2_risk_function.R")
#################

equa_W_sum2 <- function(beta,Ts,event,XB, Q) {
  
  p = ncol(XB)
  n = length(Ts)
  X = as.matrix(XB, ncol = p) #Only covariates

  # the maximum values
  Z2 = matrix(0, nrow = n, ncol = p)
  #som1 = vector()                  # denominator
  #som2 = matrix(0, nrow = n, ncol = p) #numerator
  
  for (i in 1:n) {
    qi = Q[i,]
    i1 = which(qi == max(qi))
    if(length(i1)!= 1){i1 = i1[1]}
    l1 = qi[i1]
    qi[i1] = 0
    i2 = which(qi == max(qi))
    if(length(i2)!= 1){i2 = i2[1]}
    qi[i1] = l1
    # normalize the probabilities
    sum_q= qi[i1]+qi[i2]
    
    qix = ( (qi/sum_q)*X )[c(i1,i2),]
    Z2[i,] = colSums(qix)
  }
    #### les sommes q*exp(beta x) 
    
    eZ2beta = exp(Z2%*% beta)
    Z2eZ2beta = Z2*matrix(rep(exp(Z2%*% beta),p), ncol = p)
    
   # qixbeta = ( (qi/p)*eXbeta)[c(i1,i2),]
  #  som1[i] = sum( qixbeta)
    ###  les sommes pour tilde_g   
   # qiXeXbeta = ( (qi/p)*XeXbeta)[c(i1,i2),]
    #som2[i,] = colSums(qiXeXbeta)
    
  
  ## Estimating equation
  s = 0
  for (i in 1:n) {
 #at risk   
    risk = GetRiskSet(Ts[i], Ts)
    nrisk = length(risk)
    if(nrisk==1){
      t2R = eZbeta [risk] 
      t3R = ZeZbeta[risk,] 
      
    }else if(nrisk!=1){
      
      t2R = sum(eZ2beta [risk] )
      t3R = colSums(Z2eZ2beta [risk,] ) }
    
    s = s + event[i]*(Z2[i,] - (t3R/t2R))
    
  }
  s 
  return(s = s)  
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