
#source("2_risk_function.R")
################### cumulative function

Funct_lambda2<-function(lambda0,Ts){
  #lambda0<- rep(0.1,length(event))
  n = length(Ts)
  #lambda2 = vector()
  lambda2=numeric(n)
  
  for (i in 1:n) {
    lambda2[i] =  sum(lambda0[which(Ts<=Ts[i])] )  
  }
  return(lambda2=lambda2)
}
########
### proba aposteriories (pi)

Functio_prob<-function(beta0,lambda2,Ts,event,XB, Q){
  
  XB=as.matrix(XB)
  p=ncol(XB)
  X = XB # covariates
  
  eXbeta0 = as.vector(exp(X %*% beta0))
  
  n = length(Ts)
  m = nrow(XB)
  
  #pi(ij)
  
  prob = matrix(0,n,m)
  for (i in 1:n) {
    
    numerateur = vector()
    for (j in 1:m) {
      numerateur[j] = Q[i,j] * ( eXbeta0[j])^event[i] * exp(-lambda2[i]* eXbeta0[j] )
    }
    
    denominateur = sum(numerateur)
    if (denominateur <= 0 || !is.finite(denominateur)) denominateur <- .Machine$double.eps
    prob[i,] = numerateur/denominateur
  }
 # prob
  return( prob = prob)
}
#############
################ estimations ###################################
#lambda0 ( baseline hazard function)
Function_lambda0<-function(prob,beta0,Ts,event,XB){
  XB=as.matrix(XB)
  p = ncol(XB)
  n = length(Ts)
  
  X = XB # covariates
  eXbeta0 = as.vector(exp(X%*% beta0))
  
  som1 = as.vector(prob%*%eXbeta0) # sum with j
  
  
  #lambda0_1 = vector()
  lambda0_1 <- numeric(n)
  
  for (i in 1:n) {
    
    if(event[i] == 0){lambda0_1[i] = 0
    }else{
      
      risq = GetRiskSet (Ts[i], Ts)
      nrisk = length(risq)
      
      if(nrisk==1){
        t2R = som1[risq]  # denominator
      }else if(nrisk!=1){
        
        t2R = sum( som1[risq] )
      }
        if (!is.finite(t2R) || t2R <= 0) t2R <- .Machine$double.eps # avoid the form log(0)
      lambda0_1[i] = (1/t2R)
    }
  }
 # lambda0_1
  return( lambda0_1 = lambda0_1)
}

######################
#  estimating equation ########################

equa_estimate <- function(beta,prob,Ts,event,XB) {
  
  XB=as.matrix(XB)
  p = ncol(XB)
  n = length(Ts)
  
  X = XB # covariates
  eXbeta = as.vector(exp(X %*% beta))
  #XeXbeta = X*matrix(rep(exp(X%*% beta),p), ncol = p)
  XeXbeta = X*matrix(rep(eXbeta,p), ncol = p)
  
  #######les sommes pour Z
  
  Z = prob%*%X 
  
  #### values for q*exp(beta x) pour tilfe_f
  
  som1 = as.vector(prob%*%eXbeta)
  
  ###  values for q*x*exp(beta x)
  
  som2 = prob%*%XeXbeta
  
  
  ## Estimating equation
  s <- rep(0, p)
  
  for (i in 1:n) {
    
    risk = GetRiskSet(Ts[i], Ts)
    nrisk = length(risk)
    
    if(nrisk==1){
      t2R = som1[risk] 
      t3R = som2[risk,] 
    }else if(nrisk!=1){
      
      t2R = sum( som1[risk] )
      t3R = colSums(som2[risk,] ) }
    
    s= s+ event[i]*(Z[i,] - t3R/t2R)
    
  }
  #s 
  
  return(as.numeric(s))
  
}

#####################
# ###########solve the equation 

coxph_estimate<- function(prob,Ts,event,XB,beta_ini,maxiter = 50){
  
  f <- function(x){
   val<- equa_estimate  (beta= x,prob,Ts,event,XB)
   if (!is.numeric(val) || length(val) != length(x)) {
     stop("equa_estimate must return a numeric vector of same length as x")
   }
   return(val)
  }
   fit_manual <- nleqslv( beta_ini,f, method = c("Broyden", "Newton"))
   beta0 <- as.numeric(fit_manual$x)
  iterations <- fit_manual$iter
   converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1) )
 # converged <- identical(fit$termcd, 1L) || identical(fit$termcd, 2L)
 # fit_manual = multiroot(f,start = beta_ini, maxiter = maxiter)
  #beta_hat = fit_manual$root
  #iterations =  fit_manual$iter
 # converge = as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  #converge = is.finite(fit_manual$estim.precis) && !is.na(fit_manual$estim.precis) && (fit_manual$estim.precis < 1e-6) && (iterations < maxiter)
  if (iterations  == maxiter) {
    cat("WARNING! NOT CONVERGENT FOR NEWTON!", "\n")
    converge = FALSE
  }
  return(list( beta0 = beta0, converge = converge, iterations = iterations) )
}
################

### iterrations ###########################################

#valeurs initials

Func_itteration<-function(beta0,lambda0,Ts,event,XB, Q,tol= 1e-6, maxits = 500){
  
  XB=as.matrix(XB)
  p = ncol(XB)
  n = length(Ts)
  
  beta_ini = rep(0, p)
  it = 0
  converge = FALSE
  
  while ((!converge) && (it < maxits)){ 
    
    lambda0_old = lambda0
    beta0_old = beta0
    
    #expectation  
    
    lambda2 = Funct_lambda2(lambda0_old,Ts)
    
    prob = Functio_prob(beta0_old,lambda2,Ts,event,XB, Q)
    
    
    #maximization
    lambda0 = Function_lambda0 (prob,beta0_old,Ts,event,XB)
    estime = coxph_estimate (prob,Ts,event,XB,beta_ini = beta0_old , maxiter = 50)
    beta0 = estime$beta0
    #beta_ini = beta0.old
    
    d_beta   <- sqrt(sum((beta0 - beta0_old)^2))
    d_lambda <- sqrt(sum((lambda0 - lambda0_old)^2))
    converge <- (d_beta + d_lambda) < tol
    #converge = sqrt (sum ((beta0.old -beta0))^2 + sum ( (lambda0.old - lambda0)^2 ) ) < tol  
    
    if (it == maxits && !converge) {
      cat("WARNING! NOT CONVERGENT WITH EM!", "\n")
      converge = FALSE
    }
    if (any(is.na(beta0))) {
      cat("WARNING! beta0 NOT AVAILABLE!", "\n")
      converge <- FALSE
    }
    
    it = it + 1
    
  }
  return(list(beta0=beta0, lambda0=lambda0, prob=prob, converge= converge, it=it))
}


