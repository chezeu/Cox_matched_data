


setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")
source("2_risk_function.R")

################### cumulative function




#death_times= sort(unique(Ts*event))
#death_times=death_times[death_times!=0]

### proba aposteriories (pi)

Functio_prob<-function(beta0,lambda0,Ts,event,XB,Q){
  #beta0<- c(0.1,0.1) 
  #lambda0<- rep(0.1,length(event))
  #lambda0[which(lambda0==0)]=0
  
  X = as.matrix(XB, ncol = p) # covariates
  eXbeta0 = exp(X%*% beta0)
  n = length(Ts)
  m = nrow(XB)
  
  #pi(ij)
  
  prob = matrix(0,n,m)
  for (i in 1:n) {
    
    numerator = vector()
    for (j in 1:m) {
      numerator[j] = Q[i,j] * (eXbeta0[j])^event[i] * exp( - sum(lambda0[which(Ts<=Ts[i])]) * eXbeta0[j] )
    }
    
    denominator = sum(numerator)
    
    prob[i,] = numerator/denominator
  }
  prob
  return( prob = prob)
}
#############
################ estimations ###################################
#lambda0 ( baseline risk function)

Function_lambda0<-function(prob,beta0,Ts,event,XB){
  
  p = ncol(XB)
  n = length(Ts)
  
  X = as.matrix(XB, ncol = p) # covariates
  eXbeta0 = exp(X%*% beta0)
  
  som = prob%*%eXbeta0 # somme sur les j
  
  lambda0 = vector()
  
  for (i in 1:n) {
    if(event[i]==0){lambda0[i] = 0
      }else{ 
      lambda0[i] = 1/ (as.numeric( Ts[i] <= Ts)%*%som )
      }
  }
  return( lambda0 = lambda0)
}

######################
#  estimating equation ########################

equa_estimate <- function(beta,prob,Ts,event,XB) {
  
  p = ncol(XB)
  n = length(Ts)
  
  X = as.matrix(XB, ncol = p) # covariates
  eXbeta = exp(X%*% beta)
  XeXbeta = X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  #######les sommes pour Z
  
  Z = prob%*%X 
  
  #### values for q*exp(beta x) 
  
  denominator = prob%*%eXbeta
  
  ###  values for q*x*exp(beta x)
  
  numerator = prob%*%XeXbeta
  
  dat1 = cbind(Ts, Z)[which(event==1),]
  
  ## Estimating equation
  s = matrix(0,nrow=nrow(dat1), ncol = p)
  
  for (i in 1:nrow(dat1)) {
    ts = dat1[i, 1]
    Z1 = dat1[i,2:ncol(dat1)]
    
    risk = GetRiskSet(ts, Ts)
    nrisk = length(risk)
    
    if(nrisk==1){
      denom =  denominator[risk] 
      num = matrix(numerator, ncol = p)[risk,] 
    }else if(nrisk!=1){
      
      denom = sum( denominator[risk] )
      num = colSums(matrix(numerator, ncol = p)[risk,] ) 
      }
    
    s[i,] = (Z1 - num/denom)
    
  }
  
  s <- colSums(s) 
  
  return(H_EM = s)
  
}

#####################
# ###########solve the equation 

coxph_estimate<- function(prob,Ts,event,XB,beta_ini,maxiter = 20){
  
  f <- function(x){
    equa_estimate  (beta= x,prob,Ts,event,XB)
  }
 # fit_manual <- nleqslv( c(beta_ini[1],beta_ini[2]),f, method = c("Broyden", "Newton"))
 # beta0 <- fit_manual$x
  #iterations <- fit_manual$iter
 # converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1) )
  
  #xstart <- matrix(rnorm(20,0,1), ncol = 2)
  #Zero <-  searchZeros(xstart,f)
  
  fit_manual = multiroot(f,start = c(beta_ini[1],beta_ini[2]), maxiter = maxiter)
  beta0 = fit_manual$root
  iterations =  fit_manual$iter
  converge = as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
 
   if (converge==0 ) {
    cat("WARNING! NOT CONVERGENT FOR NEWTON!", "\n")
  }
  return(list( beta0 = beta0, converge = converge, iterations = iterations) )
}
################

### iterrations ###########################################

#valeurs initials

Func_itteration<-function(beta0,lambda0,Ts,event,XB,Q,tol= 1e-6, maxits = 500){
  
  p = ncol(XB)
  n = length(Ts)
  beta_ini = rep(0, p)
  it = 0
  converge = FALSE
  
  while ((!converge) & (it < maxits)){ 
    
    lambda0.old = lambda0
    beta0.old = beta0
   
    #expectation  
    
    prob = Functio_prob(beta0.old,lambda0.old,Ts,event,XB, Q)
    
    #maximization
    lambda0 = Function_lambda0 (prob,beta0,Ts,event,XB)
    estime = coxph_estimate (prob,Ts,event,XB,beta_ini, maxiter = 20)
    beta0 = estime$beta0
    beta_ini = beta0.old
    
    converge = (sqrt (sum ((beta0.old -beta0))^2 + sum ( (lambda0.old - lambda0)^2 ) ) < tol)
  
    if (it == maxits ) {
      cat("WARNING! NOT CONVERGENT FOR EM !", "\n")
      cat("maxits=", "\n")
      print(it)
      converge = FALSE
    }else if(is.na(beta0[1]) & is.na(beta0[2])){
      cat("WARNING! beta0 NOT AVAILABLE!", "\n")
      converge = FALSE
    }
    
    it = it + 1
  }
  
  return(list(beta0=beta0, lambda0=lambda0, prob=prob, converge= converge))
}


