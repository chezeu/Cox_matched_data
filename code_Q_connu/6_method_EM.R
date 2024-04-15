

###1. Finding the risk set R(t) given some time t

GetRiskSet <- function(time_of_interest, time_vector, event_vector) {
  
  return(which(((time_vector == time_of_interest & event_vector == 1) | (time_vector > time_of_interest))))
  
}


#####
#number of observations before time t

observe<- function(time_of_interest, time_vector, event_vector) {
  return(which( (time_vector <= time_of_interest) ))
}
####

#cumulative function

Funct_lambda2<-function(lambda0,Ts,event){
  #lambda0<- rep(0.5,length(event))
  
  n = length(Ts)
  l =  which(event==0)
  lambda0[l] = 0
  
  lambda2 = vector()
  for (i in 1:n) {
    vi = observe(Ts[i], Ts, event)
    lambda2[i] =  sum(lambda0[vi])   
  }
  return(lambda2=lambda2)
}
########
### proba aposteriorie (pi)

Functio_prob<-function(beta0,lambda0,Ts,event,XB, Q){
  #beta0<- c(0.1,0.1) 
  #lambda0<- rep(0.5,length(event))
  
  lambda2 = Funct_lambda2(lambda0,Ts,event)
  X = as.matrix(XB, ncol = p) # covariates
  
  eXbeta0 = exp(X%*% beta0)
  XeXbeta0 = X*matrix(rep(exp(X%*% beta0),p), ncol = p)
  
  n = length(Ts)
  m = nrow(XB)
  
  #pi(ij)
  
  prob = matrix(0,n,m)
  for (i in 1:n) {
    
    numerateur = vector()
    for (j in 1:m) {
      numerateur[j] = Q[i,j] * ( lambda0[i]*eXbeta0[j])^event[i] * exp(-lambda2[i]* eXbeta0[j] )
    }
    
    denominateur = sum(numerateur)
    
    prob[i,] = numerateur/denominateur
  }
  prob
  return( prob = prob)
}
#############
################ estimations ###################################

#lambda0
Function_lambda0<-function(prob,beta0,Ts,event,XB, Q){
  
  p = ncol(XB)
  n = length(Ts)
  
  X = as.matrix(XB, ncol = p) # covariates
  eXbeta0 = exp(X%*% beta0)
  
  som1 = prob%*%eXbeta0 # somme sur les j
  
  s = vector()
  for (i in 1:n) {
    
    Yi = as.numeric( Ts >= Ts[i])# at risk
    risq = which(Yi==1)
    nrisk = length(Yi)
    
    if(nrisk==1){
      t2R = som1[risq] 
    }else if(nrisk!=1){
      
      t2R = sum( som1[risq] )
    }
    
    s[i] = (1/t2R)*event[i]
  }
  
  return( lambda0_1 = s)
}

######################
#  estimating equation ########################

equa_estimate <- function(beta,prob,Ts,event,XB, Q) {
  
  p = ncol(XB)
  n = length(Ts)
  
  X = as.matrix(XB, ncol = p) # covariates
  eXbeta = exp(X%*% beta)
  XeXbeta = X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  #######les sommes pour Z
  
  Z = prob%*%X 
  
  #### les sommes q*exp(beta x) pour tilfe_f
  
  som1 = prob%*%eXbeta
  
  ###  les sommes pour tilde_g
  
  som2 = prob%*%XeXbeta
  
  dat1 = cbind(Ts, Z)[which(event==1),]
  
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
      
      t2R = sum( som1[risk] )
      t3R = colSums(matrix(som2, ncol = p)[risk,] ) }
    
    s[i,] = (Z1R - t3R/t2R)
    
  }
  
  s <- colSums(s)
  
  return(s=s)
  
}

#####################
# ###########solve the equation 

coxph_estimate<- function(beta,prob,Ts,event,XB, Q,beta_ini,maxiter = 50){
  
  f <- function(x){
    equa_estimate  (beta= x,prob,Ts,event,XB, Q)
  }
  fit_manual <- nleqslv( c(beta_ini[1],beta_ini[2]),f, method = c("Broyden", "Newton"))
  beta0 <- fit_manual$x
  iterations <- fit_manual$iter
  converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1) )
  
  #xstart <- matrix(rnorm(20,0,1), ncol = 2)
  #Zero <-  searchZeros(xstart,f)
  #fit_manual <- multiroot(f,start = rep(0,p), maxiter = maxiter)
  # beta0 <- fit_manual$root
  #iterations <- fit_manual$iter
  #converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list( beta0 = beta0, converge = converge, iterations = iterations) )
}
################

### iterrations ###########################################

#valeurs initials

Func_itteration<-function(beta,beta0,lambda0,Ts,event,XB, Q,tol= 1e-6, maxits = 500){
  
  p = ncol(XB)
  n = length(Ts)
  
  it = 0
  converge = FALSE
  
  while ((!converge) & (it < maxits)){ 
    
    lambda0.old = lambda0
    beta0.old = beta0
    beta_ini = beta0.old
    
    #expectation  
    
    lambda2 = Funct_lambda2(lambda0.old,Ts,event)
    
    prob = Functio_prob(beta0.old,lambda0.old,Ts,event,XB, Q)
    
    
    #maximization
    lambda0 = Function_lambda0 (prob,beta0, Ts,event,XB, Q)
    estime = coxph_estimate (beta,prob,Ts,event,XB, Q,beta_ini, maxiter = 50)
    beta0 = estime$beta0
    
    converge = (sum ((beta0.old -beta0))^2 + sum ( (lambda0.old - lambda0)^2 ) ) < tol  
    
    
    if (it == maxits) {
      cat("WARNING! NOT CONVERGENT!", "\n")
      converge = FALSE
    }else if(is.na(beta0[1]) & is.na(beta0[2])){
      cat("WARNING! beta0 NOT AVAILABLE!", "\n")
      converge = FALSE
    }
    
    it = it + 1
    
  }
  
  return(list(beta0=beta0, lambda0=lambda0, prob=prob, converge= converge))
}

