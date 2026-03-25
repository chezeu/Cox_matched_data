


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
 # X = XB # covariates
  
  eXbeta0 = as.vector(exp(XB %*% beta0))
  
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
Function_lambda0 <- function(prob, beta0, Ts, event, XB, RiskSets){
  
  XB <- as.matrix(XB)
  n  <- length(Ts)
  
  eXbeta0 = as.vector(exp(XB %*% beta0))
 # som1 = as.vector(prob %*% eXbeta0)
  
  prob_sparse <- Matrix(prob, sparse = TRUE) #avoid som1
  
  lambda0_1 <- numeric(n)
  
  for (i in which(event==1)) {
    
      #risq <- RiskSets[[i]]   #   plus de which()
      
    #  t2R <- sum(som1[risq])
      
    risk <- RiskSets[[i]]
    risk_vec <- numeric(n)
    risk_vec[risk] <- 1
    
    pi_sum <- as.numeric(crossprod(prob_sparse, risk_vec))
    t2R <- sum(pi_sum * eXbeta0)
    
     if (!is.finite(t2R) || t2R <= 0)
     
        t2R <- .Machine$double.eps
      
     
     lambda0_1[i] = 1 / t2R
  }
  
  return(lambda0_1)
}

equa_estimate <- function(beta, prob, Ts, event, XB, RiskSets) {
  
  XB <- as.matrix(XB)
  n  <- length(Ts)
  p  <- ncol(XB)
  
  eXbeta = as.vector(exp(XB %*% beta))
  #XeXbeta = XB * matrix(rep(eXbeta, p), ncol = p)
  XeXbeta <- XB * eXbeta
  

  # Z = prob%*% XB
    prob_sparse <- Matrix(prob, sparse = TRUE)
   Z <- prob_sparse %*% XB
   
 # som1 = as.vector(prob_sparse %*% eXbeta)
 # som2 = prob_sparse %*% XeXbeta
  
  s <- numeric(p)
  
  for (i in which(event==1)) {
 
      
      risk <- RiskSets[[i]]   #  plus de which()
      
      pi_sum <- colSums(prob_sparse[risk, , drop = FALSE])
      t2R <- sum(pi_sum * eXbeta) # AVOID SUM1
      
      if (!is.finite(t2R) || t2R <= 0)
        t2R <- .Machine$double.eps
      
      #t2R <- sum(som1[risk])
      
      t3R <- as.numeric(pi_sum %*% XeXbeta)# avoid sum2
      
     # t3R <- apply(som2[risk,],2,sum)
      s <- s + (Z[i,] - t3R/t2R)
    }
  
  
  return(as.numeric(s))
}


coxph_estimate <- function(prob, Ts, event, XB,
                           beta_ini, RiskSets, maxiter = 50){
  
  f <- function(x){
    equa_estimate(x, prob, Ts, event, XB, RiskSets)
  }
  
  fit_manual <- nleqslv(beta_ini, f,
                        method =  "Newton") 
  
  beta0 <- as.numeric(fit_manual$x)
  iterations <- fit_manual$iter
  converge <- as.numeric(iterations < maxiter)
  
  if (iterations  == maxiter) {
    cat("WARNING! NOT CONVERGENT FOR NEWTON: NEED A GOOD INITIAL VALUE!", "\n")
    converge = FALSE
  }
  return(list(beta0=beta0,
              converge=converge,
              iterations=iterations))
}

 
Func_itteration <- function(beta0, lambda0, Ts, event, XB, Q, 
                            tol=1e-6, maxits=500){
  
  XB <- as.matrix(XB)
  p  <- ncol(XB)
  n  <- length(Ts)
  
  ##  PRECOMPUTE RISK SETS ONCE
  RiskSets <- lapply(1:n, function(i) which(Ts >= Ts[i]))
  
  beta_ini = rep(0, p)
  it = 0
  converge = FALSE
  
  while ((!converge) && (it < maxits)){ 
    
    lambda0_old = lambda0
    beta0_old   = beta0
    
    # E-step
    lambda2 = Funct_lambda2(lambda0_old, Ts)
    
    prob = Functio_prob(beta0_old, lambda2, Ts, event, XB, Q)
    
    # M-step
    lambda0 = Function_lambda0(prob, beta0_old, Ts, event, XB, RiskSets)
    
    estime = coxph_estimate(prob, Ts, event, XB, 
                            beta_ini = beta0_old,
                            RiskSets = RiskSets,
                            maxiter = 50)
    
    beta0 = estime$beta0
    
    d_beta   <- sqrt(sum((beta0 - beta0_old)^2))
    d_lambda <- sqrt(sum((lambda0 - lambda0_old)^2))
    
    converge <- (d_beta + d_lambda) < tol
    
    it = it + 1
  }
  
  return(list(beta0=beta0,
              lambda0=lambda0,
              prob=prob,
              converge=converge,
              it=it))
}
