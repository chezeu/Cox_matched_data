


setwd( "C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")
source("1_data_generate.R")


################### cumulative function
surv_data = Generate_data(m,n,beta)
mf = Matrix_function(n,m,C,beta,surv_data)
data_true = mf$data_true
data_naive = mf$data_naive

Ts = mf$data_true$Time
event = mf$data_true$delta
Q = mf$Q
XB = mf$XB

alpha0=0.1
gamma0=0.1

Funct_lambda<-function(gamma0,alpha0,Ts){
  #lambda0<- rep(0.1,length(event))
  
  n = length(Ts)
  lambda0 = vector()
  lambda2 = vector()
  for (i in 1:n) {
    lambda0[i] = (alpha0*(Ts[i])^(alpha0 -1))/( gamma0^alpha0)
    lambda2[i] = (Ts[i]/gamma0)^alpha0
  }
  return(list(lambda0= lambda0,lambda2=lambda2))
}

sf=Funct_lambda(gamma0,alpha0,Ts)
lambda0=sf$lambda0
lambda2=sf$lambda2
########
### proba aposteriories (pi)

Functio_prob<-function(beta0,lambda0,lambda2,Ts,event,XB,Q){
  #beta0<- c(0.1,0.1) 
  #alpha0<- 0.1
  #gamma0<- 0.1
  
  X = as.matrix(XB, ncol = p) # covariates
  eXbeta0 = exp(X%*% beta0)
  
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

prob= Functio_prob(beta0,lambda0,lambda2,Ts,event,XB,Q)
#############

######################
#  estimating equation ########################

s = c(beta,gamma0,alpha0)
solv_equations <- function(s, prob,Ts,event,XB) {
  
  p = ncol(XB)
  n = length(Ts)
  beta = s[1:2]
  gamma0 = s[3]
  alpha0= s[4]
  X = as.matrix(XB, ncol = p) # covariates
  eXbeta = exp(X%*% beta)
  XeXbeta = X*matrix(rep(exp(X%*% beta),p), ncol = p)
  

  D = sum(event)
  Z = prob%*%X 
  
  #### values for q*exp(beta x) pour tilfe_f
  
  som1 = prob%*%eXbeta
  
  ###  values for q*x*exp(beta x)
  
  som2 = prob%*%XeXbeta
  
  dat1 = cbind(Ts, Z)[which(event==1),]
  
  ## Estimating equations
  
  L_1 = (D/alpha0)+ sum(event*log(Ts/gamma0))- sum(som1*log(Ts/gamma0)*((Ts/gamma0))^alpha0)
  
  L_2 = -( (D* alpha0)/gamma0) + ((alpha0/ gamma0)* sum( som1* (Ts/gamma0)^alpha0))
  
  L_3 =colSums (Z*event)-colSums (som2* (Ts/gamma0)^alpha0 )
  
  s = c(L_1,L_2,L_3 )
  
  return(s=s)
  
}
vi= solv_equations (s, prob,Ts,event,XB)
#####################
# ###########solve the equation 

coxph_estimate2<- function(prob,Ts,event,XB,maxiter = 20){
  s=c(beta,gamma0,alpha0)
  f <- function(x){
    solv_equations (s, prob,Ts,event,XB)
  }
   fit_manual <- nleqslv( c(0,0,0,0),f, method = c("Broyden", "Newton"))
   beta0 <- fit_manual$x
  iterations <- fit_manual$iter
   converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1) )
  
  #xstart <- matrix(rnorm(20,0,1), ncol = 2)
  #Zero <-  searchZeros(xstart,f)
  #fit_manual = multiroot(f,start = c(0,0,0,0), maxiter = maxiter)
  #beta0 = fit_manual$root
  #iterations =  fit_manual$iter
  #converge = as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list( beta0 = beta0, converge = converge, iterations = iterations) )
}

fit_estimate = coxph_estimate2(prob,Ts,event,XB,maxiter = 20)
coef_estimate_s[i,] = fit_estimate$beta0
converge_estimate[i] = as.numeric(fit_estimate$converge)

################

### iterrations ###########################################

#valeurs initials

Func_itteration<-function(beta0,lambda0,Ts,event,XB, Q,tol= 1e-6, maxits = 500){
  
  p = ncol(XB)
  n = length(Ts)
  beta_ini = rep(0, p)
  it = 0
  converge = FALSE
  
  while ((!converge) & (it < maxits)){ 
    
    lambda0.old = lambda0
    beta0.old = beta0
    
    #expectation  
    
    lambda2 = Funct_lambda2(lambda0.old,Ts,event)
    
    prob = Functio_prob(beta0.old,lambda0.old,lambda2,Ts,event,XB, Q)
    
    
    #maximization
    lambda0 = Function_lambda0 (prob,beta0.old,Ts,event,XB)
    estime = coxph_estimate (prob,Ts,event,XB,beta_ini,maxiter = 20)
    beta0 = estime$beta0
    beta_ini = beta0.old
    
    converge = sqrt (sum ((beta0.old -beta0))^2 + sum ( (lambda0.old - lambda0)^2 ) ) < tol  
    
    
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




death_times= sort(unique(Ts*event))
death_times=death_times[death_times!=0]



