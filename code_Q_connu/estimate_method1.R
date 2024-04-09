library(survival)
library(plyr)
library(dplyr)
library(rootSolve) 
library(nleqslv)

library(mvtnorm)
library(glmnet)
library(matrixStats)
library(clue)

library(clue)
library(matrixStats)
library(klaR)
library(ludic)

library(doParallel)
library(RecordLinkage)

#########  multinomiale Q

Matrix_function_sum<-function(n,m,C,beta){

  Q<-matrix(0,n,m)
  #Lvec<-matrix(0,n,m)
  for (i in 1:n) {
    if(i!=1 & i!=n){ 
      pro = rep(0,m)
      pro[i] = C[1] #Probability of True links
      pro[i+1]= C[3]
      pro[i-1]=C[3]
      Q[i,]<-pro
      #Lvec[i,] <- rmultinom(1, size = 1, pro = pro)
    } 
    if(i==1)  {
      pro = rep(0,m)
      pro[i] =  C[1] #Probability of True links
      pro[i+1]= C[2]
      pro[i+2]= C[3]
      Q[i,]<-pro
      #Lvec[i,] <- rmultinom(1, size = 1, pro = pro)
    } 
    if(i==n){ 
      pro = rep(0,m)
      pro[i] = C[1] #Probability of True links
      pro[i-1]=C[2]
      pro[i-2]=C[3]
      Q[i,]<-pro
      # Lvec[i,] <- rmultinom(1, size = 1, pro= pro)
    }  
  }
  
  return(Q=Q)
}

##################

#1. Finding the risk set R(t) given some time t
GetRiskSet <- function(time_of_interest, time_vector, event_vector) {
  
  return(which(((time_vector == time_of_interest & event_vector == 1) | (time_vector > time_of_interest))))
  
}
######################
# ####################initialisation########################

#cumulative function
Funct_lambda2<-function(lambda0,surv_data){
  #lambda0<- rep(0.5,length(event))
 
  data_true <- surv_data[1:n,]
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  Ts <-as.matrix( data_A[,1]) 
  event <- data_A[,2]
  
  l<- which(event==0)
  lambda0[l]<-0
  
  observe<- function(time_of_interest, time_vector, event_vector) {
    return(which( (time_vector <= time_of_interest) ))
  }
  
  lambda2<-vector()
  for (i in 1:n) {
    vi<- observe(Ts[i], Ts, event)
    lambda2[i]<- sum(lambda0[vi])   
  }
  return(lambda2=lambda2)
}

### proba aposteriorie

Functio_prob<-function(beta0,lambda0,surv_data){
  #beta0<- c(0.1,0.1) 
  #lambda0<- rep(0.5,length(event))
  
lambda2<-Funct_lambda2(lambda0,surv_data)

Q <-Matrix_function_sum(n,m,C,beta)

XB <- surv_data[,3:ncol(surv_data)]
X <- as.matrix(XB, ncol = p) #Only covariates

data_true <- surv_data[1:n,]
data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
Ts <-as.matrix( data_A[,1]) 
event <- data_A[,2]


eXbeta0 <- exp(X%*% beta0)
XeXbeta0 <- X*matrix(rep(exp(X%*% beta0),p), ncol = p)


#pi(ij)

prob<-matrix(0,n,m)
for (i in 1:n) {
  numerateur<-vector()
  for (j in 1:m) {
    numerateur[j]<- Q[i,j] * ( lambda0[i]*eXbeta0[j])^event[i] * exp(-lambda2[i]* eXbeta0[j] )
  }
  denominateur<-sum(numerateur)
  prob[i,]<- numerateur/denominateur
}
prob
return( list(eXbeta0=eXbeta0,XeXbeta0=XeXbeta0,prob=prob))
}

################ estimations ###################################

#lambda0
Function_lambda0<-function(beta0,lambda0,surv_data){
  
  init<-Functio_prob(beta0,lambda0,surv_data)
  prob<- init$prob
  eXbeta0<-init$eXbeta0
  
  data_true <- surv_data[1:n,]
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  Ts <-as.matrix( data_A[,1]) 
  event <- data_A[,2]
  
  som1<-prob%*%eXbeta0
  
  s <- vector()
  for (i in 1:n) {
    
    ts <- Ts[i]
    risk <- GetRiskSet(ts, Ts, event)
    nrisk <- length(risk)
    if(nrisk==0){ 
      t2R<-1
    }else if(nrisk==1){
      t2R <- som1[risk] 
    }else if(nrisk!=1){
      
      t2R <- sum( som1[risk] )
    }
    
    s[i] <- (1/t2R)*event[i]
  }
  
  return(s=s)
}

# equation estimente ########################
equa_estimate <- function(beta0,lambda0,surv_data) {

  # Generate data
  XB <- surv_data[,3:ncol(surv_data)]
  X <- as.matrix(XB, ncol = p) #Only covariates
  
  data_true <- surv_data[1:n,]
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  Ts <-as.matrix( data_A[,1]) 
  event <- data_A[,2]
  
  init<-Functio_prob(beta0,lambda0,surv_data)
  prob<- init$prob
  eXbeta0<-init$eXbeta0
  XeXbeta0<-init$XeXbeta0
  
  #######les sommes pour Z

  Z <- prob%*%X 

  #### les sommes q*exp(beta x) pour tilfe_f
  
  som1<- prob%*%eXbeta0

  ###  les sommes pour tilde_g
  
  som2<-prob%*%XeXbeta0

  dat1 <- cbind(Ts, Z)[which(event==1),]
  
  ## Estimating equation
  s <- matrix(0,nrow=nrow(dat1), ncol = p)
  for (i in 1:nrow(dat1)) {
    ts <- dat1[i, 1]
    Z1R <- dat1[i,2:ncol(dat1)]
    
    risk <- GetRiskSet(ts, Ts, event)
    nrisk <- length(risk)
    if(nrisk==1){
      t2R <- som1[risk] 
      t3R <- matrix(som2, ncol = p)[risk,] 
      
    }else if(nrisk!=1){
      
      t2R <- sum( som1[risk] )
      t3R <- colSums(matrix(som2, ncol = p)[risk,] ) }
    
    s[i,] <- (Z1R - t3R/t2R)
    
  }
  
  s <- colSums(s)
  
  return(s=s)
  
}

# ###########solve the equation 
coxph_estimate<- function(beta0,lambda0,surv_data,maxiter = 50){
  f <- function(x){
    equa_estimate  (beta0=x,lambda0,surv_data)
  }
  fit_manual <- nleqslv( c(0,0),f, method = c("Broyden", "Newton"))
  beta0 <- fit_manual$x
  iterations <- fit_manual$iter
  converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1) )
  
  #xstart <- matrix(rnorm(20,0,1), ncol = 2)
  #Zero <-  searchZeros(xstart,f)
  #fit_manual <- multiroot(f,start = rep(0,p), maxiter = maxiter)
 # beta0 <- fit_manual$root
  #iterations <- fit_manual$iter
#converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list( beta0 = beta0,converge =converge) )
}



###iterrations###########################################

#valeurs initials
Func_itteration<-function(beta0,lambda0,surv_data,tol= 1e-6,maxits = 500){
  
  #set.seed (23051986)
  #surv_data <- Generate_data(m,n,beta)
  
  # Generate data
  #XB <- surv_data[,3:ncol(surv_data)]
  #X <- as.matrix(XB, ncol = p) #Only covariates
  
  #data_true <- surv_data[1:n,]
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  Ts <-as.matrix( data_A[,1]) 
  event <- data_A[,2]
  
  it = 0
  converge = FALSE
  
while ((!converge) & (it < maxits)){ 

  #expectation  

  lambda2<-Funct_lambda2(lambda0,surv_data)
init<-Functio_prob(beta0,lambda0,surv_data)
prob<- init$prob

lambda0.old<-lambda0
beta0.old<-beta0

lambda0<-Function_lambda0(beta0=beta0.old,lambda0=lambda0.old,surv_data)
#maximization
estime<-coxph_estimate(beta0=beta0.old,lambda0=lambda0.old,surv_data,maxiter =50)
beta0<-estime$beta0

it = it + 1
converge <-  (abs(beta0.old -beta0)/ (beta0.old+0.01))[1] < tol &&
                (abs(beta0.old -beta0)/ (beta0.old+0.01))[2] < tol

if (it == maxits) {
  cat("WARNING! NOT CONVERGENT!", "\n")
  converge = FALSE
}else if(is.na(beta0[1]) & is.na(beta0[2])){
  cat("WARNING! beta0 NOT AVAILABLE!", "\n")
  converge = FALSE
}

}
  
return(list(beta0=beta0, lambda0=lambda0,prob=prob, converge= converge))
}
#############################################################
library(simsalapar)
library(doParallel)

doOne2d_estimate<- function(n,m,C,beta0,lambda0,surv_data){
  
  #surv_data <- Generate_data(m,n,beta)
  data_true <- surv_data[1:n,]
  
  # Theoretical estimating equation for true data
  fit_true <- coxph(Surv(Time,delta)~.,data = data_true)
  coef_true <- as.vector(fit_true$coefficients)
  var_true <- diag(fit_true$var)
  
  fit_estimate  <- Func_itteration(beta0,lambda0,surv_data,tol= 1e-6,maxits = 100)
  coef_estimate  <- fit_estimate$beta0
  lambda0_estimate<-fit_estimate$lambda0
  converge_estimate <- fit_estimate$converge
  
  return(list( coef_true=coef_true, coef_estimate =coef_estimate,
               converge_estimate =converge_estimate,lambda0_estimate=lambda0_estimate ))
}
