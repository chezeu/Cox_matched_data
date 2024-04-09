
library(survival)
library(plyr)
library(dplyr)
library(rootSolve) 
library(nleqslv)

library(mvtnorm)
library(glmnet)
library(matrixStats)
library(clue)

#########  multinomiale Q

Matrix_function<-function(n,m,C,beta, surv_data){
  
  #surv_data <- Generate_data(m,n,beta)
  X <- surv_data[,3:ncol(surv_data)]
  XB <- as.matrix(X, ncol = p) #Only covariates
  
  ###### True-link data: The first n individual from B
  data_true <- surv_data[1:n,]
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  
  
  Q<-matrix(0,n,m)
  Lvec<-matrix(0,n,m)
  for (i in 1:n) {
    if(i!=1 & i!=n){ 
      pro = rep(0,m)
      pro[i] = C[1] #Probability of True links
      pro[i+1]= C[3]
      pro[i-1]=C[3]
      Q[i,]<-pro
      Lvec[i,] <- rmultinom(1, size = 1, pro = pro)} 
    if(i==1)  {
      pro = rep(0,m)
      pro[i] =  C[1] #Probability of True links
      pro[i+1]= C[2]
      pro[i+2]= C[3]
      Q[i,]<-pro
      Lvec[i,] <- rmultinom(1, size = 1, pro = pro)} 
    if(i==n){ 
      pro = rep(0,m)
      pro[i] = C[1] #Probability of True links
      pro[i-1]=C[2]
      pro[i-2]=C[3]
      Q[i,]<-pro
      Lvec[i,] <- rmultinom(1, size = 1, pro= pro) }  
  }
  Lvec
  X_naive <-Lvec %*% XB
  Z<- X_naive
  colnames(Z)<- c("Z1","Z2")
  # Naive data
  data_naive <- cbind( data_A, Z)
  trueLink_id <- diag(Lvec)
  
  return(list( data_true=data_true, data_A= data_A,XB=XB, data_naive=data_naive,Q=Q, Lvec= Lvec))
}

################## ##################################################


#1. Finding the risk set R(t) given some time t
GetRiskSet <- function(time_of_interest, time_vector, event_vector) {
  
  return(which(((time_vector == time_of_interest & event_vector == 1) | (time_vector > time_of_interest))))
  
}
#################
# equation naive

equa_naive <- function(beta, data_naive, XB, p) {
 
  n = nrow(data_naive)
  
  Ts <-as.matrix( data_naive[,1]) 
  event <- data_naive[,2] 
  Z <- as.matrix(data_naive[,3:(p+2)])
  X <- as.matrix(XB[,1:p])
  
  ezbeta <- exp(Z%*% beta)
  zezbeta <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p) 
  
  s<-0
  for(i in 1:n){ 
    
    Yi<- as.numeric( (Ts == Ts[i] ) | (Ts > Ts[i]))# at risk 
    risq <- which(Yi==1)
    if(length(risq)==1){  
      num_naive<- zezbeta[risq,]
      denum_naive<- ezbeta[risq]   
    }else if (length(risq)!=1){
      num_naive<-colSums( zezbeta[risq,])
      denum_naive<- sum(ezbeta[risq])  }  
    s<- s+ event[i]* ( Z[i,]- num_naive/denum_naive)
  }
  s
  return(H_naive=s) ##naive estimator 
}

# solve the equation 
coxph_equa_naive <- function(data_naive, XB, p,maxiter = 20){
  f <- function(x){
    equa_naive (beta=x, data_naive=data_naive, XB = XB, p = p )
  }
  
  fit_manual <- multiroot(f,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  iterations <- fit_manual$iter
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list(coef = beta, converge = converge,iterations=iterations))
}

library(simsalapar)
library(doParallel)


doOne2d_equat_estimat<- function(n,m,C,beta,surv_data){
  # Fixed parameters
  p = 2
  beta = matrix(c(0.5,-0.5), nrow = 2)
  
  # Generate data
  data <-Matrix_function(n,m,C,beta,surv_data)
  data_true <- data$data_true
  data_naive <- data$data_naive
  
  XB = data$XB
  #Q <- data$Q
  
 # Ts <-as.matrix( data_naive[,1]) 
  #event <- data_naive[,2] 
  #Z <- as.matrix(data_naive[,3:(p+2)])
 # X <- as.matrix(XB[,1:p])
  
 # ezbeta <- exp(Z%*% beta)
#  zezbeta <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
  
  # Theoretical estimating equation for true and naive data
  fit_true <- coxph(Surv(Time,delta)~.,data = data_true)
  coef_true <- as.vector(fit_true$coefficients)
  var_true <- diag(fit_true$var)
  
  fit_naive <- coxph(Surv(Time,delta)~.,data = data_naive)
  coef_naive <- as.vector(fit_naive$coefficients)
  var_naive <- diag(fit_naive$var)

  fit_naive2 <- coxph_equa_naive(data_naive, XB, p,maxiter = 20)
  coef_naive2 <- fit_naive2$coef
  converge_naive2 <- fit_naive2$converge
  
  return(list( coef_true=coef_true, coef_naive2=coef_naive2, converge_naive2=converge_naive2))
}
