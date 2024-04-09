
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
  
  return( Q=Q)
}

##################

#1. Finding the risk set R(t) given some time t
GetRiskSet <- function(time_of_interest, time_vector, event_vector) {
  
  return(which(((time_vector == time_of_interest & event_vector == 1) | (time_vector > time_of_interest))))
  
}
######################

#################

equa_W_sum1 <- function(beta,n,m,C,surv_data) {
  # Generate data
  
  XB <- surv_data[,3:ncol(surv_data)]
  X <- as.matrix(XB, ncol = p) #Only covariates
  
  data_true <- surv_data[1:n,]
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  Ts <-as.matrix( data_A[,1]) 
  event <- data_A[,2]
   
  Q <-Matrix_function_sum(n,m,C,beta)
  
  eXbeta <- exp(X%*% beta)
  XeXbeta <- X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  #######les sommes pour Z
  
  Z <- Q%*%X
  
  #### les sommes q*exp(beta x) pour tilfe_f
  
  som1<-  Q%*%eXbeta
  
  ###  les sommes pour tilde_g
  
  som2<- Q%*%XeXbeta
  
  
  dat1 <- cbind(Ts,Z)[which(event==1),]
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
      
      t2R <- sum(som1[risk] )
      t3R <- colSums(matrix(som2, ncol = p)[risk,] ) }
    
    s[i,] <- (Z1R - (t3R/t2R))
    
  }
  
  s <- colSums(s)
  
  return(s=s)
  
}

#################

equa_W_sum2 <- function(beta,n,m,C,surv_data) {
  
  # Generate data
  
  XB <- surv_data[,3:ncol(surv_data)]
  X <- as.matrix(XB, ncol = p) #Only covariates
  
  data_true <- surv_data[1:n,]
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  Ts <-as.matrix( data_A[,1]) 
  event <- data_A[,2]
  
  Q <-Matrix_function_sum(n,m,C,beta)
  
  eXbeta <- exp(X%*% beta)
  XeXbeta <- X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  
  
  # somme des deux valeurs max
  Z<-matrix(0, nrow = n, ncol = p)
  for (k in 1:n) {
    qk <- Q[k,]
    k1<- which(qk==max(qk))
    l1<-max(qk)
    qk[k1]<-0
    k2<-which(qk==max(qk))
    if(length(k2)!=1){k2<-k2[1]}
    qk[k1]<-l1
    qkx <- (qk*X)[c(k1,k2),]
    Z[k,]<- colSums(qkx)
  }  
  Z
  
  #ezbeta <- exp(Z%*% beta)
#  zezbeta <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
  
  #### les sommes q*exp(beta x) pour tilfe_f
  som1<- vector()
  for (j in 1:n) {
    qj <- Q[j,]
    j1<- which(qj==max(qj))
    l1<-max(qj)
    qj[j1]<-0
    j2<-which(qj==max(qj))
    if(length(j2)!=1){j2<-j2[1]}
    qj[j1]<-l1
    qjx <- (qj*eXbeta)[c(j1,j2),]
    som1[j] <- sum( qjx)
  }  
  som1 
  
  ###  les sommes pour tilde_g
  som2<-matrix(0, nrow = n, ncol = p)
  
  for (i in 1:n) {
    qi <- Q[i,]
  i1<- which(qi==max(qi))
    l1<-max(qi)
    qi[i1]<-0
    i2<-which(qi==max(qi))
    if(length(i2)!=1){i2<-i2[1]}
    qi[i1]<-l1
    qix <- (qi*XeXbeta)[c(i1,i2),]
    som2[i,] <- colSums( qix)
  }  
  som2
  
  dat1 <- cbind(Ts,Z)[which(event==1),]
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
      
      t2R <- sum(som1[risk] )
      t3R <- colSums(matrix(som2, ncol = p)[risk,] ) }
    
    s[i,] <- (Z1R - (t3R/t2R))
    
  }
  
  s <- colSums(s)
  return(H_w_sum2=s)  
}

# ###########solve the equations 
coxph_w_sum1 <- function(surv_data, p,maxiter = 20){
  f <- function(x){
    equa_W_sum1  (beta=x,n,m,C,surv_data)
  }
  
  fit_manual <- multiroot(f,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  iterations <- fit_manual$iter
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list(coef = beta, converge = converge,iterations=iterations))
}

coxph_w_sum2 <- function(surv_data, p,maxiter = 20){
  f <- function(x){
    equa_W_sum2  (beta=x,n,m,C,surv_data)
  }
  
  fit_manual <- multiroot(f,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  iterations <- fit_manual$iter
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list(coef = beta, converge = converge,iterations=iterations))
}
#############################
library(simsalapar)
library(doParallel)


doOne2d_w_sum1<- function(n,m,surv_data,C){
  # Fixed parameters
  p = 2
  beta = matrix(c(0.5,-0.5), nrow = 2)
  
  # surv_data<-Generate_data(m,n,beta)
  
  XB <- surv_data[,3:ncol(surv_data)]
  #X <- as.matrix(XB, ncol = p) #Only covariates
  
  data_true <- surv_data[1:n,]
  #data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  #Ts <-as.matrix( data_A[,1]) 
  #event <- data_A[,2]
  
  #Q <-Matrix_function(n,m,C,beta)
  
  #Z <- Q%*%X

  # Theoretical estimating equation for true and naive data
  fit_true <- coxph(Surv(Time,delta)~.,data = data_true)
  coef_true <- as.vector(fit_true$coefficients)
  var_true <- diag(fit_true$var)

  fit_w_sum1 <- coxph_w_sum1(surv_data, p,maxiter = 20)
  coef_w_sum1 <- fit_w_sum1$coef
  converge_w_sum1 <- fit_w_sum1$converge
  
  return(list( coef_true=coef_true, coef_w_sum1=coef_w_sum1, converge_w_sum1=converge_w_sum1))
}


doOne2d_w_sum2<- function(n,m,surv_data,C){
  p = 2
  beta = matrix(c(0.5,-0.5), nrow = 2)
  
  # surv_data<-Generate_data(m,n,beta)
  
  XB <- surv_data[,3:ncol(surv_data)]
  #X <- as.matrix(XB, ncol = p) #Only covariates
  
  data_true <- surv_data[1:n,]
  #data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
 # Ts <-as.matrix( data_A[,1]) 
  #event <- data_A[,2]

  
  #Q <-Matrix_function_sum(n,m,C,beta)

  # Theoretical estimating equation for true and naive data
  fit_true <- coxph(Surv(Time,delta)~.,data = data_true)
  coef_true <- as.vector(fit_true$coefficients)
  var_true <- diag(fit_true$var)
  
  fit_w_sum2 <- coxph_w_sum2(surv_data, p,maxiter = 20)
  coef_w_sum2 <- fit_w_sum2$coef
  converge_w_sum2 <- fit_w_sum2$converge
  
  return(list( coef_true=coef_true, coef_w_sum2=coef_w_sum2, converge_w_sum2=converge_w_sum2))
}
