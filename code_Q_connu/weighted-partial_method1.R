
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

Matrix_function<-function(n,m,C,beta){
  
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

#################

equa_w_partial  <- function(beta, surv_data) {
    # Generate data
 # surv_data<-Generate_data(m,n,beta)
  
    XB <- surv_data[,3:ncol(surv_data)]
    X <- as.matrix(XB, ncol = p) #Only covariates
    
    data_true <- surv_data[1:n,]
    data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
    Ts <-as.matrix( data_A[,1]) 
    event <- data_A[,2]
    Q <-Matrix_function(n,m,C,beta)
    
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
    h<-which(Q[i,]!=0)
    pi<-min(Q[i,h])
    
    risk <- GetRiskSet(ts, Ts, event)
    nrisk <- length(risk)
    if(nrisk==1){
      t2R <- som1[risk] 
      t3R <- matrix(som2, ncol = p)[risk,] 
      
    }else if(nrisk!=1){
      l<-which(Q[i,][risk]!=0)
      pk<-  min(Q[i,l] )
      t2R <- sum(som1[risk]/pk )
      t3R <- colSums(matrix(som2, ncol = p)[risk,]/pk ) 
      }
    
    s[i,] <- (1/pi)*(Z1R - (t3R/t2R))
    
  }
  
  s <- colSums(s)
  
  return(s=s)
  
}

#################

# ###########solve the equation 
coxph_w_partial <- function(surv_data, p,maxiter = 20){
  f <- function(x){
    equa_w_partial(beta=x, surv_data)
  }
  #fit_manual <- nleqslv( c(0,0),f, method = c("Broyden", "Newton"))
  #beta0 <- fit_manual$x
  #iterations <- fit_manual$iter
  #converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1) )
 
  #xstart <- matrix(rnorm(20,0,1), ncol = 2)
  #Zero <-  searchZeros(xstart,f)
  fit_manual <- multiroot(f,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  iterations <- fit_manual$iter
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list(coef = beta, converge = converge,iterations=iterations))
}

library(simsalapar)
library(doParallel)


doOne2d_w_partial <- function(n,m,surv_data,C){

  # surv_data<-Generate_data(m,n,beta)

  data_true <- surv_data[1:n,]
 
   # Theoretical estimating equation for true and naive data
  fit_true <- coxph(Surv(Time,delta)~.,data = data_true)
  coef_true <- as.vector(fit_true$coefficients)
  var_true <- diag(fit_true$var)
  
  fit_w_partial  <- coxph_w_partial(surv_data, p,maxiter = 20)
  coef_w_partial  <- fit_w_partial$coef
  converge_w_partial  <- fit_w_partial$converge
  
  return(list( coef_true=coef_true, coef_w_partial =coef_w_partial ,
               converge_w_partial =converge_w_partial ))
}


#applications
nsim=10

beta=c(0.5,-0.5)
C=c(0.7,0.1,0.1)
m=70
n=50
p=2

coef_w_partial <- matrix(0,nrow = nsim, ncol = p)
coef_true<-matrix(0,nrow = nsim, ncol = p)
converge_w_partial  <- vector()

test_w_partial <- function(nsim, C){
  for (i in 1:nsim) {
    
    surv_data<-Generate_data(m,n,beta)
    g <- doOne2d_w_partial  (n,m,surv_data,C)
    
    coef_true[i,] <- g$coef_true
    coef_w_partial [i,] <- g$coef_w_partial 
    converge_w_partial [i] <- g$converge_w_partial 
  }
  
  L2 <- which(converge_w_partial !=0)
  coef_w_partial<- coef_w_partial [L2,]
  
  return(list( coef_w_partial = coef_w_partial , coef_true= coef_true))
}


#results

result_w_partial  <- test_w_partial (nsim,C)
coef_w_sum1<-colMeans(result_w_partial $coef_w_partial )

coef_true<-colMeans(result_w_sum1$coef_true)


#################################
