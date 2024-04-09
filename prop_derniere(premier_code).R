
library(survival)
library(plyr)
library(dplyr)
library(rootSolve) 
library(nleqslv)

library(mvtnorm)
library(glmnet)
library(matrixStats)
library(clue)

## generer les bases de donn?es######################################################
Generate_data <- function(m,n,rho,beta){
  # set.seed(42)
  ##covarites data
  X1 <- rnorm(m,0,1)
  X2 <- rbinom(m,size = 1, prob = 0.7)
  X <- as.matrix(cbind(X1, X2))
  
  U <- runif(m,0,1) # Time to event variable
  Tt = -log(U)/exp((X%*%beta))  #lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = rho)
  c=1.777907
  
  Time <- pmin(Tt,c)
  delta <- as.numeric(Tt<=c)
  surv_data <- data.frame(Time,delta,X)
  colnames(surv_data)=c('Time','delta',paste("X", 1:p, sep = ""))
  
 # true_data <- surv_data[1:n,]
  
 # database_A <- true_data[,1:2]  #base de donnee de A
 # XB <- surv_data[,3:dim(surv_data)[2]] #base de donnee de covariables
  return(surv_data=surv_data )
}

#########  multinomiale Q

Matrix_function<-function(n,m,C,rho,beta){
  
  surv_data <- Generate_data(m,n,rho,beta)
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

equa_estimat <- function(beta, data_naive, XB, p,C) {
  
  Ts <-as.matrix( data_naive[,1]) 
  event <- data_naive[,2] 
  Z <- as.matrix(data_naive[,3:(p+2)])
  X <- as.matrix(XB[,1:p])
  
  ezbeta <- exp(Z%*% beta)
  eXbeta<-exp(X%*% beta)
  
  zezbeta <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
  XeXbeta<- X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  donnees <-Matrix_function(n,m,C,rho,beta)
  Q <- donnees$Q
  
  
  #######les sommes pour tilde_x
  
  som0 <- matrix(0,nrow=n, ncol = p)
  for(i in 1:n){ 
    qi <- Q[i,]
    qix <- (qi*X)[-i,]
    som0[i,] <- colSums(qix)
  }
  som0
  diagonal <- diag(Q)
  tilde_x <- (Z/diag(Q))-som0/diag(Q)   
  
  
  #### les sommes q*exp(beta x) pour tilfe_f
  som1<- vector()
  for (j in 1:n) {
    qj <- Q[j,]
    qjx <- (qj*eXbeta)[-j]
    som1[j] <- sum( qjx)
  }  
  som1 
  
  
  ###  les sommes pour tilde_g
  
  som2<-matrix(0, nrow = n, ncol = p)
  for (k in 1:n) {
    qk <- Q[k,]
    qkx <- (qk*XeXbeta)[-k,]
    som2[k,]<- colSums(qkx)
  }  
  som2
  
  dat1 <- cbind(Ts,diagonal, tilde_x)[which(event==1),]
  ## Estimating equation
  s <- matrix(0,nrow=nrow(dat1), ncol = p)
  for (i in 1:nrow(dat1)) {
    ts <- dat1[i, 1]
    qi <- dat1[i,2]
    t1R <- dat1[i,3:ncol(dat1)]
    
    risk <- GetRiskSet(ts, Ts, event)
    nrisk <- length(risk)
    if(nrisk==1){
      t2R <- (ezbeta -som1)[risk] 
      t3R <- matrix(  zezbeta - som2, ncol = p)[risk,] 
      
    }else if(nrisk!=1){
      
      t2R <- sum( (ezbeta -som1)[risk] )
      t3R <- colSums(matrix(  zezbeta - som2, ncol = p)[risk,] ) }
    
    s[i,] <- (t1R - t3R / t2R)
    
  }
  
  s <- colSums(s)
  
  return(s=s)
  
}


#2. Proposed estimating equations of huan
aee1 <- function(beta, data_naive, XB, p, alpha, Q=1) {
  N = nrow(XB)
  n = nrow(data_naive)
  
  Ts <- data_naive[,1] 
  event <- data_naive[,2] 
  Z <- as.matrix(data_naive[,3:(p+2)])
  X <- as.matrix(XB[,1:p])
  
  muZ <- exp(Z%*% beta)
  psiZ <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
  
  
  
  if (Q == 1){
    mean.muX <- mean(exp(X%*%beta)) 
    mean.psiX <- colMeans((X*matrix(rep(exp(X %*% beta),p),ncol =p)))
    
    temp1 <- Z/alpha - (1/alpha-1)*matrix(rep(colMeans(X),n), nrow = n, byrow = TRUE)
    
    dat <- cbind(Ts, temp1)[which(event==1),]
    ## Estimating equation
    ee <- apply(dat, 1, 
                function(df){
                  df <- matrix(df, nrow = 1)
                  ts <- df[, 1]
                  t1R <- df[, 2:ncol(df)]
                  
                  risk <- GetRiskSet(ts, Ts, event)
                  nrisk <- length(risk)
                  
                  t2R <- sum(muZ[risk])-(1-alpha)*nrisk*mean.muX
                  t3R <- colSums(matrix(psiZ[risk,], ncol = p)) - (1-alpha)*nrisk*mean.psiX
                  
                  return((t1R - t3R / t2R))
                  
                }) %>% 
      matrix(nrow = p)%>%
      rowSums()
  }
  else if(Q>1){
    
    alpha.Q <- matrix(0, ncol = 1, nrow = n)
    meanX.Q <- matrix(0, ncol = p, nrow = n)
    mean.muX.Q <- matrix(0, ncol = 1, nrow = n)
    mean.psiX.Q <- matrix(0, ncol = p, nrow = n)
    
    for (q in 1:Q) {
      alpha.Q[which(data_naive$block == q),] <- alpha[q]
      
      X.q <- X[which(XB$block == q),1:p]
      n.q <- length(which(data_naive$block == q))
      
      meanX.Q[which(data_naive$block == q),] <- matrix(rep(colMeans(X.q),n.q), nrow = n.q, byrow = TRUE)
      mean.muX.Q[which(data_naive$block == q),] <- mean(exp(X.q%*%beta))
      mean.psiX.Q[which(data_naive$block == q),] <- matrix(rep(colMeans((X.q*matrix(rep(exp(X.q %*% beta),p),ncol =p))),n.q), nrow = n.q, byrow = TRUE)
    }
    
    temp1 <- Z/matrix(rep(alpha.Q,p),ncol = p) - matrix(rep(1/alpha.Q-1,p),ncol = p)*meanX.Q
    
    temp2 <- muZ/alpha.Q - (1/alpha.Q-1)*mean.muX.Q
    
    temp3 <- psiZ/matrix(rep(alpha.Q,p),ncol = p) - matrix(rep(1/alpha.Q-1,p),ncol = p)*mean.psiX.Q
    
    dat <- cbind(Ts, temp1)[which(event==1),]
    ## Estimating equation
    ee <- apply(dat, 1, 
                function(df){
                  df <- matrix(df, nrow = 1)
                  ts <- df[, 1]
                  t1R <- df[, 2:ncol(df)]
                  
                  risk <- GetRiskSet(ts, Ts, event)
                  
                  t2R <- sum(temp2[risk,])
                  t3R <- colSums(matrix(temp3[risk,], ncol = p))
                  
                  return((t1R - t3R / t2R))
                  
                }) %>% 
      matrix(nrow = p)%>%
      rowSums()
    
  }
  
  return(ee)
  
}

#ee1 <- aee1(beta, data_naive, XB, p, alpha)


###################

#s <- matrix(0,nrow=nrow(dat0), ncol = p)
#for (i in 1:nrow(dat0)) {
#  ts <- dat0[i, 1]
#  t1R <- dat0[i, 2:ncol(dat0)]
  
 # risk <- GetRiskSet(ts, Ts, event)
 # nrisk <- length(risk)
  
 # t2R <- sum(muZ[risk])-(1-alpha)*nrisk*mean.muX
 # t3R <- colSums(matrix(psiZ[risk,], ncol = p)) - (1-alpha)*nrisk*mean.psiX
  
 #s[i,] <- (t1R - t3R / t2R)
  
#}

#s <- colSums(s)
#s
################



################################# ecriture de l'equation estimante ###################

#s <- equa_estimat (beta, data_naive, XB, p,Q)
####################### equation naive#########################################################
# equation naive

equa_naive <- function(beta, data_naive, XB, p) {
  m= nrow(XB)
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

#ee_n <-equa_naive (beta, data_naive, XB, p) 
  
######cox model
#fit_true <- coxph(Surv(Time,delta)~.,data = true_data)
#coef_true <- as.vector(fit_true$coefficients)
#var_true <- diag(fit_true$var)


###### naive
#dat <- Generate_data(p2,m,n,rho,beta)
#XB <- dat$XB
#database_A <- dat$database_A
#true_data <- dat$true_data
#donnees <-Matrix_function(n,m,C,XB)
#data_naive <- donnees$data_naive

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

#observed_naive <-coxph_equa_naive(data_naive, XB, p) 
#coef_naive1 <- observed_naive$coef
#conve_naive1 <- observed_naive$converge
#iterations_naive1 <- observed_naive$iterations
###################

#dat <- Generate_data(p2,m,n,rho,beta)
#XB <- dat$XB
#database_A <- dat$database_A
#true_data <- dat$true_data
#donnees <-Matrix_function(n,m,C,XB)
#data_naive <- donnees$data_naive
#Q <- donnees$Q
#matrice_trans <- donnees$Lvec

coxph_equa_estimat2 <- function(data_naive, XB, p, C,maxiter = 20){
  f2 <- function(x){
    equa_estimat (beta=x, data_naive=data_naive, XB = XB, p = p, C=C )
  }
  fit_manual <- multiroot(f2,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  f.root <- fit_manual$f.root
  iterations2 <-  fit_manual$iter
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  return(list(coef = beta, converge = converge,iterations=iterations2, f.root = f.root   ))
}

#observed2 <- coxph_equa_estimat2(data_naive, XB, p, Q)
#coef_estimat <- observed2$coef
#conve2 <- observed2$converge
#iterations2 <- observed2$iterations
####################################################

coxph_aee1 <- function(data_naive, XB, p, alpha, Q=1, maxiter = 20){
  f <- function(x){
    aee1(beta = x, data_naive=data_naive, XB = XB, p = p, alpha = alpha, Q = 1)
  }
  fit_manual <- multiroot(f,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  f.root <- fit_manual$f.root
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list(coef = beta, converge = converge, f.root = f.root))
}

  
#observed1 <- coxph_aee1 (data_naive, XB, p, alpha,  maxiter = 20)
#coef_estimat1 <- observed1$coef
#conve1 <- observed1$converge
#iterations1 <- observed1$iterations
################################################

library(simsalapar)
library(doParallel)

##################################################

doOne2d_equat_estimat<- function(n,C,alpha){
  # Fixed parameters
  m = 2*n
  rho = 0.75
  p = 2
  beta = matrix(c(0.5,-0.5), nrow = 2)
  
  # Generate data
  data <-Matrix_function(n,m,C,rho,beta)
  data_true <- data$data_true
  data_naive <- data$data_naive
  
  XB = data$XB
  Q <- data$Q
  
  Ts <-as.matrix( data_naive[,1]) 
  event <- data_naive[,2] 
  Z <- as.matrix(data_naive[,3:(p+2)])
  X <- as.matrix(XB[,1:p])
  
  ezbeta <- exp(Z%*% beta)
  eXbeta<-exp(X%*% beta)
  
  zezbeta <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
  XeXbeta<- X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  
  muZ <- exp(Z%*% beta)
  psiZ <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
  
  
  mean.muX <- mean(exp(X%*%beta))
  mean.psiX <- colMeans((X*matrix(rep(exp(X %*% beta),p),ncol =p)))
  
  
  # Theoretical estimating equation for true, naive and validation data
  fit_true <- coxph(Surv(Time,delta)~.,data = data_true)
  coef_true <- as.vector(fit_true$coefficients)
  var_true <- diag(fit_true$var)
  
  fit_naive <- coxph(Surv(Time,delta)~.,data = data_naive)
  coef_naive <- as.vector(fit_naive$coefficients)
  var_naive <- diag(fit_naive$var)
  
  
  #---------------------- Adjusted estimating equation with true alpha ------------------
  
  fit_estimat <- coxph_equa_estimat2(data_naive, XB, p, C,maxiter = 20)
  coef_estimat <- fit_estimat$coef
  converge_estimat <- fit_estimat$converge
  
  fit_naive2 <- coxph_equa_naive(data_naive, XB, p,maxiter = 20)
  coef_naive2 <- fit_naive2$coef
  converge_naive2 <- fit_naive2$converge
  
  fit_aee1 <- coxph_aee1(data_naive, XB, p, alpha,Q=1, maxiter = 20)
  coef_aee1 <- fit_aee1$coef
  converge_aee1<- fit_aee1$converge
  
  return(list ( coef_true=coef_true, coef_naive2=coef_naive2, converge_naive2=converge_naive2,  coef_aee1=  coef_aee1,
                coef_estimat=coef_estimat,converge_estimat =converge_estimat,  converge_aee1= converge_aee1))
}


################## monte carlo ##############################################################
nsim=1000

beta=c(0.5,-0.5)
C=c(0.8,0.2,0.1)
m=1000
n=500
rho=0.75
alpha=0.8
p=2
alpha_q=0.8



coef_estimat <- matrix(0,nrow = nsim, ncol = p)
converge_estimat <- vector()

test_estimat <- function(nsim, C,alpha){
  for (i in 1:nsim) {
    g <- doOne2d_equat_estimat(n,C,alpha)
    
    coef_estimat[i,] <- g$coef_estimat
    converge_estimat[i] <- g$converge_estimat
  
  }
  L1 <- which(converge_estimat!=0)
  coef_estimat<- coef_estimat[L1,]
  result_estimat <- colMeans(coef_estimat)
  
  return(list(coef_estimat=coef_estimat,L1=L1))
}

result_estimat <- test_estimat(1000,C,alpha)
coef_estimat1 <- colMeans(result_estimat$coef_estimat)
L1 <- result_estimat$L1

#############################################################
Q=1

coef_aee1 <- matrix(0,nrow = nsim, ncol = p)
converge_aee1 <- vector()

test_aee1 <- function(nsim, C,alpha){
  for (i in 1:nsim) {
    g <- doOne2d_equat_estimat(n,C,alpha)
    
    coef_aee1[i,] <- g$coef_aee1
    converge_aee1[i] <- g$converge_aee1
  }
 
  L3 <- which(converge_aee1!=0)
  coef_aee1<- coef_aee1[L3,]
  result_aee1 <- colMeans(coef_aee1)
  
  return(coef_aee1=coef_aee1)
}

result_aee1 <- test_aee1(1000,C,alpha)
coef_aee11<-colMeans(result_aee1)


#######################################################

coef_naive<- matrix(0,nrow = nsim, ncol = p)
converge_naive <- vector()

test_naive <- function(nsim, C,alpha){
  for (i in 1:nsim) {
    g <- doOne2d_equat_estimat(n,C,alpha)
    
    coef_naive[i,] <- g$coef_naive2
    converge_naive[i] <- g$converge_naive2
  }
  
  L2 <- which(converge_naive!=0)
  coef_naive<- coef_naive[L2,]
  result_naive <- colMeans(coef_naive)
  
  return(coef_naive=coef_naive)
}

result_naive <- test_naive(1000,C,alpha)
coef_naive1 <-colMeans(result_naive)




#############

h<-doOne2d_equat_estimat(n,C,alpha)
coef_true<-h$coef_true
#################################
H<-cbind(result_naive,result_aee1 ,result_estimat$coef_estimat)
boxplot(result_naive, main=" beta naive")
boxplot(result_estimat, main=" beta estimate ")
boxplot(result_aee1, main=" beta huan ")

###############################################################################


R <- 10 ## nombre de generation de donn?es

m <- 10 # taille de la base B
n <- 5# taille de la base A
p2 <- 0.8# proba de la loi binomiale
rho <- 0.75 # proportion of non censoring data
beta <- c ( 0.5,-0.5)  #parameter
p<-2 ## nombre de covariables
alpha<-0.8 #Probability of True links
C<-c(0.8,0.2,0.1)## posterior probability values

dat <- Generate_data(p2,m,n,rho,beta)
XB <- dat$XB
database_A <- dat$database_A
true_data <- dat$true_data
donnees <-Matrix_function(n,m,C,XB)
data_naive <- donnees$data_naive
Q <- donnees$Q
matrice_trans <- donnees$Lvec





BETA_naive <-matrix(0,nrow = R,ncol = p)
BETA_estimat <-matrix(0,nrow = R,ncol = p)
BETA_huang <- matrix(0,nrow = R,ncol = p)
conve2 <- vector()
conve_naive1<- vector()
conve1 <- vector()


evalution <- function(C,m,n,R){
   
  for (i in 1:R) {
       
    dat <- Generate_data(p2,m,n,rho,beta)
    XB <- dat$XB
    database_A <- dat$database_A
    true_data <- dat$true_data
    donnees <-Matrix_function(n,m,C,XB)
    data_naive <- donnees$data_naive
    data_naive <- donnees$data_naive
    Q <- donnees$Q
    matrice_trans <- donnees$Lvec
    
    observed_naive <-coxph_equa_naive(data_naive, XB, p) 
    BETA_naive[i,] <- observed_naive$coef
    conve_naive1[i] <- observed_naive$converge
    
    observed2 <- coxph_equa_estimat2(data_naive, XB, p, Q)
    BETA_estimat[i,]<- observed2$coef
    conve2[i] <- observed2$converge
    
    observed1 <- coxph_aee1 (data_naive, XB, p, alpha )
    BETA_huang [i,]<- observed1$coef
    conve1[i] <- observed2$converge
    
  }
  return(list(BETA_naive=BETA_naive, BETA_estimat=BETA_estimat,  BETA_huang= BETA_huang,
              conver_naive=conve_naive1,  conver_estimat= conve2,conve_huang=conve1))
}




result <- evalution(C,m,n,R)

X1<- result$BETA_naive
X2<- result$BETA_estimat
X3 <- result$BETA_huang
conver_naive <- result$ conver_naive
conver_estimat <- result$ conver_estimat
conver_huang <- result$conve_huang
base_result <- cbind(X1,X2,X3,conver_naive ,conver_estimat,conver_huang)
colnames(base_result)[1:2] <- c("x1_naive","x2_naive")
colnames(base_result)[3:4] <- c("x1_estimat", "x2_estimat")
colnames(base_result)[5:6] <- c("x1_huang", "x2_huang")



base_result <- data.frame(base_result)
w <- which(base_result[,8]!=0)
w1 <- which(base_result[,9]!=0)
base_result <- base_result[w,]
base_result_huang <- base_result[w1,]

BOX_naive <- boxplot(base_result[,1:2])
BOX_estimat <-  boxplot(base_result[,3:4])
BOX_huang <-  boxplot(base_result_huang[,5:6])

BETA_naive <- colMeans(base_result[,1:2])
BETA_estimat  <-  colMeans(base_result[,3:4])
BETA_huang  <-  colMeans(base_result_huang[,5:6])


fit_true <- coxph(Surv(Time,delta)~.,data = true_data)
BETA_true <- as.vector(fit_true$coefficients)

################################################################




