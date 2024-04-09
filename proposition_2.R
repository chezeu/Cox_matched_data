
library(survival)
library(plyr)
library(dplyr)
library(rootSolve) 
library(nleqslv)

m <- 100 # taille de la base B
n <- 50# taille de la base A
p2 <- 0.8# proba de la loi binomiale
rho <- 0.75 # proportion of non censoring data
beta <- c ( 0.5,-0.5)  #parameter
p<-2 ## nombre de covariables
alpha<-0.8 #Probability of True links
C<-c(0.8,0.2,0.1)## posterior probability values


## generer les bases de données######################################################
Generate_data <- function(p2,m,n,rho,beta){
  # set.seed(42)
  ##covarites data
  X1 <- rnorm(m,0,1)
  X2 <- rbinom(m,size = 1, prob = p2)
  X <- as.matrix(cbind(X1, X2))
  
  U <- runif(m,0,1) # Time to event variable
  Tt = -log(U)/exp((X%*%beta))  #lambda=1
  # Constant right censored
  c = quantile(Tt[,1], probs = rho)
  
  Time <- pmin(Tt,c)
  delta <- as.numeric(Tt<=c)
  surv_data <- data.frame(Time,delta,X)
  
  true_data <- surv_data[1:n,]
  
  database_A <- true_data[,1:2]  #base de donnee de A
  XB <- surv_data[,3:dim(surv_data)[2]] #base de donnee de covariables
  return(list(surv_data=surv_data,true_data = true_data, database_A  = database_A ,XB  =XB ))
}

#########  multinomiale

Matrix_function<-function(n,m,C){
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
      pro[i] = alpha #Probability of True links
      pro[i+1]= C[2]
      Q[i,]<-pro
      Lvec[i,] <- rmultinom(1, size = 1, pro = pro)} 
    if(i==n){  pro = rep(0,m)
    pro[i] = alpha #Probability of True links
    pro[i-1]=C[2]
    Q[i,]<-pro
    Lvec[i,] <- rmultinom(1, size = 1, pro= pro) }  
  }
  Lvec
  X_naive <-Lvec %*%as.matrix( XB)
  Z<- X_naive
  colnames(Z)<- c("Z1","Z2")
  # Naive data
  data_naive <- cbind( database_A, Z)
  return(list(Z=Z, data_naive=data_naive,Q=Q, Lvec= Lvec))
}

####################### equation naive#########################################################

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
    num_naive<-colSums( zezbeta[risq,])
    denum_naive<- sum(ezbeta[risq])    
    s<- s+ event[i]* ( Z[i,]- num_naive/denum_naive)
  }
  s
  return(H_naive=s) ##naive estimator 
}

################################# ecriture de l'equation estimante ###################

equa_estimat2 <- function(beta, data_naive, XB, p,Q) {
  m= nrow(XB)
  n = nrow(data_naive)
  
  Ts <-as.matrix( data_naive[,1]) 
  event <- data_naive[,2] 
  Z <- as.matrix(data_naive[,3:(p+2)])
  X <- as.matrix(XB[,1:p])
  
  ezbeta <- exp(Z%*% beta)
  eXbeta<-exp(X%*% beta)
  
  zezbeta <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
  XeXbeta<- X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  donnees <-Matrix_function(n,m,C)
  Q <- donnees$Q
  
  
  #### autre ecriture de tilde_f
  produit3 <- vector()
  for (j in 1:n) {
    qj <- Q[j,]
      qjx <- (qj*eXbeta)[-j]
      L1 <- which(qj!=0)
    eX <-  qjx[L1]
    produit3[j] <- 1/(qj[j]*prod(eX))
  }  
  produit3
  tilde_f2 <- ezbeta* produit3
  
  
  ### une autre ecriture de tilde_g
  
  produit5<-vector()
  som5<-matrix(0,nrow = n,ncol = p)
  for (k in 1:n) {
    qk <- Q[k,]
    qkx <- (qk*eXbeta)[-k]
    L2 <- which(qk !=0)
    eX <-  qkx[L2]
   Xk <- X[k,]
   xk2 <- X[L2,]
   a <- which(xk2[,1]==Xk[1] & xk2[,2]==Xk[2]) 
   xk2 <- xk2[-a,]
  produit5[k]<- 1/ (qk[k]*prod(eX))
  if(is.null(nrow(xk2)) ) { som5[k,]<- xk2 * eXbeta[k]       
  }else if(1<nrow(xk2)){ 
    som5[k,]<-  colSums(xk2) * eXbeta[k]}
  }  
  produit5  
  som5
  tilde_g2 <- zezbeta*produit5 - som5
  
  ######## determinons tilde_x
  
  tilde_x <- matrix(0,nrow=n, ncol = p)
  for(i in 1:n){ 
    qi <- Q[i,]
    qix <- qi*X
    tilde_x[i,] <- (1/qi[i])*Z[i,]  - colSums(qix)*(1/qi[i])
  }
  tilde_x
  
  
  ## reecriture de l'equation generale
  
  s2<-0
  for(i in 1:n){ 
    Yi <- as.numeric( (Ts == Ts[i] ) | (Ts > Ts[i]))# at risk 
    risq <- which(Yi==1)
    
    num <-colSums( tilde_g2[risq,])
    denum <- sum(tilde_f2[risq])    
    s2<- s2+ event[i]*( tilde_x[i,]- num/denum )
  }
  s2
  return(H=s2) ##naive estimator 
}
#H2<-equa_estimat(beta, data_naive, XB, p,Q)

################### naive equation solving ##################################################

dat <- Generate_data(p2,m,n,rho,beta)
XB <- dat$XB
database_A <- dat$database_A
true_data <- dat$true_data
donnees <-Matrix_function(n,m,C)
data_naive <- donnees$data_naive

coxph_equa_naive <- function(data_naive, XB, p,maxiter = 50){
  f <- function(x){
    equa_naive (beta=x, data_naive=data_naive, XB = XB, p = p )
  }
  
  fit_manual <- multiroot(f,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  iterations <- fit_manual$iter
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list(coef = beta, converge = converge,iterations=iterations))
}

observed_naive <-coxph_equa_naive(data_naive, XB, p) 
coef_naive1 <- observed_naive$coef
conve_naive1 <- observed_naive$converge
iterations_naive1 <- observed_naive$iterations

############# writting the estimate equation solving ###########################################

dat <- Generate_data(p2,m,n,rho,beta)
XB <- dat$XB
database_A <- dat$database_A
true_data <- dat$true_data
donnees <-Matrix_function(n,m,C)
data_naive <- donnees$data_naive
Q <- donnees$Q
matrice_trans <- donnees$Lvec

coxph_equa_estimat2 <- function(data_naive, XB, p, Q,maxiter = 20){
  f2 <- function(x){
    equa_estimat2 (beta=x, data_naive=data_naive, XB = XB, p = p,Q=Q )
  }
  init2 <- rep(0,p)
  g2<- nleqslv( init2 ,f2, method="Newton")
  beta2 <-g2$x
  converge <- as.numeric((g2$iter < maxiter)& (g2$termcd==1) & !is.nan(g2$fvec))
  iterations2 <- g2$iter
  return(list(coef = beta2, converge = converge,iterations=iterations2  ))
}

observed2 <- coxph_equa_estimat2(data_naive, XB, p, Q)
coef_estimat <- observed2$coef
conve2 <- observed2$converge
iterations2 <- observed2$iterations


################## monte carlo ##############################################################
R <- 1000 ## nombre de generation de données

evalution <- function(C,m,n,R){
  
  BETA_naive <-matrix(0,nrow = R,ncol = p)
  BETA_estimat <-matrix(0,nrow = R,ncol = p)
  conve2 <- matrix(0,nrow = R, ncol = p)
  conve_naive1<- vector()
  
  for (i in 1:R) {
    
    p2 <- 0.8# proba de la loi binomiale
    rho <- 0.75 # proportion of non censoring data
    beta <- c ( 0.5,-0.5)  #parameter
    p<-2 ## nombre de covariables
    
    dat <- Generate_data(p2,m,n,rho,beta)
    XB <- dat$XB
    database_A <- dat$database_A
    true_data <- dat$true_data
    donnees <-Matrix_function(n,m,C)
    data_naive <- donnees$data_naive
    data_naive <- donnees$data_naive
    Q <- donnees$Q
    matrice_trans <- donnees$Lvec
    
    observed_naive <-coxph_equa_naive(data_naive, XB, p) 
    coef_naive1 <- observed_naive$coef
    conve_naive1[i] <- observed_naive$converge
    
    observed2 <- coxph_equa_estimat2(data_naive, XB, p, Q)
    coef_estimat <- observed2$coef
    conve2[i,] <- observed2$converge
    
    BETA_naive[i,]  <- coef_naive1
    BETA_estimat[i,]  <- coef_estimat
    
  }
  return(list(BETA_naive=BETA_naive, BETA_estimat=BETA_estimat, 
              conver_naive=conve_naive1,  conver_estimat= conve2))
}

result <- evalution(C,m,n,R)

X1<- result$BETA_naive
X2<- result$BETA_estimat
conver_naive <- result$ conver_naive
conver_estimat <- result$ conver_estimat
base_result <- cbind(X1,X2,conver_naive ,conver_estimat)
colnames(base_result)[1:2] <- c("x1_naive","x2_naive")
colnames(base_result)[3:4] <- c("x1_estimat", "x2_estimat")
colnames(base_result)[6:7] <- c("conv_x1_estimat", "conv_x2_estimat")


base_result <- data.frame(base_result)
w <- which((base_result[,6] & base_result[,7])!=0)
base_result <- base_result[w,]

BOX_naive <- boxplot(base_result[,1:2])
BOX_estimat <-  boxplot(base_result[,3:4])
BETA_naive <- colMeans(base_result[,1:2])
BETA_estimat  <-  colMeans(base_result[,3:4])

fit_true <- coxph(Surv(Time,delta)~.,data = true_data)
BETA_true <- as.vector(fit_true$coefficients)

################################################################


