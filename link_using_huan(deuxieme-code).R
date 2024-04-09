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



Generate_surv_data <- function(m,n,beta,p){
  # set.seed(42)
  ##covarites data
  X1 <- rnorm(m,0,1)
  X2 <- rbinom(m,size = 1, prob = 0.7)
  X <- as.matrix(cbind(X1, X2))
  
  U <- runif(m,0,1) # Time to event variable
  Tt = -log(U)/exp((X%*%beta))  #lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = rho), (rho=0.75)
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

#################### donnees binaires 

# Make error for one binary vector x
makeError <- function(x,error){
  #x is a column of B
  #error is a proportion of error
  nE  = round(length(x)*error)
  index = sample(1:length(x), nE)
  x[index] = 1-x[index]
  return(x)
}


generate_link_data <- function(m, n, K, prevalence, error, min_prev = 0.01){
  #prevalence = rep(c(0.05,0.1,0.2,0.3,0.4), K/5)
  surv_data <- Generate_surv_data (m,n,beta,p)
  
  # First database B
  datB = matrix(0, nrow = m, ncol = K+1)
  datB[,K+1] = 1:m #id
  
  conditionB = TRUE
  while (conditionB){
    datB[,1:K] = sapply(prevalence, function(x){rbinom(n=m, size = 1,prob = x)})
    conditionB = (sum(colSums(datB[,1:K]/m) >= min_prev) < K)
  }
  
  datB = data.frame(datB)
  colnames(datB)=c(paste("R", 1:K, sep = ""),"id") 
  
  conditionA = TRUE
  while (conditionA) {
    # Second database A
    idBA <- sample(1:m,n) #ident in A appearing in B
    datA <- datB[idBA,]
    
    # Make error for A
    datA[,1:K]= apply(datA[,1:K], MARGIN = 2, FUN = makeError, error = error)
    
    conditionA = (sum(colSums(datA[,1:K]) >= 1) < K)
  }
  
  databaseB <- matrix(0,m,(p+K+1))
  databaseB[,1:2] <- as.matrix(surv_data[,3:(p+2)])
  databaseB [,3: ncol(databaseB) ]<-as.matrix( datB)
  databaseB <- data.frame(databaseB)
  colnames(databaseB)=c("x1","x2",paste("R", 1:K, sep = ""),"id") 
  
  databaseA <- matrix(0,n,(p+K+1))
  databaseA[,1:2] <- as.matrix(surv_data[idBA,1:2])
  databaseA [,3: ncol(databaseA) ]<-as.matrix( datA)
  databaseA <- data.frame(databaseA)
  colnames(databaseA)=c("Time","delta",paste("R", 1:K, sep = ""),"id") 
  
  
  datA <- databaseA[,-c(1,2)]
  datB <-  databaseB[,-c(1,2)]
  true_data <- surv_data[idBA,]
  
  return(list(databaseA=databaseA, databaseB = databaseB, prev = prevalence,datA=datA, datB= datB,true_data=true_data))
}


##########################################"" Record Linkage#######################

compare_binary <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid( XB.k,XA.k)
    
    gamma.k = as.numeric(temp[,1]==temp[,2])
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  
  
  return(comp_mat)
}


#T_id <- expand.grid(datB[,K+1],datA[,K+1] )


EM_binary <- function( comp_mat,datA, datB, K, e =0.02, tol = 1e-5, maxits = 500){
 # Starting point
  nA =  nrow(datA)
  N <- nrow(comp_mat)
  p = nA/N
  
  prev = (colMeans(datB[,1:K]))
         
  u = ((1-e)*(1-prev) + e*prev)*(1-prev) + prev*((1-e)*prev+e*(1-prev))
  m = rep(1 - e,K)
  
  # initializations
  comp_mat  <- comp_mat[,1:K]
  
  g = rep(0,N) # probability of being in Match  for each pair l
  it = 0
  converge = FALSE
  
  
  while ((!converge) & (it < maxits)){ 
    
    p.old = p
    m.old = m
    u.old = u
    ### E
    # Compute expectation
    
    m_mat = matrix(rep(m,N), nrow =N, byrow = TRUE)
    u_mat = matrix(rep(u,N), nrow =N, byrow = TRUE)
    
    
    probM = rowProds(m_mat^(comp_mat)*(1-m_mat)^(1-comp_mat))
    probU = rowProds(u_mat^(comp_mat)*(1-u_mat)^(1-comp_mat))
    
    
    g = p*probM/(p*probM+(1-p)*probU)
    
    ### Maximization
    g_mat = matrix(rep(g,K),ncol = K)
    
    p = sum(g)/N
    m = colSums(g_mat*comp_mat)/sum(g)
    u = colSums((1-g_mat)*comp_mat)/sum(1-g)
    
    if (length(which(m > 0.99999)) > 0) {
      m[which(m > 0.99999)] <- 0.99999
    }
    if (length(which(m < 1e-05)) > 0) {
      m[which(m < 1e-05)] <- 1e-05
    }
    if (length(which(u > 0.99999)) > 0) {
      u[which(u > 0.99999)] <- 0.99999
    }
    if (length(which(u < 1e-05)) > 0) {
      u[which(u < 1e-05)] <- 1e-05
    }
    
    
    it = it + 1
    
    converge <- (abs(p.old - p)/p.old < tol) && 
      all(abs(m.old - m)/m.old < tol) && 
      all(abs(u.old - u)/u.old < tol)
    
    if (it == maxits) {
      cat("WARNING! NOT CONVERGENT!", "\n")
      converge = FALSE
    }
  }
  
  comp_mat.1 <- as.data.frame(comp_mat)
  comp_mat.1 $g <- g
  
  
  nA <- nrow(datA)
  nB <-  nrow(datB)
  matrice_norm <- t( matrix(g,ncol=nA,nrow=nB))
  #chaque ligne correspond aux probabilites pour chaque elt de la base B
  
 # h <- which(  matrice_norm<=10e-9) 
 # matrice_norm[h] <- 0
  
  vect_sum <- apply(matrice_norm,1, sum)
  vect_sum <-matrix( rep(vect_sum,nA), ncol = nB,nrow = nA)
  
  ##matrice de prob normalisee
  
  Q <- matrice_norm/vect_sum
  
  return(list(g=g, p=p, m=m, u = u, comp_mat=comp_mat,comp_mat.1=comp_mat.1,Q=Q, it = it, converge = converge))
}

generate_naive_data <- function(comp_mat.1,databaseB, databaseA){
  
  surv_data <- Generate_surv_data (m,n,beta,p)
    
  X <- surv_data[,3:ncol(surv_data)]
  XB <- as.matrix(X, ncol = p) #Only covariates
  
  n_idA <-databaseA$id
  #data_true <- surv_data[n_idA,]
  
 g <- comp_mat.1$g
   
  i0 <- 1
  nA <- nrow(databaseA)
  nB <- nrow(databaseB)
  k <- nB
  j <- 1
  N<- nA*nB
  a <- vector()
  
  while(k<=N){
    b <- vector()
    for (i in i0:(j*nB)) { b[i] <- g[i] }
    prob_max <- max(b[i0:(j*nB)])
    position <- which( b[i0:(j*nB)]== prob_max)
    if(i0==1){
      if(length(position)==1){ a[j] <- position }
      if(length(position)!=1) { a[j] <- sample(position,1) }
    }
    if(i0!=1){
      if(length(position)==1){ a[j] <- position + (j-1)*nB}
      if(length(position)!=1) { a[j] <- sample(position,1)+(j-1)*nB }
    }
    i0 <- i0+nB
    j <- j+1
    k <- k+nB
    
  }
  
  
  
  idA_B <- expand.grid (datB[,K+1],datA[,K+1])
  comp_mat.1 <- data.frame(comp_mat.1)
  comp_mat.1$idA <- idA_B$Var2
  comp_mat.1$idB <- idA_B$Var1
  Naive_max <-comp_mat.1[a,]
  lB <- Naive_max$idB
  lA <- Naive_max$idA
  
  data_naive <- databaseA[,1:2]
  covar_naive_B <- databaseB[lB,]
  
  data_naive$x1 <- covar_naive_B$x1
  data_naive$x2 <-  covar_naive_B$x2
  
  
  f <- which( Naive_max$idA != Naive_max$idB)
  r <-  (ncol(Naive_max)-2)
  FP <- Naive_max[f,(r:(r+2))]
  
  return(list(data_naive=data_naive, Naive_max=Naive_max,   XB=  XB, f=f, FP=FP))
}


######################
Matrix_function <-function(n,m,C,beta,databaseB, databaseA, estimate, comp_mat){
  
  surv_data <- Generate_surv_data(m,n,beta,p)
  X <- surv_data[,3:ncol(surv_data)]
  XB <- as.matrix(X, ncol = p) #Only covariates
  
  ###### True-link data: The first n individual from B
  data_true <- surv_data[1:n,]
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  
  
  if(estimate==TRUE){ 
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
        Q[i,]<-pro
        Lvec[i,] <- rmultinom(1, size = 1, pro = pro)} 
      if(i==n){  pro = rep(0,m)
      pro[i] = C[1] #Probability of True links
      pro[i-1]=C[2]
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
    
  }else{
    
    C=c(0,0,0)
    EM <- EM_binary (comp_mat, datA, datB, K, e = 0.02, tol = 1e-5, maxits = 500)
    Q <- EM$Q
   
    
  }
  return(list( data_A= data_A,XB=XB,Q=Q))
}


#func <- Matrix_function (n,m,C,beta,databaseB, databaseA, estimate=FALSE,comp_mat)
#Q <- func$Q
#XB <- func$XB
######################################################### cox model ###########


#1. Finding the risk set R(t) given some time t
GetRiskSet <- function(time_of_interest, time_vector, event_vector) {
  
  return(which(((time_vector == time_of_interest & event_vector == 1) | (time_vector > time_of_interest))))
  
}
#################

equa_estimat <- function(beta, data_naive, XB, p,Q) {
  
  Ts <-as.matrix( data_naive[,1]) 
  event <- data_naive[,2] 
  Z <- as.matrix(data_naive[,3:(p+2)])
  X <- as.matrix(XB[,1:p])
  
  ezbeta <- exp(Z%*% beta)
  eXbeta<-exp(X%*% beta)
  
  zezbeta <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
  XeXbeta<- X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
 
  
  #######les sommes pour tilde_x
  
  som0 <- matrix(0,nrow=n, ncol = p)
  diagonal <- vector()
  for(i in 1:n){ 
    qi <- Q[i,]
    # L <- which(qi<=10^-3)
    # qi[L] <- 0
    qii <- max(qi)
    ii <- which(qi==qii)
    qix <- (qi*X)[-ii,]
    
    som0[i,] <- colSums(qix)
    diagonal[i] <- qii
  }
  som0
  diagonal
  tilde_x <- (Z/ diagonal)-som0/diagonal  
  
  
  #### les sommes q*exp(beta x) pour tilfe_f
  som1<- vector()
  for (j in 1:n) {
    qj <- Q[j,]
    #L <- which(qj<=10^-3)
    #qj[L] <- 0
    qjj <- max(qj)
    l <- which(qj==qjj)
    qjx <- (qj*eXbeta)[-l]
    som1[j] <- sum( qjx)
  }  
  som1 
  
  
  ###  les sommes pour tilde_g
  
  som2<-matrix(0, nrow = n, ncol = p)
  
  for (k in 1:n) {
    qk <- Q[k,]
    # L <- which(qk<=10^-3)
    # qk[L] <- 0
    qkk <- max(qk)
    l <- which(qk==qkk)
    qkx <- (qk*XeXbeta)[-l,]
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

#esti_H <- equa_estimat (beta, data_naive, XB,p, Q) 



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

#naive_H <- equa_naive (beta, data_naive, XB, p) 
##################################################################
################### resoudre les equations 

coxph_equa_naive <- function(data_naive, XB, p,maxiter = 100){
  f <- function(x){
    equa_naive (beta=x, data_naive=data_naive, XB = XB, p = p )
  }
  
  fit_manual <- multiroot(f,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  iterations <- fit_manual$iter
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list(coef = beta, converge = converge,iterations=iterations))
}


coxph_equa_estimat2 <- function(data_naive, XB, p,Q,maxiter = 100){
  f2 <- function(x){
    equa_estimat (beta=x, data_naive=data_naive, XB = XB, p = p, Q=Q )
  }
  fit_manual <- multiroot(f2,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  f.root <- fit_manual$f.root
  iterations2 <-  fit_manual$iter
  #fit_manual <- nleqslv( c(0,0),f2, method="Newton")
  #beta <- fit_manual$x
  #iterations2 <- fit_manual$iter
  #f.root <- fit_manual$fvec
 # set.seed(14)
  # xstart <- matrix(rnorm(20,0,10), ncol = 2)
  #Zero <-  searchZeros(xstart,f2)
  
  
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  return(list(coef = beta, converge = converge,iterations=iterations2, f.root = f.root   ))
}


################################################

library(simsalapar)
library(doParallel)

################################################## estimation des parametres ##

doOne2d_equat_estimat<- function(n,Q,databaseB, databaseA,true_data){
  # Fixed parameters
  #m = 2*n
  #rho = 0.75
  p = 2
  beta = matrix(c(0.5,-0.5), nrow = 2)
  
  
  Naive <- generate_naive_data(comp_mat.1,databaseB, databaseA)
  naive_max <- Naive$Naive_max
  data_naive<- Naive$data_naive
  XB <- as.matrix(Naive$XB)
  
  Ts <-as.matrix( data_naive [,1]) 
  event <-  data_naive[,2]
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
  fit_true <- coxph(Surv(Time,delta)~.,data = true_data)
  coef_true <- as.vector(fit_true$coefficients)
  var_true <- diag(fit_true$var)
  
  fit_naive <- coxph(Surv(Time,delta)~.,data = data_naive)
  coef_naive <- as.vector(fit_naive$coefficients)
  var_naive <- diag(fit_naive$var)
  
  
  #---------------------- Adjusted estimating equation with true alpha ------------------
  
  fit_estimat <- coxph_equa_estimat2(data_naive, XB, p, Q,maxiter = 100)
  coef_estimat <- fit_estimat$coef
  converge_estimat <- fit_estimat$converge
  
  fit_naive2 <- coxph_equa_naive(data_naive, XB, p,maxiter = 100)
  coef_naive2 <- fit_naive2$coef
  converge_naive2 <- fit_naive2$converge
  
  return(list (  coef_true= coef_true, coef_naive= coef_naive , coef_naive2=coef_naive2, converge_naive2=converge_naive2,
                 coef_estimat=coef_estimat,converge_estimat =converge_estimat  ))
}

#parameters <- doOne2d_equat_estimat (n,Q,databaseB, databaseA)

################################################### simulations


m=500
n=200
K=30
p=2

error=0.02# pas de convergence lorsque error=0.1 
beta=c(0.5,-0.5)
prevalence = rep(c(0.05,0.1,0.2,0.3,0.4), K/5)

surv_data <- Generate_surv_data (m,n,beta,p)

gen <- generate_link_data (m, n, K, prevalence, error, min_prev = 0.01)
databaseA <- gen$databaseA
databaseB <- gen$databaseB
prev <- gen$prev
datB <- gen$ datB
datA <- gen$ datA
true_data <- gen$true_data

comp_mat<- compare_binary(datA, datB, K)

EM <- EM_binary ( comp_mat,datA, datB, K, e = 0.02, tol = 1e-5, maxits = 500)
Q <- EM$Q
comp_mat.1 <- EM$comp_mat.1

Naive <- generate_naive_data(comp_mat.1,databaseB, databaseA)
Naive_max <- Naive$Naive_max
data_naive<- Naive$data_naive
XB <- Naive$XB
f <- Naive$f
FP <- Naive$FP

nsim=500

coef_estimat <- matrix(0,nrow = nsim, ncol = p)
converge_estimat <- vector()

test_estimat <- function(nsim, Q){
  for (i in 1:nsim) {
    g <- doOne2d_equat_estimat (n,Q,databaseB, databaseA,true_data)
    
    coef_estimat[i,] <- g$coef_estimat
    converge_estimat[i] <- g$converge_estimat
    
  }
  L1 <- which(converge_estimat!=0)
  coef_estimat<- coef_estimat[L1,]
  result_estimat <- colMeans(coef_estimat)
  
  return(list(coef_estimat=coef_estimat,L1=L1))
}

result_estimat <- test_estimat(500,Q)
coef_estimat1 <- colMeans(result_estimat$coef_estimat)
nombre_cvgence <- length (result_estimat$L1)
#############################################################

nsim=500

beta=c(0.5,-0.5)
m=400
n=200
p=2

coef_naive2<- matrix(0,nrow = nsim, ncol = p)
converge_naive2 <- vector()

test_naive <- function(nsim, Q){
  for (i in 1:nsim) {
    g <- doOne2d_equat_estimat (n,Q,databaseB, databaseA,true_data)
    
    coef_naive2[i,] <- g$coef_naive2
    converge_naive2[i] <- g$converge_naive2
    
  }
  L1 <- which(converge_naive2!=0)
  coef_naive2<- coef_naive2[L1,]
  result_naive <- colMeans(coef_naive2)
  
  return(list(coef_estimat=coef_naive2,L1=L1))
}

result_naive <- test_naive(500,Q)
coef_naive2 <- colMeans(result_naive$coef_estimat )
nombre_cvgence_naive <-length( result_naive$L1)



























