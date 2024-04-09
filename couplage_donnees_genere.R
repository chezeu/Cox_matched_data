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
  
  
  return(list(databaseA=databaseA, databaseB = databaseB, prev = prevalence))
}



m=40
 n=20
 K=10
p=2

 error=10/100
beta=c(0.5,-0.5)
prevalence = rep(c(0.05,0.1,0.2,0.3,0.4), K/5)

surv_data <- Generate_surv_data (m,n,beta,p)

gen <- generate_link_data (m, n, K, prevalence, error, min_prev = 0.01)
databaseA <- gen$databaseA
databaseB <- gen$databaseB
prev <- gen$prev


##########################################"" Record Linkage#######################
library(reclin2)

linkage_data <- function(databaseB, databaseA) {
  
  surv_data <- Generate_surv_data (m,n,beta,p)
  gen <- generate_link_data (m, n, K, prevalence, error, min_prev = 0.01)
  databaseA <- gen$databaseA
  databaseB <- gen$databaseB
  prev <- gen$prev
  
  X <- surv_data[,3:ncol(surv_data)]
  XB <- as.matrix(X, ncol = p) #Only covariates
  
  ###### True-link data: The first n individual from B
  data_true <- surv_data[1:n,]
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  
  pp <- pair_minsim(databaseA, databaseB, 
                    on =  c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","id"), minsim = 0)
  #pp <- print(pp)
  
  # show the pairs of comparison en laissant la variable 
  #p <- compare_pairs(p, on= c("fname_c1","lname_c1","by","bm","bd"), inplace = TRUE, 
  #             comparators = list(fname_c1 = jaro_winkler(), lname_c1  = jaro_winkler() ))
  
  pp <- compare_pairs(pp, on =  c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","id"), inplace = TRUE)
  #p <- print(p)
  
  ## classification
  model <- problink_em ( ~R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+id, data=pp)
  pp <- predict(model, pp, type ="mpost", add = TRUE, binary = TRUE)
  
  matrice_norm <- pp$mpost
  a <- nrow(databaseA)
  b <-  nrow(databaseB)
  matrice_norm <- t( matrix(matrice_norm,ncol=a,nrow=b))
  #chaque ligne correspond aux probabilites pour chaque elt de la base B
  
  vect_sum <- apply(matrice_norm,1, sum)
  vect_sum <-matrix( rep(vect_sum,a), ncol = b,nrow = a)
  
  ##matrice de prob normalisee
  
  Q <- matrice_norm/vect_sum
  
  matrice_norm <- matrix(t(Q), ncol = 1)
  pp$matrice_norm <- matrice_norm      #proba normalisee
  
  # Select pairs with a mpost > 0.5
  pp <- select_threshold(pp, "selected", "mpost", 0.5, inplace = TRUE)
  L <- which(pp$selected==TRUE)
  Naive_link  <- pp[L,]
  
  #pp <- select_threshold(pp, "selected_2", "matrice_prob", 0.5, inplace = TRUE)
  #L2 <- which(pp$selected_2==TRUE)
  #Naive_true <- pp[L2,]
  
  ##creat a linked data
  #linked_data <- link(pp)
  return(list( Q=Q,Naive_link =Naive_link,pp=pp,databaseB=databaseB,databaseA=databaseA  ))
}



Link <- linkage_data (databaseB, databaseA)
databaseB <-Link$databaseB
databaseA <- Link$databaseA
Q <- Link$Q
pp <- Link$pp
Naive_link <- Link$Naive_link


generate_naive_data <- function(pp,databaseB, databaseA){
  
  surv_data <- Generate_surv_data (m,n,beta,p)
  gen <- generate_link_data (m, n, K, prevalence, error, min_prev = 0.01)
  databaseA <- gen$databaseA
  databaseB <- gen$databaseB
  prev <- gen$prev
  
  X <- surv_data[,3:ncol(surv_data)]
  XB <- as.matrix(X, ncol = p) #Only covariates
  data_true <- surv_data[1:n,]
  
  mpost <- pp$mpost
  i0 <- 1
  n <- n
  m <- m
  k <- m
  j <- 1
  K <- n*m
  a <- vector()
  
  while(k<=K){
    b <- vector()
    for (i in i0:(j*m)) { b[i] <- mpost[i] }
    prob_max <- max(b[i0:(j*m)])
    position <- which( b[i0:(j*m)]== prob_max)
    if(i0==1){
      if(length(position==1)){ a[j] <- position }
      if(length(position!=1)) { a[j] <- sample(position,1) }
    }
    if(i0!=1){
      if(length(position)==1){ a[j] <- position + (j-1)*m}
      if(length(position)!=1) { a[j] <- sample(position,1)+(j-1)*m }
    }
    i0 <- i0+m
    j <- j+1
    k <- k+m
    
  }
  Naive_max <-pp[a,]
  lB <- Naive_max$.y
  lA <- Naive_max$.x
  
  data_naive <- databaseA[,1:2]
  covar_naive_B <- databaseB[lB,]
  
  data_naive$x1 <- covar_naive_B$x1
  data_naive$x2 <-  covar_naive_B$x2
  
  
  return(list(data_naive=data_naive, Naive_max=Naive_max))
}

Naive <- generate_naive_data(pp,databaseB, databaseA)
naive_max <- Naive$Naive_max
data_naive<- Naive$data_naive

######################
Matrix_function <-function(n,m,C,beta,databaseB, databaseA, estimate){
  
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
    Link <- linkage_data (databaseB, databaseA)
    Q <- Link$Q
    
  }
  return(list( data_true=data_true, data_A= data_A,XB=XB,Q=Q))
}


func <- Matrix_function (n,m,C,beta,databaseB, databaseA, estimate=FALSE)
Q <- func$Q
XB <- func$XB
data_true <- func$data_true




#1. Finding the risk set R(t) given some time t
GetRiskSet <- function(time_of_interest, time_vector, event_vector) {
  
  return(which(((time_vector == time_of_interest & event_vector == 1) | (time_vector > time_of_interest))))
  
}
#################

equa_estimat <- function(beta, data_naive, XB, p,C=c(0,0,0),estimate=FALSE) {
  
  
  Ts <-as.matrix( data_naive [,3]) 
  event <-  data_naive[,4]
  Z <- as.matrix(data_naive[,1:2])
  X <- as.matrix(XB[,1:p])
  
  ezbeta <- exp(Z%*% beta)
  eXbeta<-exp(X%*% beta)
  
  zezbeta <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
  XeXbeta<- X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  donnees <-Matrix_function(n,m,C,beta ,databaseB, databaseA, estimate=FALSE)
  Q <- donnees$Q
  h <- which(Q<=10e-3) 
  Q[h] <- 0
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

esti_beta <- equa_estimat (beta, data_naive, XB, p,C,estimate=FALSE) 



#2. Proposed estimating equations of huan
aee1 <- function(beta, data_naive, XB, p, alpha,Q=1) {
  N = nrow(XB)
  n = nrow(data_naive)
  
  Ts <-as.matrix( data_naive [,3]) 
  event <-  data_naive[,4]
  Z <- as.matrix(data_naive[,1:2])
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

huan_beta <- aee1(beta, data_naive, XB, p, alpha,Q=1)


equa_naive <- function(beta, data_naive, XB, p) {
  m= nrow(XB)
  n = nrow(data_naive)
  
  Ts <-as.matrix( data_naive[,3]) 
  event <- data_naive[,4] 
  Z <- as.matrix(data_naive[,1:p])
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

naive_beta <- equa_naive (beta, data_naive, XB, p) 

###################
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


coxph_equa_estimat2 <- function(data_naive, XB, p, C,estimate,maxiter = 20){
  f2 <- function(x){
    equa_estimat (beta=x, data_naive=data_naive, XB = XB, p = p, C=C, estimate )
  }
  fit_manual <- multiroot(f2,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  f.root <- fit_manual$f.root
  iterations2 <-  fit_manual$iter
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  return(list(coef = beta, converge = converge,iterations=iterations2, f.root = f.root   ))
}


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

################################################

library(simsalapar)
library(doParallel)

##################################################

doOne2d_equat_estimat<- function(n,C=c(0,0,0),alpha,estimate=FALSE){
  # Fixed parameters
  m = 2*n
  #rho = 0.75
  p = 2
  beta = matrix(c(0.5,-0.5), nrow = 2)
  
  func <- Matrix_function (n,m,C,beta,databaseB, databaseA, estimate=FALSE)
  Q <- func$Q
  XB <- func$XB
  data_true <- func$data_true
  
  tim_data <- generate_naive_data(pp,databaseB, databaseA)
  data_naive <- tim_data$data_naive
  
  Ts <-as.matrix( data_naive [,3]) 
  event <-  data_naive[,4]
  Z <- as.matrix(data_naive[,1:2])
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
  
  fit_estimat <- coxph_equa_estimat2(data_naive, XB, p, C=c(0,0,0),estimate=FALSE,maxiter = 20)
  coef_estimat <- fit_estimat$coef
  converge_estimat <- fit_estimat$converge
  
  fit_naive2 <- coxph_equa_naive(data_naive, XB, p,maxiter = 20)
  coef_naive2 <- fit_naive2$coef
  converge_naive2 <- fit_naive2$converge
  
  fit_aee1 <- coxph_aee1(data_naive, XB, p, alpha,Q=1, maxiter = 20)
  coef_aee1 <- fit_aee1$coef
  converge_aee1<- fit_aee1$converge
  
  return(list ( coef_true=coef_true, coef_naive2=coef_naive2, converge_naive2=converge_naive2,  coef_aee1=  coef_aee1,
                converge_aee1= converge_aee1,    coef_estimat=coef_estimat,converge_estimat =converge_estimat  ))
}

do <- doOne2d_equat_estimat (n,C,alpha,estimate=FALSE)

###################################################

nsim=100

beta=c(0.5,-0.5)
C=c(0.8,0.2,0.1)
m=40 
n=20
rho=0.75
alpha=0.75
p=2
alpha_q=0.75

coef_estimat <- matrix(0,nrow = nsim, ncol = p)
converge_estimat <- vector()

test_estimat <- function(nsim, C,alpha){
  for (i in 1:nsim) {
    g <- doOne2d_equat_estimat (n,C,alpha,estimate=FALSE)
    
    coef_estimat[i,] <- g$coef_estimat
    converge_estimat[i] <- g$converge_estimat
    
  }
  L1 <- which(converge_estimat!=0)
  coef_estimat<- coef_estimat[L1,]
  result_estimat <- colMeans(coef_estimat)
  
  return(list(coef_estimat=coef_estimat,L1=L1))
}

result_estimat <- test_estimat(100,C,alpha)
coef_estimat1 <- colMeans(result_estimat$coef_estimat)
nombre_cvgence <- result_estimat$L1
#############################################################













######################  compare #############################

# Fellei-Sunter using 3 categorical comparison 

compare3 <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid(XA.k, XB.k)
    if (k == (K+1)){
      gamma.k = as.numeric(temp[,1]==temp[,2])
    }else{
      gamma.k = temp[,1]+temp[,2]
    }
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  return(comp_mat)
}

com3 <- compare3(datA, datB, K)

compare4 <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid(XA.k, XB.k)
    if (k == (K+1)){
      gamma.k = as.numeric(temp[,1]==temp[,2])
    }else{
      gamma.k = 2*temp[,1]+temp[,2]
    }
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  return(comp_mat)
}

com4 <- compare4(datA, datB, K)


compare_binary <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid(XA.k, XB.k)
    
    gamma.k = as.numeric(temp[,1]==temp[,2])
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  return(comp_mat)
}


com_binary <- compare_binary(datA, datB, K)

########################### EM ##################################

EM_binary <- function(comp_mat, datA, datB, K, e = 0.01, tol = 1e-5, maxits = 200){
  
  # Starting point
  nB =  nrow(datB)
  N <- nrow(comp_mat)
  p = nB/N
  
  prev = (colMeans(datA[,1:K]))
  
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
  
  
  
  return(list(g=g, p=p, m=m, u = u, it = it, converge = converge))
  
}


EM_binary_obs <- function(comp_mat, datA, datB, K){
  
  M = comp_mat[comp_mat[,K+1]==1,1:K]
  U = comp_mat[comp_mat[,K+1]==0,1:K]
  
  m = colMeans(M)
  u = colMeans(U)
  
  N = nrow(comp_mat)
  p = nrow(M)/N
  
  
  # initializations
  comp_mat  <- comp_mat[,1:K]
  
  m_mat = matrix(rep(m,N), nrow =N, byrow = TRUE)
  u_mat = matrix(rep(u,N), nrow =N, byrow = TRUE)
  
  
  probM = rowProds(m_mat^(comp_mat)*(1-m_mat)^(1-comp_mat))
  probU = rowProds(u_mat^(comp_mat)*(1-u_mat)^(1-comp_mat))
  
  
  g = p*probM/(p*probM+(1-p)*probU)
  
  
  return(list(g=g))
}


check <- function(param, N){
  if (length(which(param < 1/N)) > 0) {
    param[which(param < 1/N)] <- 1/N
  }
  
  if (length(which(param > (1-1/N))) > 0) {
    param[which(param > (1-1/N))] <- 1-1/N
  }
  return(param)
}

check1 <- function(m,N){
  min.par = 1/N
  if (length(which(m > (1-min.par))) > 0) {
    m[which(m > (1-min.par))] <- 1-min.par
  }
  if (length(which(m < min.par)) > 0) {
    m[which(m < min.par)] <- min.par
  }
  
  return(m)
  
}

start_3 <- function(datA, datB,K){
  e = 0.01 #is q in your documents
  nA = nrow(datA)
  nB = nrow(datB)
  p = 1/nA
  
  # This is the empirical prevalence that I use for estimate starting point
  # To compute the true, you replace it with true parameters
  prev =  colMeans(datA[,1:K])
  
  p0M = (1-e)*(1-prev) #(0,0)
  p1M = rep(e,K) #(0,1) + (1,0)
  p2M = 1-p1M-p0M # (1,1)
  
  p0U = (1-prev)*((1-e)*(1-prev)+e*prev) #(0,0)
  p1U =   (1-prev)*((1-e)*prev + e*(1-prev)) + prev*((1-e)*(1-prev) + e*prev) #(0,1) + (1,0)
  p2U = 1- p1U - p0U # (1,1)
  
  return(c(p, p2M, p1M, p0M, p2U, p1U, p0U ))
  
}

E3 <- function(m22, m21, m20, u22,u21, u20, N, c22, c21, c20){
  m22.mat = matrix(rep(m22,N), nrow =N, byrow = TRUE)
  m21.mat = matrix(rep(m21,N), nrow =N, byrow = TRUE)
  m20.mat = matrix(rep(m20,N), nrow =N, byrow = TRUE)
  
  u22.mat = matrix(rep(u22,N), nrow =N, byrow = TRUE)
  u21.mat = matrix(rep(u21,N), nrow =N, byrow = TRUE)
  u20.mat = matrix(rep(u20,N), nrow =N, byrow = TRUE)
  
  
  p2M = rowProds(m22.mat^(c22)*m21.mat^(c21)*m20.mat^(c20))
  p2U = rowProds(u22.mat^(c22)*u21.mat^(c21)*u20.mat^(c20))
  
  return(cbind(p2M, p2U))
}


M3 <- function(g,K2,c22,c21,c20){
  N = length(g)
  g2.mat = matrix(rep(g,K2),ncol = K2)
  
  m22 = colSums(g2.mat*c22)/sum(g)
  m22 = check1(m22,N)
  m21 = colSums(g2.mat*c21)/sum(g)
  m21 = check1(m21,N)
  m20 = 1-m22-m21 #check 0,1
  
  if (length(which(m20 < 1/N)) > 0){
    m22[which(m20 < 1/N)] = m22[which(m20 < 1/N)] - 1/(2*N)
    m21[which(m20 < 1/N)] = m21[which(m20 < 1/N)] - 1/(2*N)
    m20[which(m20 < 1/N)] = m20[which(m20 < 1/N)] + 1/N 
  }
  
  u22 = colSums((1-g2.mat)*c22)/sum(1-g)
  u22 = check1(u22,N)
  u21 = colSums((1-g2.mat)*c21)/sum(1-g)
  u21 = check1(u21,N)
  u20 = 1-u22-u21 #check 0,1
  
  if (length(which(u20 < 1/N)) > 0){
    u22[which(u20 < 1/N)] = u22[which(u20 < 1/N)] - 1/(2*N)
    u21[which(u20 < 1/N)] = u21[which(u20 < 1/N)] - 1/(2*N)
    u20[which(u20 < 1/N)] = u20[which(u20 < 1/N)] + 1/N 
  }
  return(list(m22 = m22,m21 = m21,m20=m20,u22=u22,u21=u21,u20=u20))
}

EM3 <- function(comp_mat, datA, datB, K, tol = 1e-6, maxits = 500){
  
  # Starting point
  #start <- find_start(comp_mat, datA, datB, K)
  start <- start_3(datA, datB, K)
  
  p = start[1]
  
  m22 = start[2:(K+1)]
  m21 = start[(K+2):(2*K+1)]
  m20 = start[(2*K+2):(3*K+1)]
  
  
  u22 = start[(3*K+2):(4*K+1)]
  u21 = start[(4*K+2):(5*K+1)]
  u20 = start[(5*K+2):(6*K+1)]
  
  
  # initializations
  comp_mat  <- comp_mat[,1:K]
  N <- nrow(comp_mat)
  c20 <- as.numeric(comp_mat==0)
  c21 <- as.numeric(comp_mat==1)
  c22 <- as.numeric(comp_mat==2)
  
  loglikelihood <- function(p,N,m22, m21, m20, u22,u21, u20, c22, c21, c20){
    
    p2 = E3(m22, m21, m20, u22,u21, u20, N, c22, c21, c20)
    pM = p2[,1]
    pU = p2[,2]
    
    ll = sum(log(p*pM+(1-p)*pU)) #likelihood not complete likelihood
    return(ll)
  }
  
  iter <- 0
  converge = FALSE
  diff <- tol + 1
  old.ll = loglikelihood(p,N,m22, m21, m20, u22,u21, u20, c22, c21, c20)
  ll = old.ll
  
  
  while (!converge && iter < maxits){ 
    
    p.old = p
    
    m22.old = m22 #length K2
    m21.old = m21 #length K2
    m20.old = m20
    
    u22.old = u22 #length K2
    u21.old = u21 #length K2
    u20.old = u20
    
    m.old = c(m20.old,m21.old,m22.old)
    u.old = c(u20.old,u21.old,u22.old)
    
    ############################ E step
    p2 = E3(m22, m21, m20, u22,u21, u20, N, c22, c21, c20)
    
    pM = p2[,1]
    pU = p2[,2]
    
    
    # Expectations
    g = p*pM/(p*pM+(1-p)*pU)
    
    ############################### M step
    p = sum(g)/N
    
    max2 = M3(g,K,c22,c21,c20)
    m22 = max2$m22
    m21 = max2$m21
    m20 = max2$m20
    u22 = max2$u22
    u21 = max2$u21
    u20 = max2$u20
    
    m = c(m20,m21,m22)
    u = c(u20,u21,u22)
    
    ### Stopping 
    converge <- (abs(p.old - p)/p.old < tol) && 
      all(abs(m.old - m)/m.old < tol) &&
      all(abs(u.old - u)/u.old < tol)
    
    new.ll = loglikelihood(p,N,m22, m21, m20, u22,u21, u20, c22, c21, c20)
    
    diff <- new.ll - old.ll
    
    old.ll <- new.ll
    ll <- c(ll, old.ll)
    
    iter = iter + 1
    
    
    if (iter == maxits) {
      cat("WARNING! NOT CONVERGENT!", "\n")
      converge = FALSE
    }
    
  }
  
  
  return(list(g=g, loglikelihood = ll, iter = iter, converge = converge))
}



EM4 <- function(comp_mat, datA, datB, K, tol = 1e-6, maxits = 300){
  
  # Starting point
  #start <- find_start(comp_mat, datA, datB, K)
  start <- start_4(datA, datB, K)
  
  p = start[1]
  
  m23 = start[2:(K+1)]
  m22 = start[(K+2):(2*K+1)]
  m21 = start[(2*K+2):(3*K+1)]
  m20 = start[(3*K+2):(4*K+1)]
  
  u23 = start[(4*K+2):(5*K+1)]
  u22 = start[(5*K+2):(6*K+1)]
  u21 = start[(6*K+2):(7*K+1)]
  u20 = start[(7*K+2):(8*K+1)]
  
  # initializations
  comp_mat  <- comp_mat[,1:K]
  N <- nrow(comp_mat)
  c20 <- as.numeric(comp_mat==0)
  c21 <- as.numeric(comp_mat==1)
  c22 <- as.numeric(comp_mat==2)
  c23 <- as.numeric(comp_mat==3)
  
  loglikelihood <- function(p,N,m23,m22, m21, m20,u23, u22,u21, u20, c23, c22, c21, c20){
    
    p2 = E4(m23, m22, m21, m20, u23, u22,u21, u20, N, c23, c22, c21, c20)
    pM = p2[,1]
    pU = p2[,2]
    
    ll = sum(log(p*pM+(1-p)*pU)) #likelihood not complete likelihood
    return(ll)
  }
  
  iter <- 0
  converge = FALSE
  diff <- tol + 1
  old.ll = loglikelihood(p,N,m23,m22, m21, m20, u23, u22,u21, u20, c23, c22, c21, c20)
  ll = old.ll
  
  
  while (!converge && iter < maxits){ 
    
    p.old = p
    
    m23.old = m23
    m22.old = m22 #length K2
    m21.old = m21 #length K2
    m20.old = m20
    
    u23.old = u23
    u22.old = u22 #length K2
    u21.old = u21 #length K2
    u20.old = u20
    
    m.old = c(m20.old,m21.old,m22.old, m23.old )
    u.old = c(u20.old,u21.old,u22.old, u23.old)
    ############################ E step
    p2 = E4(m23,m22, m21, m20, u23, u22,u21, u20, N, c23, c22, c21, c20)
    
    pM = p2[,1]
    pU = p2[,2]
    
    
    # Expectations
    g = p*pM/(p*pM+(1-p)*pU)
    
    ############################### M step
    p = sum(g)/N
    
    max2 = M4(g,K,c23,c22,c21,c20)
    m23 = max2$m23
    m22 = max2$m22
    m21 = max2$m21
    m20 = max2$m20
    
    u23 = max2$u23
    u22 = max2$u22
    u21 = max2$u21
    u20 = max2$u20
    
    
    m = c(m20,m21,m22, m23)
    u = c(u20,u21,u22, u23)
    ### Stopping 
    converge <- (abs(p.old - p)/p.old < tol) && 
      all(abs(m.old - m)/m.old < tol) &&
      all(abs(u.old - u)/u.old < tol)
    
    new.ll = loglikelihood(p,N,m23,m22, m21, m20, u23, u22,u21, u20, c23, c22, c21, c20)
    
    diff <- new.ll - old.ll
    
    old.ll <- new.ll
    ll <- c(ll, old.ll)
    
    iter = iter + 1
    
    
    if (iter == maxits) {
      cat("WARNING! NOT CONVERGENT!", "\n")
      converge = FALSE
    }
  }
  
  
  return(list(g=g, loglikelihood = ll, iter = iter, converge=converge))
}

############################    method


FS <- function(datA, datB, K,tol=1e-6, maxits = 500){
  nA = nrow(datA)
  nB =  nrow(datB)
  
  comp_mat <- compare_binary(datA, datB, K)
  
  fit <- EM_binary(comp_mat, datA, datB, K,tol = tol, maxits = maxits)
  g = fit$g
  converge = fit$converge
  
  
  # Gmat = matrix(g, nrow = nB, byrow = TRUE)
  # Gmat[Gmat<0.5] = 0
  # opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  # predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # temp = (predict_cate[,1]-1)*nA+ predict_cate[,2]
  # # Percentage of correct link
  # nPredict = sum((g[temp]>=0.5))
  # nTruePredict = sum((predict_cate[,2] == datB[,K+1])&(g[temp]>=0.5))
  # 
  # TPR_FS11= nTruePredict/nB
  # PPV_FS11 = nTruePredict/nPredict
  
  threshold = 0.5
  index = which(g>= threshold)
  TPR_FS5 = sum(comp_mat[index,K+1])/nB
  if (length(index)==0){
    PPV_FS5 = 1
  }else{
    PPV_FS5 = sum(comp_mat[index,K+1])/length(index)
  }
  
  return(c(TPR_FS5, PPV_FS5, converge))
}

FS3 <- function(datA, datB, K, tol=1e-6, maxits = 500){
  nA = nrow(datA)
  nB = nrow(datB)
  comp_mat <- compare3(datA, datB, K=K)
  
  
  ## Using EM with the above estimated starting point
  fit = EM3(comp_mat, datA, datB, K, tol=tol, maxits = maxits)
  g = fit$g
  converge = fit$converge
  
  
  
  
  # Gmat = matrix(g, nrow = nB, byrow = TRUE)
  # Gmat[Gmat<0.5] = 0
  # opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  # predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # temp = (predict_cate[,1]-1)*nA+ predict_cate[,2]
  # # Percentage of correct link
  # nPredict = sum((g[temp]>=0.5))
  # nTruePredict = sum((predict_cate[,2] == datB[,K+1])&(g[temp]>=0.5))
  # 
  # TPR_FS11= nTruePredict/nB
  # PPV_FS11 = nTruePredict/nPredict
  
  
  
  threshold = 0.5
  index = which(g>= threshold)
  TPR_FS5 = sum(comp_mat[index,K+1])/nB
  if (length(index)==0){
    PPV_FS5 = 1
  }else{
    PPV_FS5 = sum(comp_mat[index,K+1])/length(index)
  }
  
  
  
  return(c(TPR_FS5, PPV_FS5, converge))
}


bayesian <- function(datA, datB, K){
  nB = nrow(datB)
  nA = nrow(datA)
  bayes = recordLink(datB[,1:K], datA[,1:K], eps_plus =0.01, eps_minus = 0.01,use_diff = FALSE)
  g = as.vector(t(bayes))
  
  # Gmat = bayes
  # Gmat[Gmat<0.5] = 0
  # opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  # predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # temp = (predict_cate[,1]-1)*nA+ predict_cate[,2]
  # # Percentage of correct link
  # nPredict = sum((g[temp]>=0.5))
  # nTruePredict = sum((predict_cate[,2] == datB[,K+1])&(g[temp]>=0.5))
  # 
  # TPR_bayes11= nTruePredict/nB
  # PPV_bayes11 = nTruePredict/nPredict
  
  #############################
  threshold = 0.5
  upper = which(bayes >= threshold)
  indexB = upper%%nB
  indexB[which(indexB==0)] = nB
  indexA = ceiling(upper/nB)
  idB = datB[indexB, K+1]
  idA = datA[indexA, K+1]
  
  TPR_bayes5 = sum(idA==idB)/nB
  if (length(idB) == 0){
    PPV_bayes5 = 1
  }else{
    PPV_bayes5 = sum(idA==idB)/length(idB)
  }
  
  
  
  
  return(c(TPR_bayes5, PPV_bayes5, 1))
}


############################## simulation  #########################

prev = c(0.05,0.1,0.2,0.3,0.4)
doOne1 <- function(nA, nB, K, prev, error){
  library(tictoc)

  prevalence = rep(prev,  K/5)
  data = generate_data(nA = nA, nB = nB, K = K, prevalence = prevalence, error = error)
  
  datA = data$dataA
  datB = data$dataB
  
  tic()
  FS <- FS(datA, datB, K,  tol=1e-6, maxits = 500)
  temp = toc(quiet = TRUE)
  time_FS <- temp$toc-temp$tic
  
  tic()
  FS3 <- FS3(datA, datB, K,  tol=1e-6, maxits = 500)
  temp = toc(quiet = TRUE)
  time_FS3 <- temp$toc-temp$tic
  
  #tic()
  #FS4 <- FS4(datA, datB, K,  tol=1e-6, maxits = 500)
  #temp = toc(quiet = TRUE)
  #time_FS4 <- temp$toc-temp$tic
  
  tic()
  Bayesian <- bayesian(datA, datB, K)
  temp = toc(quiet = TRUE)
  time_Bayesian <- temp$toc-temp$tic
  
  results <- c(FS, FS3, Bayesian, time_FS, time_FS3,  time_Bayesian)
  return(results)
}


vlis1 <- function(nsim=5) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500 ),
    nB = list(type="frozen", value = 200),
    K = list(type="grid", value = c(30,40,50)),
    error = list(type="grid", value = c(0.02,0.04,0.06)),
    prev = list(type="frozen", value = 0.2)
  )
  return(vList)
}

vlis2 <- function(nsim=5) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500 ),
    nB = list(type="frozen", value = 200),
    K = list(type="frozen", value = 40),
    error = list(type="frozen", value = 0.04),
    prev = list(type="grid", value = c(0.1,0.2,0.3))
  )
  return(vList)
}

vlis3 <- function(nsim=5) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="grid", value = c(400,800,1200) ),
    nB = list(type="frozen", value = 200),
    K = list(type="frozen", value = 40),
    error = list(type="frozen", value = 0.04),
    prev = list(type="frozen", value = 0.2)
  )
  return(vList)
}


#####################################
runSims <- function(vList=vlis1(),doOne = doOne1, seedList=NULL){
  res <- simsalapar::doForeach(vList,  doOne = doOne,  cluster=makeCluster(8, type="PSOCK"), seed = seedList)
  return(res)
}

#setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Binary\\Programs")
#setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Binary\\Programs")


nsim = 5

res1 <- runSims(vlis1(nsim), doOne = doOne1, seedList = 1:nsim)
save(res1, file="0722_100s_res1_binary.RData")
         
#res2 <- runSims(vlis2(nsim), doOne = doOne1)
#save(res2, file="0324_1000s_res2_binary.RData")

#res3 <- runSims(vlis3(nsim), doOne = doOne1)
#save(res3, file="0324_1000s_res3_binary.RData")




























