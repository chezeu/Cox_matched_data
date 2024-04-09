library(survival)
library(plyr)
library(dplyr)
library(rootSolve) 
library(nleqslv)

library(mvtnorm)
library(glmnet)
library(matrixStats)
library(clue)

#Exponential data and exponential error
generate_data_exp_exp = function(nA, nB, K, lambdaK, error, lambdaE, round ){
  if (round == FALSE){
    XA = matrix(rexp(nA*K,rate = lambdaK),ncol = K)
    #S= sample(1:nA,nB)
    XA1 = XA[1:nB,]
    XA2  = XA[-(1:nB),]
    
    X = matrix(rbinom(nB*K, 1, error), ncol =K)
    XB = (1-X)* XA1 + X*(XA1 + matrix(rexp(nB*K, rate = lambdaE),ncol =K))
    
    datA = matrix(0,nrow = nA, ncol = K+1)
    datB = matrix(0,nrow = nB, ncol = K+1)
    datA[,1:K] = XA
    datA[,K+1] = 1:nA #id
    datB[,1:K] = XB
    datB[,K+1] = 1:nB #id
    
    
  }else if (round == TRUE){
    XA = matrix(ceiling(rexp(nA*K,rate = lambdaK)),ncol = K)
    #S= sample(1:nA,nB)
    XA1 = XA[1:nB,]
    XA2  = XA[-(1:nB),]
    
    X = matrix(rbinom(nB*K, 1, error), ncol =K)
    XB = (1-X)* XA1 + X*(XA1 + matrix(ceiling(rexp(nB*K, rate = lambdaE)), ncol =K))
    
    datA = matrix(0,nrow = nA, ncol = K+1)
    datB = matrix(0,nrow = nB, ncol = K+1)
    datA[,1:K] = XA
    datA[,K+1] = 1:nA #id
    datB[,1:K] = XB
    datB[,K+1] = 1:nB #id
  }
  
  #####
  return(list(dataA=datA, dataB = datB))
}
generate <- generate_data_exp_exp(nA, nB, K, lambdaK, error, lambdaE, round )

 datA <- generate$dataA
datB <- generate$dataB
###
compare_abs <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid(XA.k, XB.k)
    if (k == (K+1)){
      gamma.k = as.numeric(temp[,1]==temp[,2])
    }else{
      gamma.k = abs(temp[,1]-temp[,2])
    }
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  return(comp_mat)
}

compare <- compare_abs (datA, datB, K)
  


################################### Record linkage ##############################################

library(RecordLinkage)
library(plyr)
library(dplyr)
data("RLdata500")

# enlever les doublons 
generate_data_baseLinked <- function(n,m,K,error, lambdaE){
pairs=compare.dedup(RLdata500,identity=identity.RLdata500,
                    blockfld=list(c(5,6),c(6,7),c(5,7)))#prendre comme variable de blockline list
pairs=emWeights(pairs) #calculer le poids de record linkage par EM
pairs=emClassify(pairs, threshold.upper=24, threshold.lower=-7)

possibles <- getPairs(pairs, show="possible")
links=getPairs(pairs,show="links", single.rows=TRUE)
link_ids <- links[, c("id1", "id2")]

a <- link_ids[,1]
database <- RLdata500[-a,]
######### ajouter des erreurs pour la seconde base sur les variable numerique 

#### nB=m, nA=n

PA = database[,5:7] 
#S= sample(1:nA,nB)
PA1 = PA[1:n,]
PA2  = PA[-(1:n),]

P = matrix(rbinom(m*K, 1, error), ncol =K)

donnee_A = (1-P)* PA1 + P*(PA1 + matrix(ceiling(rexp(n*K, rate = lambdaE)), ncol =K))

databaseA <- database[1:n,]
databaseA[,5:7] <- donnee_A
databaseA$idA <- 1:nrow(databaseA) 

databaseB <- database[1:m,]

databaseB$idB <- 1:nrow(databaseB)

true_link <- databaseB[1:n,]

return(list( true_link=true_link, databaseA=databaseA, databaseB=databaseB ))
 }

 gen <- generate_data_baseLinked(n,m,K,error, lambdaE)



########################################################

Generate_data <- function(m,n,beta,p){
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

surv_data <- Generate_data (m,n,beta,p)

#############################################################
#blocking
#p <- pair_blocking (databaseA, databaseB, "fname_c1" , FALSE) 
#######################"
databaseA <- gen$databaseA
databaseB <- gen$databaseB
library(reclin2)
attach(databaseA)
attach(databaseB)

linkage_data <- function(databaseB, databaseA) {

  surv_data <- Generate_data(m,n,beta,p)
  X <- surv_data[,3:ncol(surv_data)]
  XB <- as.matrix(X, ncol = p) #Only covariates
  
  ###### True-link data: The first n individual from B
  data_true <- surv_data[1:n,]
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  
  databaseA$Time <- data_A$Time
  databaseA$delta <- data_A$delta
  dat <- data.frame(XB)
  databaseB$x1 <- dat$X1
  databaseB$x2 <- dat$X2
  
pp <- pair_minsim(databaseA, databaseB, 
                        on =  c("by","bm","bd"), minsim = 0)
#pp <- print(pp)

# show the pairs of comparison en laissant la variable 
#p <- compare_pairs(p, on= c("fname_c1","lname_c1","by","bm","bd"), inplace = TRUE, 
#             comparators = list(fname_c1 = jaro_winkler(), lname_c1  = jaro_winkler() ))

pp <- compare_pairs(pp, on= c("by","bm","bd"), inplace = TRUE)
#p <- print(p)

## classification
model <- problink_em ( ~ by +bm + bd, data=pp)
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
  
  surv_data <- Generate_data(m,n,beta,p)
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

data_naive <- databaseA[,9:10]
covar_naive_B <- databaseB[lB,]

data_naive$x1 <- covar_naive_B$x1
data_naive$x2 <-  covar_naive_B$x2


return(list(data_naive=data_naive, Naive_max=Naive_max))
}

Naive <- generate_naive_data(pp,databaseB, databaseA)
naive_max <- Naive$Naive_max
data_naive<- Naive$data_naive


#########  multinomiale Q
##(estimate==true) pour dire qu'on a generé la matrice de proba

Matrix_function <-function(n,m,C,beta,databaseB, databaseA, estimate){

  surv_data <- Generate_data(m,n,beta,p)
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
################## ##################################################


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

m=40
n=20
K=3
error=40/100 #pourcentage d'erreur
lambdaE=1/5# erreur ajouté à chaque element de la seconde base
round=TRUE
estimate=FALSE
lambdaK=1/50 # parametre de la loi expo

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
                coef_estimat=coef_estimat,converge_estimat =converge_estimat,  converge_aee1= converge_aee1))
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


coef_aee1 <- matrix(0,nrow = nsim, ncol = p)
converge_aee1 <- vector()

test_aee1 <- function(nsim, C,alpha){
  for (i in 1:nsim) {
    g <-  doOne2d_equat_estimat(n,C,alpha,estimate=FALSE)
    
    coef_aee1[i,] <- g$coef_aee1
    converge_aee1[i] <- g$converge_aee1
  }
  
  L3 <- which(converge_aee1!=0)
  coef_aee1<- coef_aee1[L3,]
  result_aee1 <- colMeans(coef_aee1)
  
  return(coef_aee1=coef_aee1)
}

result_aee1 <- test_aee1(100,C,alpha)
coef_aee11<-colMeans(result_aee1)


#######################################################

coef_naive<- matrix(0,nrow = nsim, ncol = p)
converge_naive <- vector()

test_naive <- function(nsim, C,alpha){
  for (i in 1:nsim) {
    g <-  doOne2d_equat_estimat(n,C,alpha,estimate=FALSE)
    
    coef_naive[i,] <- g$coef_naive2
    converge_naive[i] <- g$converge_naive2
  }
  
  L2 <- which(converge_naive!=0)
  coef_naive<- coef_naive[L2,]
  result_naive <- colMeans(coef_naive)
  
  return(coef_naive=coef_naive)
}

result_naive <- test_naive(100,C,alpha)
coef_naive1 <-colMeans(result_naive)

 




