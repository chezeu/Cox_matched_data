library(survival)
library(rootSolve) 
library(nleqslv)
setwd( "C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_EM_method")

source("1_data_generate.R")

library(mvtnorm)
library(glmnet)
library(matrixStats)
library(clue)
library(klaR)
library(ludic)
library(RecordLinkage)

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

#################
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

surv_data <- Generate_surv_data (m,n,beta,p)

generate_naive_data <- function(comp_mat.1,databaseB, databaseA,surv_data ){
  
    
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


