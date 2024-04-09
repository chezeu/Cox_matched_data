
# Make error for one binary vector x
makeError <- function(x,error){
  #x is a column of B
  #error is a proportion of error
  nE  = round(length(x)*error)
  index = sample(1:length(x), nE)
  x[index] = 1-x[index]
  return(x)
}

generate_data <- function(nA, nB, K, prevalence, error, min_prev = 0.01){
  #prevalence = rep(c(0.05,0.1,0.2,0.3,0.4), K/5)
  
  # First database A
  datA = matrix(0, nrow = nA, ncol = K+1)
  datA[,K+1] = 1:nA #id
  
  conditionA = TRUE
  while (conditionA){
    datA[,1:K] = sapply(prevalence, function(x){rbinom(n=nA, size = 1,prob = x)})
    conditionA = (sum(colSums(datA[,1:K]/nA) >= min_prev) < K)
  }
  
  datA = data.frame(datA)
  colnames(datA)=c(paste("R", 1:K, sep = ""),"id") 
  
  conditionB = TRUE
  while (conditionB) {
    # Second database B
    idAB <- sample(1:nA,nB) #ident in A appearing in B
    datB <- datA[idAB,]
    
    # Make error for B
    datB[,1:K]= apply(datB[,1:K], MARGIN = 2, FUN = makeError, error = error)
    
    conditionB = (sum(colSums(datB[,1:K]) >= 1) < K)
  }
  
  
  return(list(dataA=datA, dataB = datB, prev = prevalence))
}

library(doParallel)
##################################################


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
#############################################

library(clue)
library(matrixStats)
library(klaR)


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

#################################################


library(doParallel)
library(clue)
library(matrixStats)
library(klaR)
library(ludic)


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


#####################################################################



m=500
n=200
K=30
p=2

error=0.1# pas de convergence lorsque error=0.1 
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




