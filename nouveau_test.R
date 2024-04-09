
library(mvtnorm)
library(glmnet)
library(matrixStats)
library(clue)

# This function generate survival data
# Sample size: N
# Number of Gaussian covariates: p
# True values of coeficients: beta
# Proportion of non-censoring: rho

# # Monte-Carlo estimation of constant censoring time that yields 30% censoring
# censoring <- sapply(1:100000, FUN = function(i){
#   N = 1000
#   beta = c(0.5,-0.5)
#   rho = 0.7 #proportion of non-censoring
#   
#   
#   X1 <- rnorm(N,0,1)
#   X2 <- rbinom(N,size = 1, prob = 0.7)
#   X <- as.matrix(cbind(X1, X2))
# 
# 
#   # Time to event variable
#   U <- runif(N,0,1)
#   
#   Tt = -log(U)/exp((X%*%beta))  #lambda=1
# 
#   # Constant right censored
#   C = quantile(Tt, probs = rho)
#   return(C)
# })

# Genrate survival time with one normal covariate and one binary covariate
generate_surv_mixed = function(N, p, beta, rho ){
  
  X1 <- rnorm(N,0,1)
  X2 <- rbinom(N,size = 1, prob = 0.7)
  X <- as.matrix(cbind(X1, X2))
  
  
  # Time to event variable
  U <- runif(N,0,1)
  Tt = -log(U)/exp((X%*%beta))  #lambda=1
  
  # Constant right censored
  #C = quantile(Tt, probs = rho)
  C = 1.777907 #obtained from above Monte-Carlo estimation
  
  Time = pmin(Tt, C)
  delta = as.numeric(Tt <= C)
  #print(sum(delta)/n)
  
  # Surival data
  surv_data= data.frame(Time, delta, X)
  colnames(surv_data)=c('Time','delta',paste("X", 1:p, sep = ""))
  
  return(surv_data)
}

# This function generates two databases A and B with
# n: Number on individuals in A
# N: Number of individuals in B
# p: Number of covariates
# beta: True coefficients
# rho: Proportion of non-censoring
# alpha_q: Record linkage quality

generate_linked_mixed <- function(n, N, p, beta, rho , alpha_q){
  
  # Survival data (Database B)
  surv_data <- generate_surv_mixed(N, p, beta, rho)
  X <- surv_data[,3:ncol(surv_data)]
  XB <- as.matrix(X, ncol = p) #Only covariates
  
  ###### True-link data: The first n individual from B
  data_true <- surv_data[1:n,]
  
  
  ###### Linked data
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  
  Lvec <- sapply(1:n,FUN = function(i){
    prob = rep((1-alpha_q)/N,N)
    prob[i] = alpha_q + (1-alpha_q)/N #Probability of True links
    return(sample(1:N, size = 1, prob = prob))
  })
  Lvec 
  
  L = matrix(0, nrow = n, ncol =N)
  L[cbind(1:n, Lvec)] = 1
  
  X_naive <- L %*% XB
  
  # Naive data
  data_naive <- cbind(data_A, X_naive)
  trueLink_id <- diag(L)
  
  return(list(data_true = data_true, data_A = data_A,   XB = XB, data_naive = data_naive, trueLink_id = trueLink_id))
}

