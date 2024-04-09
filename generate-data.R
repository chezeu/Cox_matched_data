library(mvtnorm)
library(glmnet)
library(matrixStats)
library(clue)
# This function generate survival data
# Sample size: N
# Number of Gaussian covariates: p
# True values of coeficients: beta
# Proportion of non-censoring: rho

# Genrate survival time with one normal covariate and one binary covariate
generate_surv_mixed = function(N, p, beta, rho ){
  
  X1 <- rbinom(N,size = 1, prob = 0.7)
  
  # Time to event variable
  U <- runif(N,0,1)
  Tt = -log(U)/exp((X1*beta))  #lambda=1
  
  # Constant right censored
  C = 0.6577907
  
  Time = pmin(Tt, C)
  delta = as.numeric(Tt <= C)
  #print(sum(delta)/n)
  
  # Surival data
  surv_data= data.frame(Time, delta, X1)
  colnames(surv_data)=c('Time','delta',paste("X1", 1:p, sep = ""))
  
  return(list(surv_data=surv_data, duration = Tt))
}

# This function generates two databases A and B with
# n: Number on individuals in A
# N: Number of individuals in B
# p: Number of covariates
# beta: True coefficients
# rho: Proportion of non-censoring
# alpha_q: Record linkage quality

generate_linked_mixed <- function(n, N, p, beta, rho , alpha_q, informative){
  
  # Survival data (Database B)
  surv_mixed <- generate_surv_mixed(N, p, beta, rho)
  surv_data <- surv_mixed$surv_data
  X <- surv_data[,3:ncol(surv_data)]
  
  Time <- surv_mixed$duration
  XB <- as.matrix(X, ncol = p) #Only covariates
  
  temp <- generate_surv_mixed(1000000, p, beta, rho)
  qt <- quantile(temp$duration,1-alpha_q)
  
  ###### True-link data: The first n individual from B
  data_true <- surv_data[1:n,]
  
  
  ###### Linked data
  data_A <- surv_data[1:n,1:2] #Survival time and censored indicator
  
  if (informative==1){ #dependent on X
    Lvec <- sapply(1:n,FUN = function(i){
      if (X[i,1]>=qnorm(1-alpha_q)){
        j <- i
      }else{
        j <- sample(1:N,size=1)
      }
      return(j)
    })
  }else if (informative ==2){ #dependent on survival time
    Lvec <- sapply(1:n,FUN = function(i){
      if (Time[i]>=qt){
        j <- i
      }else{
        j <- sample(1:N,size=1)
      }
      return(j)
    })
  }
  
  
  
  Lvec 
  
  L = matrix(0, nrow = n, ncol =N)
  L[cbind(1:n, Lvec)] = 1
  
  X_naive <- L %*% XB
  
  # Naive data
  data_naive <- cbind(data_A, X_naive)
  trueLink_id <- diag(L)
  
  return(list(data_true = data_true, data_A = data_A,   XB = XB, data_naive = data_naive, trueLink_id = trueLink_id))
}

