library(clue)
library(matrixStats)
library(ludic)
library(dplyr)
library(simsalapar)
library(doParallel)
library(mvtnorm)
library(mixtools)

rm(list = ls())

#datA, datB: two databases A and B, B is smaller database, 
#datA, datB include K matching variables and 1 ID column
#K1, K2, K3: number of matching variables respected to binary, categorical and continuous comparisonvalue
#K1+K2 = K
# e3: noise in data B
gendata = function(nA, nB, K3, e3){
  datA = matrix(0,nrow = nA, ncol = K3+1)
  datA[,K3+1] = 1:nA
  
  # K3 date variable
  datA[,1:K3] <- matrix(sample(10:2000, size=nA*K3, replace = TRUE), ncol =K3)
  
  #datA[,1:K3] <- mvtnorm::rmvnorm(nA, mean = rep(0,K3))
  
  datA = data.frame(datA)
  colnames(datA)=c(paste("R", 1:K3, sep = ""),"id") 
  
  # Second database B
  idAB <- sample(1:nA,nB) #ident in A appearing in B
  datB <- datA[idAB,]
  
  # Make error for B
  makeError_date <- function(x,e){
    #x is a column of B
    #e is a proportion of error
    nE  = round(length(x)*e)
    index = sample(1:length(x), nE)
    x[index] = x[index] + sample(-3:3, size = nE, replace = TRUE) #+- 3 days
    return(x)
  }
  datB[,1:K3] = apply(datB[,1:K3], MARGIN = 2, FUN = makeError_date, e = e3)
  
  
  #####
  return(list(dataA=datA, dataB = datB))
}


compare_date = function(datA, datB, K3, type){
  #type = 1: distance for each variable seperately
  #type = 2: distance for all vector of matching variables
  nA = nrow(datA)
  nB = nrow(datB)
  N =  nA*nB #total record pairs
  
  # Compare one pair
  compare_one = function(b, a, K3, type){
    b = as.numeric(b)
    a = as.numeric(a)
    
    # numeric comparison part
    if (type==1){
      comp = rep(0,K3+1)
      comp[K3+1] = as.numeric(a[K3+1]==b[K3+1]) #Match indicator
      comp[1:K3] = abs(a[1:K3]-b[1:K3])
    }else if (type==2){
      comp = rep(0,2)
      comp[2]=as.numeric(a[K3+1]==b[K3+1]) #Match indicator
      comp[1]=dist(rbind(a[1:K3],b[1:K3]), method = "euclidean")
    }
    return(comp)
  }
  
  require(doParallel)
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)  
  registerDoParallel(cl) 
  comp_mat <- foreach(ib = 1:nB,.combine = rbind)%:%foreach(ia=1:nA,.combine = rbind) %dopar% {
    compare_one(datB[ib,],datA[ia,],K3, type)
  }
  stopCluster(cl)
  return(comp_mat)
}




# type = 2 -> there is only 1 comparison variable
# EM function for type 2
EM_cont2 = function(X, datA, datB, K3, int_mean, int_sd, tol = 1e-5, maxits = 200){
  #int_mean = c(mean_M, mean_U), initial values for mean..
  nA = nrow(datA)
  N = nrow(X)

  X_cont <- X[,1]
  # Fixed starting point, Normal distribution
  p = 1/nA
  mean_M = int_mean[1] 
  mean_U = int_mean[2]
  sd_M = int_sd[1]
  sd_U = int_sd[2]
  
  # initializations
  g = rep(0,N) # probability of being in Match  for each pair l
  it = 0
  converged = FALSE
  
  while ((!converged) & (it < maxits)){ 
    pOld = p
    
    mean_MOld = mean_M
    mean_UOld = mean_U
    
    sd_MOld = sd_M
    sd_UOld = sd_U
    ### E
    # Compute expectation
    probM <- dnorm(X_cont, mean = mean_M, sd = sd_M)
    probU <- dnorm(X_cont, mean = mean_U, sd = sd_U)
    g <- (p*probM)/(p*probM+(1-p)*probU)
    
    
    ### Maximization
    
    mean_M = as.vector((X_cont%*%g)/sum(g))
    mean_U = as.vector((X_cont%*%(1-g))/sum(1-g))
    sd_M = sqrt((((X_cont-mean_M)^2)%*%g)/sum(g))
    sd_U = sqrt((((X_cont-mean_U)^2)%*%(1-g))/sum(1-g))
    p = sum(g)/N
    
    if(sd_M < 0.001){
      sd_M <- 0.001
    }
    
    paraOld = c(mean_MOld,mean_UOld,sd_MOld,sd_UOld,pOld)
    paraCurrent = c(mean_M, mean_U, sd_M, sd_U, p)
      
    it = it + 1
    
    converged = max(abs(paraOld - paraCurrent)/abs(paraOld)) <= tol 
  }
  
  probM <- dnorm(X_cont, mean = mean_M, sd = sd_M)
  probU <- dnorm(X_cont, mean = mean_U, sd = sd_U)
  g <- (p*probM)/(p*probM+(1-p)*probU)
    
  return(list(g =g, p= p,mean_M = mean_M, mean_U=mean_U, sd_M=sd_M, sd_U=sd_U, it=it))
}



# nA = 200
# nB = 100
# K3 = 3
# e3 = 0.1
# data = gendata(nA, nB, K3, e3)
# datA = data$dataA
# datB = data$dataB
# 
# ## Without scale
# comp_mat <- compare_date(datA, datB, K3, type=2)
# 
# X_cont <- comp_mat[,1]
# 
# layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
# hist(X_cont, freq = FALSE,  main = "Histogram of all pairs without scale", xlab = "Distance")
# lines(density(X_cont, from = 0), col = "red")
# 
# 
# hist(X_cont[comp_mat[,2]==1], freq = FALSE,  main = "Histogram of Matched pairs without scale", xlab = "Distance")
# lines(density(X_cont[comp_mat[,2]==1], from = 0), col = "red")
# 
# hist(X_cont[comp_mat[,2]==0], freq = FALSE, main = "Histogram of Unmatched pairs without scale", xlab = "Distance")
# lines(density(X_cont[comp_mat[,2]==0], from = 0), col = "red")
# 
# 
# ###### with scale 
# datA_scale <- datA
# datB_scale <- datB
# #datA_scale[,1:K3] <- scale(datA[,1:K3])
# #datB_scale[,1:K3] <- scale(datB[,1:K3], colMeans(datA[,1:K3]),colSds(as.matrix(datA[,1:K3])))
# 
# datA_scale[,1:K3] <- scale(datA[,1:K3])
# datB_scale[,1:K3] <- scale(datB[,1:K3])

# comp_mat_scale <- compare_date(datA_scale, datB_scale, K3, type = 2)
# X_cont_scale <- comp_mat_scale[,1]
# 
# layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
# hist(X_cont_scale, freq = FALSE,  main = "Histogram of all pairs with scale", xlab = "Distance")
# lines(density(X_cont_scale, from = 0), col = "red")
# 
# 
# hist(X_cont_scale[comp_mat_scale[,2]==1], freq = FALSE,  main = "Histogram of Matched pairs with scale", xlab = "Distance")
# lines(density(X_cont_scale[comp_mat_scale[,2]==1], from = 0), col = "red")
# 
# hist(X_cont_scale[comp_mat_scale[,2]==0], freq = FALSE, main = "Histogram of Unmatched pairs with scale", xlab = "Distance")
# lines(density(X_cont_scale[comp_mat_scale[,2]==0], from = 0), col = "red")
# 
# 
# 
# ### Mixture of normal distribution
# # true paramater
# truemean_M <- mean(X_cont[comp_mat[,2]==1])
# truemean_U <- mean(X_cont[comp_mat[,2]==0])
# truesd_M <- sd(X_cont[comp_mat[,2]==1])
# truesd_U <- sd(X_cont[comp_mat[,2]==0])
# c(1/nA, truemean_M, truemean_U, truesd_M, truesd_U)
# 
# #
# int_mean <- c(0.5, mean(X_cont))
# int_sd <- c(0.2, sd(X_cont))
# fit <- EM_cont2(comp_mat, datA, datB, K3, int_mean, int_sd)
# g = fit$g
# g[g<0]=0
# Gmat = matrix(g, nrow = nB, byrow = TRUE)
# opti_cont <- solve_LSAP(Gmat, maximum = TRUE)
# predict_cont = cbind(seq_along(opti_cont), opti_cont)
# 
# # Percentage of correct link
# TPR_cont = sum(predict_cont[,2] == datB[,K3+1])/nB
# TPR_cont
# # Using package mixtool to compare
# a <- normalmixEM(X_cont, lambda = c(1/nA,1-1/nA), mu = int_mean, sigma = int_sd)
# c(a$lambda[1], a$mu, a$sigma)
# aa <- normalmixEM(X_cont, lambda = c(1/nA,1-1/nA), mu = c(truemean_M, truemean_U), sigma=c(truesd_M, truesd_U))
# c(aa$lambda[1], aa$mu, aa$sigma)

###########################
check <- function(datA, datB, K3, mean_M, sd_M){
  nB = nrow(datB)
  nA = nrow(datA)
  comp_mat <- compare_date(datA, datB, K3, type=2)
  
  X_cont <- comp_mat[,1]
  
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  hist(X_cont, freq = FALSE,  main = "Histogram of all pairs ", xlab = "Distance")
  lines(density(X_cont, from = 0), col = "red")
  
  
  hist(X_cont[comp_mat[,2]==1], freq = FALSE,  main = "Histogram of Matched pairs ", xlab = "Distance")
  lines(density(X_cont[comp_mat[,2]==1], from = 0), col = "red")
  
  hist(X_cont[comp_mat[,2]==0], freq = FALSE, main = "Histogram of Unmatched pairs ", xlab = "Distance")
  lines(density(X_cont[comp_mat[,2]==0], from = 0), col = "red")
  
  int_mean <- c(mean_M, mean(X_cont))
  int_sd <- c(sd_M, sd(X_cont))
  #fit <- EM_cont2(comp_mat, datA, datB, K3, int_mean, int_sd)
  fit <-  normalmixEM(X_cont, lambda = c(1/nA,1-1/nA), mu = int_mean, sigma = int_sd)
  g = fit$g
  g[g<0]=0
  Gmat = matrix(g, nrow = nB, byrow = TRUE)
  opti_cont <- solve_LSAP(Gmat, maximum = TRUE)
  predict_cont = cbind(seq_along(opti_cont), opti_cont)
  
  # Percentage of correct link
  TPR_cont = sum(predict_cont[,2] == datB[,K3+1])/nB
  return(TPR_cont)
}

nA = 500
nB = 100
K3 = 3
e3 = 0.1
data = gendata(nA, nB, K3, e3)
datA = data$dataA
datB = data$dataB

datA_scale <- datA
datB_scale <- datB
#datA_scale[,1:K3] <- scale(datA[,1:K3])
#datB_scale[,1:K3] <- scale(datB[,1:K3], colMeans(datA[,1:K3]),colSds(as.matrix(datA[,1:K3])))

datA_scale[,1:K3] <- scale(datA[,1:K3])
datB_scale[,1:K3] <- scale(datB[,1:K3])

check(datA, datB, K3, mean_M = 0.5, sd_M = 1)
check(datA_scale, datB_scale, K3, mean_M = 0.2, sd_M = 0.1)
