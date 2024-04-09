
makeError_date <- function(x,e){
  #x is a column of B
  #e is a proportion of error
  nE  = round(length(x)*e)
  index = sample(1:length(x), nE)
  x[index] = x[index] + sample(c(-7:-1,1:7), size = nE, replace = TRUE) #+- 7 days
  return(x)
}

generate_data = function(nA, nB, K, e, max_date){
  datA = matrix(0,nrow = nA, ncol = K+1)
  datA[,K+1] = 1:nA #id
  
  # K date variable
  datA[,1:K] <- matrix(sample(10:max_date, size=nA*K, replace = TRUE), ncol =K)
  
  datA = data.frame(datA)
  colnames(datA)=c(paste("R", 1:K, sep = ""),"id") 
  
  # Second database B
  idAB <- sample(1:nA,nB) #ident in A appearing in B
  datB <- datA[idAB,]
  
  # Make error for B
  datB[,1:K] = apply(datB[,1:K], MARGIN = 2, FUN = makeError_date, e = e)
  
  
  #####
  return(list(dataA=datA, dataB = datB))
}


library(doParallel)

compare_cont = function(datA, datB, K){
  # Compare one pair
  compare_one = function(b, a, K){
    #type = 2: vector distance for all  matching variables
    b = as.numeric(b)
    a = as.numeric(a)
    
    comp = rep(0,2)
    comp[2]=as.numeric(a[K+1]==b[K+1]) #Match indicator
    comp[1]=dist(rbind(a[1:K],b[1:K]), method = "euclidean")
    
    return(comp)
  }
  
  nA = nrow(datA)
  nB = nrow(datB)
  N =  nA*nB #total record pairs
  
  require(doParallel)
  no_cores <- 8
  cl <- makeCluster(no_cores)  
  registerDoParallel(cl) 
  comp_mat <- foreach(ib = 1:nB,.combine = rbind)%:%foreach(ia=1:nA,.combine = rbind) %dopar% {
    compare_one(datB[ib,],datA[ia,],K)
  }
  stopCluster(cl)
  return(comp_mat)
}

compare_binary = function(datA, datB, K){
  # Compare one pair
  compare_one = function(b, a, K){
    #simple binary comparison (agreement/disagreement)
    b = as.numeric(b)
    a = as.numeric(a)
    
    comp = as.numeric(a==b)
    return(comp)
  }
  
  nA = nrow(datA)
  nB = nrow(datB)
  N =  nA*nB #total record pairs
  
  require(doParallel)
  no_cores <- 8
  cl <- makeCluster(no_cores)  
  registerDoParallel(cl) 
  comp_mat <- foreach(ib = 1:nB,.combine = rbind)%:%foreach(ia=1:nA,.combine = rbind) %dopar% {
    compare_one(datB[ib,],datA[ia,],K)
  }
  stopCluster(cl)
  return(comp_mat)
}

nA = 100
nB = 50
K = 3
e = 0.1
max_date = 50
data = generate_data(nA, nB, K, e, max_date )
datA = data$dataA
datB = data$dataB
comp_mat <- compare_cont(datA, datB, K )
X <- as.vector(comp_mat[,1])+0.01

library(mixtools)
lambda = 1/nA
fit <- gammamixEM(X, lambda = c(1/nA, 1-1/nA), verb = TRUE)

library(mixR)
fit <- mixfit(X, ncomp = 2, family = "gamma", init.method = "kmeans")
