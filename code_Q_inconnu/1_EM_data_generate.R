
# Make error for one binary vector in databases A

makeError <- function(x,error){
  #x is a column of B
  #error is a proportion of error
  nE  = round(length(x)*error)
  index = sample(1:length(x), nE)
  x[index] = 1-x[index]
  return(x)
}

########################################################
# data generation of the bigger data base

Generate_data <- function(K,nA,nB,prevalence,error, min_prev = 0.01){
  #prevalence = rep(c(0.05,0.1,0.2,0.3,0.4), K/5)
  #error=0.04
  
  # First database B with matching variables
  datB = matrix(0, nrow = nB, ncol = K+1)
  datB[,K+1] = 1:nB #id
  
  conditionB = TRUE
  while (conditionB){
    datB[,1:K] = sapply(prevalence, function(x){rbinom(nB, size = 1,prob = x)})
    conditionB = (sum(colSums(datB[,1:K]/nB) >= min_prev) < K)
  }
  
  datB = data.frame(datB)
  colnames(datB)=c(paste("R", 1:K, sep = ""),"id") 
  
  # database A with matching variables
  conditionA = TRUE
  while (conditionA) {
    # Second database A
    idBA <- sample(1:nB,nA) #ident in A appearing in B
    datA <- datB[idBA,]
    
    # Make error for A
    datA[,1:K]= apply(datA[,1:K], MARGIN = 2, FUN = makeError, error = error)
    
    conditionA = (sum(colSums(datA[,1:K]) >= 1) < K)
  }
  
  ##covarites data
  p=2
  X1 = rnorm(nB,0,1)
  X2 = rbinom(nB,size = 1, prob = 0.7)
  X  = as.matrix(cbind(X1, X2))[idBA,]
  
  U  = runif(nA,0,1) # Time to event variable
  Tt = -log(U)/exp((X%*%beta))#lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = 0.7)
  c  = 1.777907
  
  Time      = pmin(Tt,c) # le vrai temps
  delta     = as.numeric(Tt<=c)
  surv_data = data.frame(Time,delta,X)
  
  colnames(surv_data) = c('Time','delta',paste("X", 1:p, sep = ""))
  XB = as.matrix(cbind(X1, X2))
  
  databaseB <- matrix(0,nB,(p+K+1))
  databaseB[,1:2] <- as.matrix(XB)
  databaseB [,3: ncol(databaseB) ]<-as.matrix( datB)
  databaseB <- data.frame(databaseB)
  colnames(databaseB)=c("x1","x2",paste("R", 1:K, sep = ""),"id") 
  
  databaseA <- matrix(0,nA,(p+K+1))
  databaseA[,1:2] <- as.matrix(surv_data[,1:2])
  databaseA [,3: ncol(databaseA) ]<-as.matrix( datA)
  databaseA <- data.frame(databaseA)
  colnames(databaseA)=c("Time","delta",paste("R", 1:K, sep = ""),"id") 
  
  return(list(surv_data=surv_data,XB=XB,databaseB=databaseB,
              databaseA=databaseA,datB=datB,datA=datA,idBA=idBA ) )
}


