

########################################################
# data generation of the bigger data base

Generate_data <- function(K,nA,nB,prevalence,censor, min_prev = 0.01){
  
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

    # Second database A
    idBA <- sample(1:nB,nA) #ident in A appearing in B
    datA <- datB[idBA,]

  ##covarites data
  p=3
  X1 = rnorm(nB,0,1)
  X2 = rbinom(nB,size = 1, prob = 0.7)
  X3 = rnorm(nB,0,2)
  X  = as.matrix(cbind(X1, X2, X3))[idBA,]
  
  U  = runif(nA,0,1) # Time to event variable
  Tt = -log(U)/exp((X%*%beta))#lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = 0.7)
  #c  = 1.777907
  c=  censor
  
  Time      = pmin(Tt,c) # le vrai temps
  delta     = as.numeric(Tt<=c)
  surv_data = data.frame(Time,delta,X)
  
  colnames(surv_data) = c('Time','delta',paste("X", 1:p, sep = ""))
  XB = as.matrix(cbind(X1, X2,X3))
  
  databaseB <- matrix(0,nB,(p+K+1))
  databaseB[,1:3] <- as.matrix(XB)
  databaseB [,4: ncol(databaseB) ]<-as.matrix( datB)
  databaseB <- data.frame(databaseB)
  colnames(databaseB)=c("x1","x2","X3",paste("R", 1:K, sep = ""),"id") 
  
  databaseA <- matrix(0,nA,(2+K+1))
  databaseA[,1:2] <- as.matrix(surv_data[,1:2])
  databaseA [,3:ncol(databaseA) ]<-as.matrix( datA)
  databaseA <- data.frame(databaseA)
  colnames(databaseA)=c("Time","delta",paste("R", 1:K, sep = ""),"id") 
  
  return(list(surv_data=surv_data,XB=XB,databaseB=databaseB,
              databaseA=databaseA,datB=datB,datA=datA,idBA=idBA ) )
}


