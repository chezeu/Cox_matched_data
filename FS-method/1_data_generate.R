


library(RecordLinkage)
library(reclin2)
#library(plyr)
#library(dplyr)

data("RLdata500")
#data("RLdata10000")

s1 <- 1:nB
datasetB <- RLdata500[s1,] #matching variables
# Second database A
idBA <- sample(1:nB,nA) #ident in A appearing in B
datasetA <- datasetB[idBA,]

var_block= "fname_c1"

matching_variables = colnames(datasetA)

########################################################
# data generation of the database

Generate_data <- function(nA,nB,idBA,beta){
  
  # First database B 
  datB = matrix(0, nrow = nB, ncol = p+1)
  datB[,p+1] = 1:nB #id
  
  ##covarites data
  p=3
  X1 = rnorm(nB,0,1)
  X2 = rbinom(nB,size = 1, prob = 0.7)
  X3 = rnorm(nB,0,2)
  X  = as.matrix(cbind(X1, X2, X3))[idBA,]
   
  #datA <- datB[idBA,]
  
  U  = runif(nA,0,1) # Time to event variable
  Tt = -log(U)/exp((X%*%beta))#lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = 0.7)
  c=  2.948162
  
  Time      = pmin(Tt,c) # le vrai temps
  delta     = as.numeric(Tt<=c)
  surv_data = data.frame(Time,delta,X)
  
  colnames(surv_data) = c('Time','delta',paste("X", 1:p, sep = ""))
  XB = as.matrix(cbind(X1, X2,X3))
  
  datB[,1:3] <- as.matrix(XB)
  datB <- data.frame(datB)
  colnames(datB)=c("x1","x2","X3","id") 
  
  datA <- matrix(0,nA,3)
  datA[,1:2] <- as.matrix(surv_data[,1:2])
  datA[,3] <- idBA 
  datA <- data.frame(datA)
  colnames(datA)=c("Time","delta","id") 
  
  return(list(surv_data=surv_data,XB=XB,datB=datB,datA=datA ) )
}

#mf = Generate_data (nA,nB,idBA,beta)

