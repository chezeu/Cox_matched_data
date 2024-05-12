
########################################################
# data generation of the bigger data base
Generate_data <- function(n0,n,beta,sigma,alpha){
  # set.seed(42)
  ##covarites data
  X1 = rnorm(n,0,1)
 # X1 = runif(n,0,1)
  X2 = rbinom(n,size = 1, prob = 0.7)
  X  = as.matrix(cbind(X1, X2))
  
  U  = runif(n,0,1) # Time to event variable
  Tt = sigma*(-log(U)/exp((X%*%beta)))^(1/alpha) #lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = 0.7)
   c  = 1.777764
  #c= 1.344032
  
  Time = pmin(Tt,c) # le vrai temps
  delta = as.numeric(Tt<=c)
  surv_data = data.frame(Time,delta,X)
  colnames(surv_data) = c('Time','delta',paste("X", 1:2, sep = ""))
  # data XB 
  #X3=rnorm(n0,0,1)
  X3 = runif(n0,0,1)
  X4 = rbinom(n0,size = 1, prob = 0.7)
  XB = rbind(X,cbind(X3,X4))
  
  return(list(surv_data=surv_data, XB=XB))
}

######### Q matrix construction 

Matrix_function<-function(C,XB=XB,surv_data){
  
  #C=c(0.8,0.1,0.1)
  ###### True-link data
  data_true = surv_data
  data_A = surv_data[,1:2] #Survival time and censored indicator
  
  n = nrow(data_A)
  m = nrow(XB)
  Q = matrix(0,n,m)
  L = matrix(0,n,m) #multinomial generation of the naive covariable
  for (i in 1:n) {
    if(i==1){Q[i,1:3] = C
    }else if(i==n){
      Q[i,c(n-2,n-1,n)] = C
    }else{
    Q[i,(i+2)] = C[2]
    Q[i,(i+3)] = C[1]
    Q[i,(i+4)] = C[3]
    }
  }
  
  Lvec <- sapply(1:n,FUN = function(i){
    return(sample(1:m, size = 1, prob = Q[i,]))
  })
  Lvec 
  L[cbind(1:n, Lvec)] = 1
  X_naive <-L %*% XB # the naive covariables 
  colnames(X_naive)<- c("Z1","Z2")
  
  # Naive data
  data_naive <- cbind( data_A, X_naive)
  
  return(list( data_true=data_true, data_A= data_A, data_naive=data_naive,Q=Q, L= L))
}


