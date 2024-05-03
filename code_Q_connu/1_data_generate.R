
########################################################
# data generation of the bigger data base
Generate_data <- function(m,n,beta,sigma,alpha){
  # set.seed(42)
  ##covarites data
  X1 = rnorm(m,0,1)
  X2 = rbinom(m,size = 1, prob = 0.8)
  X  = as.matrix(cbind(X1, X2))
  
  U  = runif(m,0,1) # Time to event variable
  Tt = sigma* (-log(U)/exp((X%*%beta)))^(1/alpha) #lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = 0.75)
  c  = 1.777907
  
  Time      = pmin(Tt,c) # le vrai temps
  delta     = as.numeric(Tt<=c)
  surv_data = data.frame(Time,delta,X)
  colnames(surv_data) = c('Time','delta',paste("X", 1:2, sep = ""))
  
  return(surv_data=surv_data )
}

######### Q matrix construction 

Matrix_function<-function(n,m,C,beta, surv_data){
  
  X  =  surv_data[,3:(p+2)]
  XB = as.matrix(X, ncol = p) #Only covariates
  
  ###### True-link data: The first n individual from B
  data_true = surv_data[1:n,]
  data_A = surv_data[1:n,1:2] #Survival time and censored indicator
  
  Q = matrix(0,n,m)
  Lvec = matrix(0,n,m) #multinomial generation of the naive covariable
  for (i in 1:n) {
    if(i!=1 & i!=n){ 
      pro = rep(0,m)
      pro[i] = C[1] #Probability of True links (the highest proba)
      pro[i+1]= C[3]
      pro[i-1]=C[3]
      Q[i,]<-pro
      Lvec[i,] <- rmultinom(1, size = 1, pro = pro)} 
    if(i==1)  {
      pro = rep(0,m)
      pro[i] =  C[1] #Probability of True links
      pro[i+1]= C[2]
      pro[i+2]= C[3]
      Q[i,]<-pro
      Lvec[i,] <- rmultinom(1, size = 1, pro = pro)} 
    if(i==n){ 
      pro = rep(0,m)
      pro[i] = C[1] #Probability of True links
      pro[i-1]=C[2]
      pro[i-2]=C[3]
      Q[i,]<-pro
      Lvec[i,] <- rmultinom(1, size = 1, pro= pro) }  
  }
  
  Lvec
  X_naive <-Lvec %*% XB # the naive covariables 
  colnames(X_naive)<- c("Z1","Z2")
  
  # Naive data
  data_naive <- cbind( data_A, X_naive)
  
  return(list( data_true=data_true, data_A= data_A,XB=XB, data_naive=data_naive,Q=Q, Lvec= Lvec))
}


