
#########  multinomiale Q

Matrix_function<-function(n,m,C,beta, surv_data){
  
  #surv_data <- Generate_data(m,n,beta)
  X  =  surv_data[,3:(p+2)]
  XB = as.matrix(X, ncol = p) #Only covariates
  
  ###### True-link data: The first n individual from B
  data_true = surv_data[1:n,]
  data_A = surv_data[1:n,1:2] #Survival time and censored indicator
  
  
  Q = matrix(0,n,m)
  Lvec = matrix(0,n,m)
  for (i in 1:n) {
    if(i!=1 & i!=n){ 
      pro = rep(0,m)
      pro[i] = C[1] #Probability of True links
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
  X_naive <-Lvec %*% XB
  colnames(X_naive)<- c("Z1","Z2")
  
  # Naive data
  data_naive <- cbind( data_A, X_naive)
  
  return(list( data_true=data_true, data_A= data_A,XB=XB, data_naive=data_naive,Q=Q, Lvec= Lvec))
}