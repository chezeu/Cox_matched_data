
########################################################
# data generation of the bigger data base
Generate_data <- function(nB,nA,beta,Q){
  ##covarites data
  p=2
  X1 = rnorm(nB,0,1)
  X2 = rbinom(nB,size = 1, prob = 0.7)
  X  = as.matrix(cbind(X1, X2))
  
  L_true = sapply(1:nA, FUN = function(j){
    return( rmultinom(1,size=1, prob = Q[j,]))
  } )
  X_true = t(L_true)%*%X
  
  U  = runif(nA,0,1) # Time to event variable
  Tt = -log(U)/exp((X_true%*%beta))#lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = 0.7)
  c  = 1.777907
  
  Time      = pmin(Tt,c) # le vrai temps
  delta     = as.numeric(Tt<=c)
  surv_data = data.frame(Time,delta,X_true) # ou bien Z
  colnames(surv_data) = c('Time','delta',paste("X", 1:p, sep = ""))
  
  L_naive =sapply(1:nA, FUN= function(i){  which.max(Q[i,])
  })  
  X_naive = X[L_naive,]
  data_naive = data.frame(Time,delta, X_naive) 
  
  XB = as.matrix(X, ncol = p) 
  return(list(surv_data=surv_data, data_naive=data_naive,XB=XB) )
}

######### Q matrix construction 
  C=c(0.7,0.2,0.1) 
  Q = matrix(0,nA,nB)
  for (i in 1:nA) {
    Q[i,i] = C[1]
    if(i!=1 & i!=nA){ 
      Q[i,i+1] = C[2]
      Q[i,i-1] = C[3]
    }else if(i==1)  {
      Q[i,i+1] = C[2]
      Q[i,i+2] = C[3]
    }else if(i==nA){ 
      Q[i,i-1] = C[2]
      Q[i,i-2] = C[3]}}
  
######################################


#Matrix_function<-function(gamma,XB,surv_data){

###### True-link data
#data_true = surv_data
# data_A = surv_data[,1:2] #Survival time and censored indicator

#n = nrow(data_A)
#m = nrow(XB)
#Q = matrix(0,n,m)
#L = matrix(0,n,m) #multinomial generation of the naive covariable
# Survival data (Database B)
#for (i in 1:n) {
# Q[i,] = rep((1-gamma)/m,m)
#Q[i,i] = gamma + (1-gamma)/m #Probability of True links
#}

#Lvec = sapply(1:n,FUN = function(i){
#  return(sample(1:m, size = 1, prob = Q[i,]))
#  })
# Lvec 

#L = matrix(0, nrow = n, ncol =m)
#L[cbind(1:n, Lvec)] = 1

#X_naive = L %*% XB

# Naive data
#data_naive = cbind(data_A, X_naive)

#return(list( data_true=data_true, data_A= data_A, data_naive=data_naive,Q=Q, L= L))
#}


