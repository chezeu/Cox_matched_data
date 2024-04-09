##### Generate la matrice de prob de couplage 
ek <- 0.02 # erreur
p1 <- 0.8 # proba de la loi binomiale
p2 <- 0.7
n <- 5 # taille de la base B
l <- 3 # taille de la base A
N <- n*l # nombre de paire d'appariement
rho <- 0.75 # proportion of non censoring data
beta <- c ( 0.25,-0.25)  #estimate time
p<-2 ## nombre de covariables
alpha<-0.8 #Probability of True links


Generate_data <- function(p1,p2,n,l,rho,beta){
  set.seed(42)
  ##covarites data
  X1 <- rbinom(n,size = 1, prob = p1)
  X2 <- rbinom(n,size = 1, prob = p2)
  X <- as.matrix(cbind(X1, X2))
  
  U <- runif(n,0,1) # Time to event variable
  Tt = -log(U)/exp((X%*%beta))  #lambda=1
  # Constant right censored
  c = quantile(Tt[,1], probs = rho)
  
  Time <- pmin(Tt,c)
  delta <- as.numeric(Tt<=c)
  surv_data <- data.frame(Time,delta,X)
  #colnames(surv_data)=c('Time','delta',paste("X", 1:p, sep = ""))
  
  n <- dim(surv_data)[1]
  d2 <- sample(1:n,l)
  true_data <- surv_data[d2,]
  
  database_A <- true_data[,1:2]  #base de donnee de A
  XB <- surv_data[,3:dim(surv_data)[2]] #base de donnee de covariables
  return(list(surv_data=surv_data,true_data = true_data, database_A  = database_A ,XB  =XB ,d2=d2))
}
#k <- dim(surv_data)[2]


generate_linked_mixed <- function(n, N, p, beta, rho , alpha){
  
  # Survival data (Database B)
 data <- Generate_data(p1,p2,n,l,rho,beta)
  X <- data$XB
  XB <- as.matrix(X, ncol = p) #Only covariates
  ind_A<-data$d2 ##indices selectionés au hazard
  ###### True-link data: The first n individual from B
  surv_data<-data$surv_data
  data_true <- data$true_data
  
  
  ###### Linked data
  data_A <- data$database_A   #Survival time and censored indicator
  Q<-matrix(0,l,n)
  Lvec<-vector()
  for (i in 1:l) {
    prob = rep((1-alpha)/(n-1),n)
  prob[i] = alpha #Probability of True links
  Q[i,]<-prob
  Lvec[i] <- sample(1:n, size = 1, prob = prob)
  }
  Lvec 
  Q 
  L = matrix(0, nrow = l, ncol =n)
  L[cbind(1:l, Lvec)] = 1
  
  X_naive <- L %*% XB
  
  # Naive data
  data_naive <- cbind(data_A, X_naive)
  trueLink_id <- diag(L)
  
  return(list(surv_data=surv_data,data_true = data_true, data_A = data_A,   XB = XB, data_naive = data_naive,
              trueLink_id = trueLink_id,Q=Q))
}



#1. Finding the risk set R(t) given some time t
GetRiskSet <- function(time_of_interest, time_vector, event_vector) {
  
  return(which(((time_vector == time_of_interest & event_vector == 1) | (time_vector > time_of_interest))))
  
}

#2. Proposed estimating equations 

linked<- generate_linked_mixed (n, N, p, beta, rho , alpha)
n = nrow(XB)
l = nrow(data_naive)
alpha<-alpha
Ts <- data_naive[,1] 
event <- data_naive[,2] 
Z <- as.matrix(data_naive[,3:(p+2)])
X <- as.matrix(XB[,1:p])
Q<-linked$Q

 tildeX<-matrix(0,l ,ncol(XB))
  colnames(tildeX)<-c("tildex1", "tildex2")
  
  for (i in 1:l){
    S<-0
    qi<-Q[i,]
    for (j in 1:n) {
      if(qi[j]!=alpha) {qij<- qi[j] ;S<-S+qij*XB[j,] }   
    }
    S
    df <-1/qi[i]*Z[i,]-S*1/qi[i]
    
    tildeX[i,]<- df}
  tildeX
  
  dat <- cbind(Ts,  tildeX)[which(event==1),]
  
  
  
  
  ###### Numerateur
  
  tilde_g<-matrix(0,l,ncol(XB))
  colnames(tilde_g)<-c("tilde_g1", "tilde_g2")
 
  for (i in 1:dim(New_data)[1] ){
    S<-0
    b<-1
    k<-d2[i]
    for (j in 1:dim(XB)[1]) {
      qij<- (1-qi[i])/ (dim(XB)[1]-1) ;Xj<- as.numeric(XB[j,]);  b<-b*qij*exp(beta %*% Xj) 
      b
      if(j !=k) { Xj<- as.numeric(XB[j,]);  S<-S+XB[j,]*exp(beta %*% Xj) }   
      S
    }
    Zi<- as.numeric(Z[i,])
    df<- Z[i,]*exp(beta%*%Zi)/b - S
    tilde_g[i,]<- df
  }
  tilde_g 
  
  
  
  
