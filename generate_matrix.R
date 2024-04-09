  

library(survival)
library(plyr)
library(dplyr)
library(rootSolve) 
library(nleqslv)

##### Generate la matrice de prob de couplage 
#<ek <- 0.03 # erreur
p2 <- 0.8# proba de la loi binomiale
m <- 1000 # taille de la base B
n <- 500# taille de la base A
N <- n*m # nombre de paire d'appariement
rho <- 0.75 # proportion of non censoring data
beta <- c (-0.5,0.5)  #parameter
p<-2 ## nombre de covariables
alpha<-0.8 #Probability of True links

## generer les bases de données
Generate_data <- function(p2,m,n,rho,beta){
 # set.seed(42)
  ##covarites data
  X1 <- rnorm(m,0,1)
  X2 <- rbinom(m,size = 1, prob = p2)
  X <- as.matrix(cbind(X1, X2))
  
  U <- runif(m,0,1) # Time to event variable
  Tt = -log(U)/exp((X%*%beta))  #lambda=1
  # Constant right censored
 # c = quantile(Tt[,1], probs = rho)
  
  Time <- pmin(Tt,c)
  delta <- as.numeric(Tt<=c)
  surv_data <- data.frame(Time,delta,X)
  #colnames(surv_data)=c('Time','delta',paste("X", 1:p, sep = ""))
  
  true_data <- surv_data[1:n,]
  
  database_A <- true_data[,1:2]  #base de donnee de A
  XB <- surv_data[,3:dim(surv_data)[2]] #base de donnee de covariables
  return(list(surv_data=surv_data,true_data = true_data, database_A  = database_A ,XB  =XB ))
}
#k <- dim(surv_data)[2]

#########  multinomiale
C<-c(0.7,0.3,1.5)## posterior probability values
Matrix_function<-function(n,m,C){
Q<-matrix(0,n,m)
Lvec<-matrix(0,n,m)
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
pro[i] = C[1] #Probability of True links
pro[i+1]= C[2]
Q[i,]<-pro
Lvec[i,] <- rmultinom(1, size = 1, pro = pro)} 
if(i==n){  pro = rep(0,m)
pro[i] = C[1] #Probability of True links
pro[i-1]=C[2]
Q[i,]<-pro
Lvec[i,] <- rmultinom(1, size = 1, pro= pro) }  
}
Lvec
X_naive <-Lvec %*%as.matrix( XB)
Z<- X_naive
colnames(Z)<- c("Z1","Z2")
# Naive data
data_naive <- cbind( database_A, Z)
return(list(Z=Z, data_naive=data_naive,Q=Q, Lvec= Lvec))
}


## equation naive

equa_naive <- function(beta, data_naive, XB, p) {
  m= nrow(XB)
  n = nrow(data_naive)
  
  Ts <-as.matrix( data_naive[,1]) 
  event <- data_naive[,2] 
  Z <- as.matrix(data_naive[,3:(p+2)])
  X <- as.matrix(XB[,1:p])
  
  ezbeta <- exp(Z%*% beta)
  zezbeta <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p) 

  s<-0
  for(i in 1:n){ 
   
    Yi<- as.numeric( (Ts == Ts[i] ) | (Ts > Ts[i]))# at risk 
    risq <- which(Yi==1)
    num_naive<-colSums( zezbeta[risq,])
    denum_naive<- sum(ezbeta[risq])    
    s<- s+ event[i]* ( Z[i,]- num_naive/denum_naive)
  }
  s
  return(H_naive=s) ##naive estimator 
}


############################## estimating equation

#equa_estimat <- function(beta, data_naive, XB, p,Q) {
#  m= nrow(XB)
#  n = nrow(data_naive)
  
#  Ts <-as.matrix( data_naive[,1]) 
#  event <- data_naive[,2] 
#  Z <- as.matrix(data_naive[,3:(p+2)])
#  X <- as.matrix(XB[,1:p])
  
# donnees <-Matrix_function(n,m,C)
#  Q <- donnees$Q
  
#  ezbeta <- exp(Z%*% beta)
#  eXbeta<-exp(X%*% beta)
  
#  zezbeta <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
 # XeXbeta<- X*matrix(rep(exp(X%*% beta),p), ncol = p)

    
    ## determiner tilde_f
  
 # t<-Q* t(matrix(rep(eXbeta,n),nrow  = m)) #  qexp(betax)
#  produit2<- vector() 
 
#for (j in 1:n) {
#  t_j<-t[j,]
#  L2 <- t_j[-j]
 # L1<-which( L2!=0) # l'ensemble pour les t_j diff de zero et diff de j
#produit2[j]<- prod(L2[L1])
#}  
#produit2  

#tilde_f<-matrix((1/diag(Q)))*  ezbeta  *  matrix(1/produit2 ) 

## determiner tilde_g

#produit<-vector()
#som<-matrix(0,nrow = n,ncol = p)
#for (k in 1:n) {
 # qk <- Q[k,]
#  t_k<-t[k,]
 # L2 <- t_k[-k]
 # L1<-which( L2!=0) # l'ensemble pour les t_k diff de zero 
#  produit[k]<- prod(t_k[L1])
#  som[k,]<-  colSums( (qk*X)[-k,]) * eXbeta[k]
  
#}  
#produit  
#som

#tilde_g <-   zezbeta* 1/ matrix( rep(diag(Q)*produit,p), ncol=p ) - som

## generale equation
  
 #     s<-0
  #  for(i in 1:n){ 
  #    Yi <- as.numeric( (Ts == Ts[i] ) | (Ts > Ts[i]))# at risk 
   #   risq <- which(Yi==1)
    #  qi <- Q[i,]
   
    # x_tilde
#    tilde_xi <- (qi[i]^-1)*Z[i,]- (qi[i]^-1)* colSums (qi[-i]*X[-i,])
    
 #   num <-colSums( tilde_g[risq,])
  #  denum <- sum(tilde_f[risq,])    
#    s<- s+ event[i]*( tilde_xi- num/denum )
#  }
#s
#return(H=s) ##naive estimator 
#}

#H<-equa_estimat(beta, data_naive, XB, p,Q)

################################# reecriture de l'equation estimante ###################

equa_estimat2 <- function(beta, data_naive, XB, p,Q) {
  m= nrow(XB)
  n = nrow(data_naive)
  
  Ts <-as.matrix( data_naive[,1]) 
  event <- data_naive[,2] 
  Z <- as.matrix(data_naive[,3:(p+2)])
  X <- as.matrix(XB[,1:p])
  
  ezbeta <- exp(Z%*% beta)
  eXbeta<-exp(X%*% beta)
  
  zezbeta <- Z*matrix(rep(exp(Z%*% beta),p), ncol = p)
  XeXbeta<- X*matrix(rep(exp(X%*% beta),p), ncol = p)
  
  donnees <-Matrix_function(n,m,C)
  Q <- donnees$Q
  

  #### autre ecriture de tilde_f
  produit3 <- vector()
  for (j in 1:n) {
    qj <- Q[j,]
    L1 <- which(qj!=0)
    qj <- qj[L1]
    eX <- eXbeta[-j]
    produit3[j] <- 1/(prod(qj)*prod(eX))
  }  
  produit3
  tilde_f2 <- ezbeta* produit3
  
  
  ### une autre ecriture de tilde_g
  
  produit5<-vector()
  som5<-matrix(0,nrow = n,ncol = p)
  for (k in 1:n) {
    qk <- Q[k,]
    L2 <- which(qk !=0)
    qk <- qk[L2]
    eX <- eXbeta[-k]
    produit5[k]<-1/ (prod(qk)*prod(eX))
    som5[k,]<-  colSums(X) * eXbeta[k]
    
  }  
  produit5  
  som5
  tilde_g2 <- zezbeta*produit5 - som5
  
  ######## determinons tilde_x
  
  tilde_x <- matrix(0,nrow=n, ncol = p)
  for(i in 1:n){ 
    qi <- Q[i,]
    qix <- qi*X
 tilde_x[i,] <- (1/qi[i])*Z[i,]  - colSums(qix)*(1/qi[i])
    }
  tilde_x
  

  ## reecriture de l'equation generale
  
  s2<-0
  for(i in 1:n){ 
    Yi <- as.numeric( (Ts == Ts[i] ) | (Ts > Ts[i]))# at risk 
    risq <- which(Yi==1)
      
    num <-colSums( tilde_g2[risq,])
    denum <- sum(tilde_f2[risq])    
    s2<- s2+ event[i]*( tilde_x[i,]- num/denum )
  }
  s2
  x1 <- s2[1]
  x2 <- s2[2]
  return(H=x1) ##naive estimator 
}
#H2<-equa_estimat(beta, data_naive, XB, p,Q)


####### naive equation solving

dat <- Generate_data(p2,m,n,rho,beta)
XB <- dat$XB
database_A <- dat$database_A
true_data <- dat$true_data
donnees <-Matrix_function(n,m,C)
data_naive <- donnees$data_naive

coxph_equa_naive <- function(data_naive, XB, p,maxiter = 20){
  f <- function(x){
    equa_naive (beta=x, data_naive=data_naive, XB = XB, p = p )
  }
  
  fit_manual <- multiroot(f,start = rep(0,p), maxiter = maxiter)
  beta <- fit_manual$root
  iterations <- fit_manual$iter
  converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
  
  return(list(coef = beta, converge = converge,iterations=iterations))
}

observed_naive <-coxph_equa_naive(data_naive, XB, p) 
coef_naive1 <- observed_naive$coef
conve_naive1 <- observed_naive$converge
iterations_naive1 <- observed_naive$iterations

### pour valider l'ecriture precedente 

#cof_cox <-coxph(Surv(Time,delta)~.,data = data_naive)
#cox_naive <- cof_cox$coefficients
#############estimate equation solving 

#dat <- Generate_data(p2,m,n,rho,beta)
#XB <- dat$XB
#database_A <- dat$database_A
#true_data <- dat$true_data
#donnees <-Matrix_function(n,m,C)
#data_naive <- donnees$data_naive
#Q <- donnees$Q

#coxph_equa_estimat <- function(data_naive, XB, p, Q,maxiter = 50){
 # f <- function(x){
#    equa_estimat (beta=x, data_naive=data_naive, XB = XB, p = p,Q=Q )
 # }
  
  #fit_manual <- multiroot(f,start = c(0,0.5))
#  g <- nleqslv( c(0,0),f, method="Newton")
  
  ## chercher la valeur initiale optimale
 # set.seed(14)
#  xstart <- matrix(runif(20,-1,1), ncol = 2)
 #Zero <-  searchZeros(xstart,f)
  
#  beta <- fit_manual$root
 # converge <- as.numeric((fit_manual$iter < maxiter)& !is.nan(fit_manual$estim.precis) & !is.na(fit_manual$estim.precis) & (fit_manual$estim.precis<1e-6))
#  iterations2 <- fit_manual$iter
 # maximum_f <- fit_manual$f.root
#  return(list(coef = beta, converge = converge,iterations2=iterations2, maximum_f=maximum_f))
#}

#observed2 <- coxph_equa_estimat(data_naive, XB, p, Q)
#coef2 <- observed2$coef
#conve2 <- observed2$converge
#iterations2 <- observed2$iterations2
#maximum_f <- observed2$maximum_f

############# rewritting the estimate equation solving 

dat <- Generate_data(p2,m,n,rho,beta)
XB <- dat$XB
database_A <- dat$database_A
true_data <- dat$true_data
donnees <-Matrix_function(n,m,C)
data_naive <- donnees$data_naive
Q <- donnees$Q
matrice_trans <- donnees$Lvec

coxph_equa_estimat2 <- function(data_naive, XB, p, Q,maxiter = 50){
  f2 <- function(x){
    equa_estimat2 (beta=x, data_naive=data_naive, XB = XB, p = p,Q=Q )
  }
  
  ## chercher la valeur initiale optimale
  
 # set.seed(14)
 # xstart <- matrix(runif(200,-1,1), ncol = 2)
 # Zero2 <-  searchZeros(xstart,f2)
  
 # if(!is.null(Zero2)){
  #xsol2 <-  matrix(Zero2$x, ncol = 2)
 # init2 <- Zero2$xstart
  #### recherche de la valeur f pour chaque beta
  #r <- nrow( xsol2)
  #imf2 <- matrix(0,nrow=r,ncol=p)
 # for (i in 1:r) {imf2[i,] <- f2( xsol2[1,])}
 # imf2
 # g2<- nleqslv( init2[1,] ,f2, method="Newton")
 # }else if(is.null(Zero2)){
 #   init2 <- rep(0,p)
 #   g2<- nleqslv( init2 ,f2, method="Newton")
 # }
 #g2 
  #fit_manual2 <- multiroot(f2, c(0,0))
 init2 <- rep(0,p)
 g2<- nleqslv( init2 ,f2, method="Newton")
  beta2 <-g2$x
  converge <- as.numeric((g2$iter < maxiter)& (g2$termcd==1) & !is.nan(g2$fvec))
  iterations2 <- g2$iter
  #image_f2 <-imf2
 # image_f2_used <- imf2[1,]
  return(list(coef = beta2, converge = converge,iterations=iterations2  ))
}

observed2 <- coxph_equa_estimat2(data_naive, XB, p, Q)
coef_estimat <- observed2$coef
conve2 <- observed2$converge
iterations2 <- observed2$iterations


############
# Theoretical estimating equation for true
fit_true <- coxph(Surv(Time,delta)~.,data = true_data)
coef_true <- as.vector(fit_true$coefficients)
var_true <- diag(fit_true$var)

################## monte carlo #####################
R <- 1000 
C <- c(0.8,0.2,0.1)

evalution <- function(C,m,n,R){

BETA_naive <-matrix(0,nrow = R,ncol = p)
BETA_estimat <-matrix(0,nrow = R,ncol = p)
conve2 <- matrix(0,nrow = R, ncol = p)
conve_naive1<- vector()

for (i in 1:R) {

p2 <- 0.8# proba de la loi binomiale
rho <- 0.75 # proportion of non censoring data
beta <- c ( 0.5,-0.5)  #parameter
p<-2 ## nombre de covariables

dat <- Generate_data(p2,m,n,rho,beta)
XB <- dat$XB
database_A <- dat$database_A
true_data <- dat$true_data
donnees <-Matrix_function(n,m,C)
data_naive <- donnees$data_naive
data_naive <- donnees$data_naive
Q <- donnees$Q
matrice_trans <- donnees$Lvec

observed_naive <-coxph_equa_naive(data_naive, XB, p) 
coef_naive1 <- observed_naive$coef
conve_naive1[i] <- observed_naive$converge

observed2 <- coxph_equa_estimat2(data_naive, XB, p, Q)
coef_estimat <- observed2$coef
conve2[i,] <- observed2$converge

BETA_naive[i,]  <- coef_naive1
BETA_estimat[i,]  <- coef_estimat

}
return(list(BETA_naive=BETA_naive, BETA_estimat=BETA_estimat, 
            conver_naive=conve_naive1,  conver_estimat= conve2))
}

result <- evalution(C,m,n,R)

  X1<- result$BETA_naive
X2<- result$BETA_estimat
conver_naive <- result$ conver_naive
conver_estimat <- result$ conver_estimat
base_result <- cbind(X1,X2,conver_naive ,conver_estimat)
colnames(base_result)[1:2] <- c("x1_naive","x2_naive")
colnames(base_result)[3:4] <- c("x1_estimat", "x2_estimat")
colnames(base_result)[6:7] <- c("conv_x1_estimat", "conv_x2_estimat")


base_result <- data.frame(base_result)
w <- which((base_result[,6] & base_result[,7])!=0)
base_result <- base_result[w,]

BOX_naive <- boxplot(base_result[,1:2])
BOX_estimat <-  boxplot(base_result[,3:4])
BETA_naive <- colMeans(base_result[,1:2])
BETA_estimat  <-  colMeans(base_result[,3:4])

fit_true <- coxph(Surv(Time,delta)~.,data = true_data)
BETA_true <- as.vector(fit_true$coefficients)

################################################################
################################ tracer la fonction estimante #######################

f2 <- function(x){
  equa_estimat2 (beta=x, data_naive=data_naive, XB = XB, p = p,Q=Q )
}
x1 <- seq(-0.2,0.2,by=0.1)
x2 <- seq(0.2,-0.2,by=-0.1)
x1lim <- range(x1)
x2lim <- range(x2)
x <- as.matrix(cbind(x1,x2),ncol=p)
a <- nrow(x)
y <- matrix(0,ncol=p,nrow =a )
for (i in 1:a) {
  y[i,] <- f2(x[i,])
}
y 
y2 <- matrix(0,ncol=p,nrow =a )
for (i in 1:a) {
  y2[i,] <- equa_naive(x[i,],data_naive = data_naive,XB,p)
}
y2 
Ylim <- range(y)
library(rgl)
library(plotly)
 
open3d()
plot3d(x1,x2,y,type="l", xlim=x1lim,ylim = x2lim,zlim=Ylim,box=TRUE)

plot3d(x1,x2,y2,type="l")

lmpoly <- lm( formula = y ~ I(x1^2) + I(x2^2) + I(x1*x2) + x1 +x2)
cat("equation de nappe: \n")
print(lmpoly)

x1.grid <- seq(x1lim[1],x1lim[2], length.out=45)
x2.grid <- seq(x2lim[1],x2lim[2], length.out=45)
x1_x2.grid <- expand.grid(x1=x1.grid,x2=x2.grid)

y.nappe <- with(x1_x2.grid,
                lmpoly$coefficients[2]*x1^2+lmpoly$coefficients[3]*x2^2+lmpoly$coefficients[4]*x1*x2+
                  lmpoly$coefficients[5]*x1+  lmpoly$coefficients[6]*x2+ lmpoly$coefficients[1])

y.mat <-matrix(y.nappe, nrow = 45, ncol = 45, byrow = T) 

##########################################

library(Ryacas)



f2 <- function(x){
  equa_estimat2 (beta=x, data_naive=data_naive, XB = XB, p = p,Q=Q )
}
x1 <- seq(-0.2,0.2,by=0.1)
x2 <- seq(0.2,-0.2,by=-0.1)
x1lim <- range(x1)
x2lim <- range(x2)
x <- as.matrix(cbind(x1,x2),ncol=p)
a <- nrow(x)
y <- vector()
for (i in 1:a) {
  y[i] <- f2(x[i,])
}
y 
range(y)
open3d()
plot3d(f2, col = colorRampPalette(c("blue", "white", "red")), 
       xlab = "X", ylab = "Y", zlab = "Sinc( r )", 
       xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2),
       aspect = c(-1.488456e+03,  6.072517e+72, 0.5))

plot3d(x1,x2,y,type="l", xlim=x1lim,ylim = x2lim,zlim=Ylim,box=TRUE)






