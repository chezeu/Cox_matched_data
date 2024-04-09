
beta=c(0.5,-0.5)
sigma=1
alpha=1
C=c(0.6,0.2,0.2)
m=15
n=10
p=2
nsim=10
lambda0<- rep(0.5,n)
beta0<- c(0.1,0.1) 

funct_test <- function(m,n,beta,beta0,p,lambda0,surv_data,nsim, C){
#applications true data

coef_true<-matrix(0,nrow = nsim, ncol = p)

#applications naive data

coef_naive<- matrix(0,nrow = nsim, ncol = p)
converge_naive <- vector()

#applications weigted sum approach 1

coef_w_sum1<- matrix(0,nrow = nsim, ncol = p)
converge_w_sum1 <- vector()

#applications weigted sum approach 2

coef_w_sum2<- matrix(0,nrow = nsim, ncol = p)
converge_w_sum2 <- vector()

#applications estimate

coef_estimate <- matrix(0,nrow = nsim, ncol = p)
converge_estimate<- vector()
lambda0_estimate<-matrix(0,nrow = nsim, ncol = n)


  for (i in 1:nsim) {
    
#generate data    
    surv_data <- Generate_data(m,n,beta)
    
#naive    
    g_naive <- doOne2d_equat_estimat(n,m,C,beta,surv_data)
    
    coef_naive[i,] <- g_naive$coef_naive2
    converge_naive[i] <- g_naive$converge_naive2

# sum 1   
    g_sum1 <- doOne2d_w_sum1 (n,m,surv_data,C)
    
    coef_w_sum1[i,] <- g_sum1$coef_w_sum1
    converge_w_sum1[i] <- g_sum1$converge_w_sum1 
    
 #sum 2   
    g_sum2 <- doOne2d_w_sum2 (n,m,surv_data,C)
    
    coef_w_sum2[i,] <- g_sum2$coef_w_sum2
    converge_w_sum2[i] <- g_sum2$converge_w_sum2
    
 #estimate
    g_estim <- doOne2d_estimate (n,m,C,beta0,lambda0,surv_data)
    
    coef_estimate [i,] <- g_estim$coef_estimate
    converge_estimate [i] <- g_estim$converge_estimate 
    lambda0_estimate[i,]<-g_estim$lambda0_estimate
  
#true    
    coef_true[i,] <- g_estim$coef_true
    
  }
  
  L_naive <- which(converge_naive!=0)
  coef_naive<- coef_naive[L_naive,]
  
  L_sum1 <- which(converge_w_sum1!=0)
  coef_w_sum1<- coef_w_sum1[L_sum1,]
  
  L_sum2 <- which(converge_w_sum2!=0)
  coef_w_sum2<- coef_w_sum2[L_sum2,]
  
  L_estim <- which(converge_estimate !=0)
  coef_estimate<- coef_estimate [L_estim,]
  
  return(list(coef_naive=coef_naive,converge_naive=converge_naive,
              coef_w_sum1= coef_w_sum1, converge_w_sum1=converge_w_sum1,
              coef_w_sum2= coef_w_sum2, converge_w_sum2=converge_w_sum2,
              coef_estimate = coef_estimate , converge_estimate= converge_estimate,
              coef_true= coef_true,lambda0_estimate=lambda0_estimate))
}




#################################

#results

result_test <- funct_test(m,n,beta,beta0,p,lambda0,surv_data,nsim, C)

coef_naive <-colMeans(result_test$coef_naive)

coef_w_sum1<-colMeans(result_test$coef_w_sum1)

coef_true<-colMeans(result_test$coef_true)

coef_w_sum2<-colMeans(result_test$coef_w_sum2)

coef_estimate<-colMeans(result_test$coef_estimate )

lambda0_estimate<-result_test$lambda0_estimate

####################
