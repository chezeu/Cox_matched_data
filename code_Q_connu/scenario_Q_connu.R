
sigma=1
alpha=1
C_sample=cbind(c(0.6,0.2,0.2),c(0.7,0.15,0.15),c(0.8,0.1,0.1))
m=70
n=40
p=2
nsim=10

beta=c(0.5,-0.5)
lambda0<- rep(0.5,n)
beta0<- c(0.1,0.1) 

 # tableau de parametres

scenarios=NULL
for (i in 1:ncol(C_sample)){
  scenarios=rbind(scenarios,c(nsim,n,m,p,sigma,alpha,C_sample[,i]))
}

colnames(scenarios)=c("nsim","n","m", "p","sigma","alpha","Prob_1","Prob_2","Prob_3")

scenarios=data.frame(scenarios)

###################

estimates_survival<- function(nsim,n,m,p,sigma,alpha,C,beta, lambda0, beta0){
  
  coef_true_s<-matrix(0,nrow = nsim, ncol = p)
  
  coef_naive_s<- matrix(0,nrow = nsim, ncol = p)
  
  coef_w_sum1_s<- matrix(0,nrow = nsim, ncol = p)
  
  coef_w_sum2_s<- matrix(0,nrow = nsim, ncol = p)
  
  coef_estimate_s <- matrix(0,nrow = nsim, ncol = p)
  
  converge_estimate_s<- vector()
  
  for (i in 1:nsim){
    
    surv_data =Generate_data(m,n,beta)
       
      g_naive <- doOne2d_equat_estimat(n,m,C,beta,surv_data)
      coef_naive_s[i,] <- g_naive$coef_naive2
      
     g_sum1 <- doOne2d_w_sum1 (n,m,surv_data,C)
      coef_w_sum1_s[i,] <- g_sum1$coef_w_sum1
      
      g_sum2 <- doOne2d_w_sum2 (n,m,surv_data,C)
      coef_w_sum2_s[i,] <- g_sum2$coef_w_sum2
     
      g_estim <- doOne2d_estimate (n,m,C,beta0,lambda0,surv_data)
      coef_estimate_s [i,] <- g_estim$coef_estimate
      coef_true_s[i,] <- g_estim$coef_true
      converge_estimate_s<-g_estim$converge_estimate 
    }
  
    
  return(list( coef_true_s1=coef_true_s[,1],coef_true_s2=coef_true_s[,2],coef_naive_s1=coef_naive_s[,1],
               coef_naive_s2=coef_naive_s[,2], coef_w_sum1_s1=coef_w_sum1_s[,1],coef_w_sum1_s2=coef_w_sum1_s[,2],
               coef_w_sum2_s1=coef_w_sum2_s[,1], coef_w_sum2_s2=coef_w_sum2_s[,2],  
               coef_estimate_s1=coef_estimate_s[,1],coef_estimate_s2=coef_estimate_s[,2]))

} 

######################################

set.seed(23)

for (i in (1:nrow(scenarios))){
  
  nsim_sample=scenarios[i,1]
  n_sample=scenarios[i,2]
  m_sample=scenarios[i,3]
  p_sample=scenarios[i,4]
  sigma_sample=scenarios[i,5]
  alpha_sample=scenarios[i,6]
  C=vector()
  C[1]=scenarios[i,7]
  C[2]=scenarios[i,8]
  C[3]=scenarios[i,9]
  
  results_sample<-estimates_survival (nsim=nsim_sample,n=n_sample,m=m_sample,
                      p=p_sample,sigma=sigma_sample,alpha=alpha_sample,
                     C=C, beta, lambda0, beta0)
  
  boxplot(results_sample)
 
}

