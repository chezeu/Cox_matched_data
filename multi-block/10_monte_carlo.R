


library(survival)
library(rootSolve)
library(nleqslv)


setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/multi-block")
source("02_data_generate.R")
source("03_block.R")
source("04_risk_function.R")
source("05_naive_method.R")
source("06_method_w_sum1.R")
source("07_method_w_sum2.R")
source("08_EM_method_EM.R")
source("09_scenarios.R")


################### nsim monte carlos
estimates_survival<- function(nsim,nA,nB,p,beta,K,prevalence){
  
  
  coef_true_s = matrix(0,nrow = nsim, ncol = p)
  
  coef_naive_s = matrix(0,nrow = nsim, ncol = p)
  converge_naive = vector()
  
  coef_w_sum1_s= matrix(0,nrow = nsim, ncol = p)
  converge_w_sum1 = vector()
  
  coef_w_sum2_s = matrix(0,nrow = nsim, ncol = p)
  converge_w_sum2 = vector()
  
  coef_estimate_s = matrix(0,nrow = nsim, ncol = p)
  converge_estimate = vector()
  
  
  for (i in 1:nsim){
    
    #data generation
    mf = multi_block ( K,nA_block, nB_block, p, prevalence, beta)
    surv_data = mf$surv_data
    XB = mf$XB[,-4]
    data_true=surv_data
    
   # databaseA=mf$databaseA
    #databaseB =mf$ databaseB
    
    datA=mf$datA # matching variables
    datB= mf$datB
    #####
    #record linkage 
    comp_mat = compare_binary (datA, datB, K, nA_block )
    linked = Matrix_Q(comp_mat,K,datA, datB,nA_block,nB_block)
    Q=linked$Q
    apply(Q,1,sum)

    naive_id = linked$naive_id
    #####naive 
    Ndata = generate_naive_data(naive_id,XB , surv_data)
    naive_data = Ndata$data_naive
    Z = as.matrix(naive_data[,-c(1,2,6)])
    
    Ts = data_true$Time
    event = data_true$delta
    
    lambda0 =  rep(0.1,length(event))
    lambda0 [which(event == 0)] = 0
    beta0 = rep(0.1,p)
    
    # times of event
    # death_times= Ts[which(event==1)]
    
    
    # Theoretical estimating equation for true and naive data
    fit_true = coxph(Surv(Time,delta)~.,data = data_true)
    coef_true_s[i,] = as.vector(fit_true$coefficients)
    
    # naive method
    fit_naive = coxph_equa_naive(Ts,event, Z, maxiter = 20)
    coef_naive_s[i,] = fit_naive$coef
    converge_naive[i] = fit_naive$converge
    
    # method with weighted average
    fit_w_sum1 = coxph_w_sum1(Ts,event,XB, Q,maxiter = 20)
    coef_w_sum1_s[i,] = fit_w_sum1$coef
    converge_w_sum1[i] = fit_w_sum1$converge
    
    #method with the maximum of proba
    fit_w_sum2 = coxph_w_sum2(Ts,event,XB, Q, maxiter = 20)
    coef_w_sum2_s[i,] = fit_w_sum2$coef
    converge_w_sum2[i] = fit_w_sum2$converge
    
    # EM method
    fit_estimate = Func_itteration(beta0,lambda0,Ts,event,XB,Q,tol= 1e-6,maxits = 300)
    coef_estimate_s[i,] = fit_estimate$beta0
    it_EM = fit_estimate$it
    
    
  }
  
  return(list( coef_true_s1 = coef_true_s[,1], coef_true_s2 = coef_true_s[,2],coef_true_s3 = coef_true_s[,3],
               coef_naive_s1 = coef_naive_s[,1], coef_naive_s2 = coef_naive_s[,2], coef_naive_s3 = coef_naive_s[,3], coef_w_sum1_s1 = coef_w_sum1_s[,1], 
               coef_w_sum1_s2 = coef_w_sum1_s[,2],coef_w_sum1_s3 = coef_w_sum1_s[,3], coef_w_sum2_s1 = coef_w_sum2_s[,1], coef_w_sum2_s2 = coef_w_sum2_s[,2],  
               coef_w_sum2_s3 = coef_w_sum2_s[,3], coef_estimate_s1 = coef_estimate_s[,1], coef_estimate_s2 = coef_estimate_s[,2], 
               coef_estimate_s3 = coef_estimate_s[,3],  converge_naive = converge_naive, 
               converge_w_sum1= converge_w_sum1,converge_w_sum2 = converge_w_sum2, converge_estimate = converge_estimate )) 
  } 
######################################



setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/multi-block")

for (i in (1:nrow(scenarios))){
  nsim=scenarios[i,1]
  K=scenarios[i,2]
  nA=scenarios[i,3]
  nB= scenarios[i,4]
  nA_block= c(scenarios[i,5],scenarios[i,6],scenarios[i,7] )
  nB_block= c(scenarios[i,8],scenarios[i,9],scenarios[i,10] )
  prevalence = rep(prevalence_sample, K/length(prevalence_sample))
  
  results_sample = estimates_survival(nsim,nA,nB,p,beta,K,prevalence) #monte carlos
  
  filename = paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/multi-block/Results/","nsim=",nsim,
                    "_nA=",nA, "_K=",K,".Rdata")
  
  save(results_sample,file = filename)
}

