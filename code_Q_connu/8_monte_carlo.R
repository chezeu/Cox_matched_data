


library(survival)
library(rootSolve)
library(nleqslv)


setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")
source("1_data_generate.R")
source("2_risk_function.R")
source("3_naive_method.R")
source("4_method_w_sum1.R")
source("5_method_w_sum2.R")
source("6_method_EM.R")
source("EM_with_N_d.R")
source("7_scenarios.R")


################### nsim monte carlos
estimates_survival<- function(nsim,nA,nB,C,p,beta,sigma,alpha){
  
  
  coef_true_s = matrix(0,nrow = nsim, ncol = p)
  
  coef_naive_s = matrix(0,nrow = nsim, ncol = p)
  converge_naive = vector()
  
  coef_w_sum1_s = matrix(0,nrow = nsim, ncol = p)
  converge_w_sum1 = vector()
   
  coef_w_sum2_s = matrix(0,nrow = nsim, ncol = p)
  converge_w_sum2 = vector()
  
  coef_estimate_s = matrix(0,nrow = nsim, ncol = p)
  converge_estimate = vector()
  
  coef_estimate_s_Nd = matrix(0,nrow = nsim, ncol = p)
  converge_estimate_Nd = vector()
  
  #Q = matrix(0,nA,nB)
  #for (i in 1:nA) {
   # Q[i,i] = C[1]
  #  if(i!=1 & i!=nA){ 
   #   Q[i,i+1] = C[2]
    #  Q[i,i-1] = C[3]
  #  }else if(i==1)  {
   #   Q[i,i+1] = C[2]
    #  Q[i,i+2] = C[3]
    #}else if(i==nA){ 
    #  Q[i,i-1] = C[2]
     # Q[i,i-2] = C[3]}}
  
  #C=c(0.8,0.1,0.1)
  
  Q = matrix(0, nA, nB)
  for (j in 1:nA) {
   vec = sample(1:nB, 3)
    Q[j,vec] = C
  }
  
  for (i in 1:nsim){
    
    #data generation
    mf = Generate_data(nB,nA,beta,Q)
    surv_data = mf$surv_data
    XB = mf$XB
    data_naive = mf$data_naive
    data_true=surv_data
    
    Ts = data_true$Time
    event = data_true$delta
  
    Z = as.matrix(data_naive[,-c(1,2)]) 

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
    fit_estimate = Func_itteration(beta0,lambda0,Ts,event,XB,Q,tol= 1e-6,maxits = 100)
    coef_estimate_s[i,] = fit_estimate$beta0
    converge_estimate[i] = as.numeric(fit_estimate$converge)
    
    # EM with N_d
    fit_estimate_Nd = Func_itteration_Nd(beta0,lambda0,Ts,event,XB, Q,tol= 1e-6,maxits = 100)
    coef_estimate_s_Nd[i,] = fit_estimate_Nd$beta0
    converge_estimate_Nd[i] = as.numeric(fit_estimate_Nd$converge)
    
    
  }
  
  return(list( coef_true_s1 = coef_true_s[,1], coef_true_s2 = coef_true_s[,2],
               coef_naive_s1 = coef_naive_s[,1], coef_naive_s2 = coef_naive_s[,2], coef_w_sum1_s1 = coef_w_sum1_s[,1], 
               coef_w_sum1_s2 = coef_w_sum1_s[,2], coef_w_sum2_s1 = coef_w_sum2_s[,1], coef_w_sum2_s2 = coef_w_sum2_s[,2],  
               coef_estimate_s1 = coef_estimate_s[,1], coef_estimate_s2 = coef_estimate_s[,2], 
               coef_estimate_s_Nd1 = coef_estimate_s_Nd[,1], coef_estimate_s_Nd2 = coef_estimate_s_Nd[,2],converge_naive = converge_naive, 
               converge_w_sum1= converge_w_sum1,converge_w_sum2 = converge_w_sum2, converge_estimate = converge_estimate,
               converge_estimate_Nd = converge_estimate_Nd ))  } 
######################################



setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")

for (i in (1:nrow(scenarios))){
  nsim=scenarios[i,1]
  nB=scenarios[i,2]
  nA=scenarios[i,3]
  C=vector()
  C[1]=scenarios[i,4]
  C[2]=scenarios[i,5]
  C[3]=scenarios[i,6]
#gamma=scenarios[i,4]
  
    results_sample = estimates_survival(nsim,nA,nB,C,p,beta,sigma,alpha) #monte carlos
  
  filename = paste0("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu/Results/","nsim=",nsim,
                    "_nA=",nA, "_prob1=", C[1], ".Rdata")
  
  save(results_sample,file = filename)
  
}

