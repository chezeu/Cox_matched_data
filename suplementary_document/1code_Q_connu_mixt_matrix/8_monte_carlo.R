
library(survival)
library(rootSolve)
library(nleqslv)

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/suplementary_document/1code_Q_connu_mixt_matrix")
source("1_data_generate.R")
source("2_risk_function.R")
source("3_naive_method.R")
source("4_method_w_sum1.R")
source("5_method_w_sum2.R")
source("6_method_EM.R")
source("6_1_method_q_sum.R")
source("7_scenarios.R")


################### nsim monte carlos
estimates_survival<- function(nsim,nA,nB,p,beta,alpha,censor){
  
  
  coef_true_s = matrix(0,nrow = nsim, ncol = p)
  
  coef_naive_s = matrix(0,nrow = nsim, ncol = p)
  converge_naive = vector()
  
  coef_w_sum1_s = matrix(0,nrow = nsim, ncol = p)
  converge_w_sum1 = vector()
  
  coef_w_sum2_s = matrix(0,nrow = nsim, ncol = p)
  converge_w_sum2 = vector()
  
  coef_estimate_s = matrix(0,nrow = nsim, ncol = p)
  converge_estimate = vector()
  
  coef_Q_sum_s = matrix(0,nrow = nsim, ncol = p)
  converge_Q_sum = vector()
  
  # matrix 
  Q = matrix(0, nA, nB)
  
  t1= round((alpha*nA)) #individuals with (1,0,0)
  t2= nA - t1 #individuals with (0.6,0.2,0.2)
  for (i in 1:t1) {
    vec = sample(1:nB, 3)
    Q[i,vec] = c(1,0,0)
  }
  if(t1!=nA){ 
    for (i in (t1+1):nA ){
      vec = sample(1:nB, 3)
      Q[i,vec] = c(0.6,0.2,0.2)
    }
  }
  ########################  #number of individuals with percentage of true matches
  #  t1=round((60*nA)/100) # individuals with (0.8,0.1,0.1)
  # t2=round((20*nA)/100) #individuals with (0.6,0.2,0.2)
  # t3= nA -(t1+t2) #individuals with (0.4,0.3,0.3)
  
  # for (i in 1:t1) {
  #   vec = sample(1:nB, 3)
  #   Q[i,vec] = c(0.8,0.1,0.1)
  # }
  # for (i in (t1+1):(t1+t2) ){
  #   vec = sample(1:nB, 3)
  #   Q[i,vec] = c(0.6,0.2,0.2)
  # }
  # for (i in (t1+t2+1):(t1+t2+t3) ){
  #   vec = sample(1:nB, 3)
  #   Q[i,vec] = c(0.4,0.3,0.3)
  # }
  #############################  ##
  
  for (i in 1:nsim){
    #data generation
    mf = Generate_data(nB,nA,beta,Q,censor)
    surv_data = mf$surv_data
    XB = mf$XB
    data_naive = mf$data_naive
    data_true=surv_data
    
    Ts = data_true$Time
    event = data_true$delta
    
    Z = as.matrix(data_naive[,3:(p+2)]) 
    
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
    
    # method with Q average
    fit_Q_sum = coxph_Q_sum(Ts,event,XB, Q,maxiter = 20)
    coef_Q_sum_s[i,] = fit_Q_sum$coef
    converge_Q_sum[i] = fit_Q_sum$converge
    
    # EM method
    fit_estimate = Func_itteration(beta0,lambda0,Ts,event,XB,Q,tol= 1e-6,maxits = 500)
    coef_estimate_s[i,] = fit_estimate$beta0
    converge_estimate[i] = as.numeric(fit_estimate$converge)
    
  }
  
  return(list( coef_true_s1 = coef_true_s[,1], coef_true_s2 = coef_true_s[,2],coef_true_s3 = coef_true_s[,3],
               coef_naive_s1 = coef_naive_s[,1], coef_naive_s2 = coef_naive_s[,2], coef_naive_s3 = coef_naive_s[,3], 
               coef_w_sum1_s1 = coef_w_sum1_s[,1],coef_w_sum1_s2 = coef_w_sum1_s[,2],coef_w_sum1_s3 = coef_w_sum1_s[,3], 
               coef_w_sum2_s1 = coef_w_sum2_s[,1], coef_w_sum2_s2 = coef_w_sum2_s[,2],coef_w_sum2_s3 = coef_w_sum2_s[,3], 
               coef_Q_sum_s1 = coef_Q_sum_s[,1],coef_Q_sum_s2 = coef_Q_sum_s[,2], coef_Q_sum_s3= coef_Q_sum_s[,3], 
               coef_estimate_s1 = coef_estimate_s[,1], coef_estimate_s2 = coef_estimate_s[,2],  coef_estimate_s3 = coef_estimate_s[,3],
               converge_naive = converge_naive,  converge_w_sum1= converge_w_sum1,
               converge_w_sum2 = converge_w_sum2, converge_Q_sum= converge_Q_sum,converge_estimate = converge_estimate ))  } 
######################################



setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/suplementary_document/1code_Q_connu_mixt_matrix")

for (i in (1:nrow(scenarios))){
  nsim=scenarios[i,1]
  nB=scenarios[i,2]
  nA=scenarios[i,3]
  alpha = scenarios[i,4]
  censor=scenarios[i,5]
  
  results_sample = estimates_survival(nsim,nA,nB,p,beta,alpha,censor) #monte carlos
  
  filename = paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/suplementary_document/1code_Q_connu_mixt_matrix/Results/","nsim=",nsim,
                    "_nA=",nA,"_nB=",nB,"_alpha=",alpha,"_censor=",censor,".Rdata")
  
  save(results_sample,file = filename)
  
}

