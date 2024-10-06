


library(survival)
library(rootSolve)
library(nleqslv)

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/FS-method")
source("7_scenarios.R")
source("1_FS.R")
source("1_data_generate.R")
source("2_risk_function.R")
source("3_naive_method.R")
source("4_method_w_sum1.R")
source("5_method_w_sum2.R")
source("6_method_EM.R")





################### nsim monte carlos
estimates_survival<- function(nsim,nA,nB,p,beta,datasetA,datasetB,var_block,matching_variables){
  
  coef_true_s = matrix(0,nrow = nsim, ncol = p)
  
  coef_naive_s = matrix(0,nrow = nsim, ncol = p)
  converge_naive = vector()
  
  coef_w_sum1_s = matrix(0,nrow = nsim, ncol = p)
  converge_w_sum1 = vector()
  
  coef_w_sum2_s = matrix(0,nrow = nsim, ncol = p)
  converge_w_sum2 = vector()
  
  coef_estimate_s = matrix(0,nrow = nsim, ncol = p)
  converge_estimate = vector()
  
  prob_match = linkage_function (datasetA,datasetB,var_block,matching_variables)
  data_block = function_block ( prob_match, datasetA, datasetB)
  blockA =data_block$blockA
  new_dataA= data_block$new_dataA
  new_dataB = data_block$new_dataB
  
  Q= Matrix_block (blockA,datasetA,datasetB,new_dataA,new_dataB,prob_match)
  
  for (i in 1:nsim){
 
    #data generation
    mf = Generate_data(nA,nB,idBA,beta)
    surv_data = mf$surv_data
    XB = mf$XB
    data_true=surv_data
    
    Ts = data_true$Time
    event = data_true$delta
    
    naive_id = apply( Q,1, which.max)
    Z = as.matrix(XB[naive_id,])
    naive_data = cbind(Ts,event,Z)
    naive_data= data.frame(naive_data)
   # Matching_naive= new_dataB[naive_id,]
    
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
    fit_estimate = Func_itteration(beta0,lambda0,Ts,event,XB,Q,tol= 1e-6,maxits = 500)
    coef_estimate_s[i,] = fit_estimate$beta0
    converge_estimate[i] = as.numeric(fit_estimate$converge)
    
  }
  
  return(list( coef_true_s1 = coef_true_s[,1], coef_true_s2 = coef_true_s[,2],coef_true_s3 = coef_true_s[,3],
               coef_naive_s1 = coef_naive_s[,1], coef_naive_s2 = coef_naive_s[,2], coef_naive_s3 = coef_naive_s[,3], coef_w_sum1_s1 = coef_w_sum1_s[,1], 
               coef_w_sum1_s2 = coef_w_sum1_s[,2],coef_w_sum1_s3 = coef_w_sum1_s[,3], coef_w_sum2_s1 = coef_w_sum2_s[,1], coef_w_sum2_s2 = coef_w_sum2_s[,2],  
               coef_w_sum2_s3 = coef_w_sum2_s[,3], coef_estimate_s1 = coef_estimate_s[,1], coef_estimate_s2 = coef_estimate_s[,2], 
               coef_estimate_s3 = coef_estimate_s[,3],  converge_naive = converge_naive, 
               converge_w_sum1= converge_w_sum1,converge_w_sum2 = converge_w_sum2, converge_estimate = converge_estimate ))  } 
######################################



setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/FS-method")
  
  results = estimates_survival(nsim,nA,nB,p,beta,datasetA,datasetB,var_block,matching_variables) #monte carlos
  
  filename = paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/FS-method/Results/","nsim=",nsim,
                    "_nA=",nA,".Rdata")
  
  save(results,file = filename)
  

