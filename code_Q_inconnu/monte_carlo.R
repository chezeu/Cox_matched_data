


library(survival)
library(rootSolve)
library(nleqslv)


setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")
source("1_data_generate.R")
source("2_matrix_generate.R")
source("3_naive_method.R")
source("4_method_w_sum1.R")
source("5_method_w_sum2.R")
source("6_method_EM.R")



################### monte carlos
estimates_survival<- function(nsim,n,m,C,p,beta,Q,lambda0,beta0){
  
  
  coef_true_s = matrix(0,nrow = nsim, ncol = p)
  
  coef_naive_s = matrix(0,nrow = nsim, ncol = p)
  converge_naive = vector()
  
  coef_w_sum1_s = matrix(0,nrow = nsim, ncol = p)
  converge_w_sum1 = vector()
  
  coef_w_sum2_s = matrix(0,nrow = nsim, ncol = p)
  converge_w_sum2 = vector()
  
  coef_estimate_s = matrix(0,nrow = nsim, ncol = p)
  converge_estimate = vector()
  
  
  for (i in 1:nsim){
    
    
    surv_data = Generate_data(m,n,beta)
    mf = Matrix_function(n,m,C,beta,surv_data)
    data_true = mf$data_true
    data_naive = mf$data_naive
    
    Ts = mf$data_true$Time
    event = mf$data_true$delta
    Q = mf$Q
    XB = mf$XB
    Z = as.matrix(data_naive[,3:(p+2)]) 
    
    lambda0 =  rep(0.1,length(event))
    beta0 = c(0.1,0.1)
    
    prob =  Functio_prob(beta0,lambda0,Ts,event,XB, Q)
    
    # Theoretical estimating equation for true and naive data
    fit_true = coxph(Surv(Time,delta)~.,data = data_true)
    coef_true_s[i,] = as.vector(fit_true$coefficients)
    
    fit_naive = coxph_equa_naive(Ts,event, Z, maxiter = 20)
    coef_naive_s[i,] = fit_naive$coef
    converge_naive[i] = fit_naive$converge
    
    fit_w_sum1 = coxph_w_sum1(Ts,event,XB, Q,maxiter = 20)
    coef_w_sum1_s[i,] = fit_w_sum1$coef
    converge_w_sum1[i] = fit_w_sum1$converge
    
    fit_w_sum2 = coxph_w_sum2(Ts,event,XB, Q, maxiter = 20)
    coef_w_sum2_s[i,] = fit_w_sum2$coef
    converge_w_sum2[i] = fit_w_sum2$converge
    
    fit_estimate = Func_itteration(beta,beta0,lambda0,Ts,event,XB, Q,tol= 1e-6,maxits = 500)
    coef_estimate_s[i,] = fit_estimate$beta0
    converge_estimate[i] = as.numeric(fit_estimate$converge)
    
    
  }
  
  return(list( coef_true_s1 = coef_true_s[,1], coef_true_s2 = coef_true_s[,2],
               coef_naive_s1 = coef_naive_s[,1], coef_naive_s2 = coef_naive_s[,2], coef_w_sum1_s1 = coef_w_sum1_s[,1], 
               coef_w_sum1_s2 = coef_w_sum1_s[,2], coef_w_sum2_s1 = coef_w_sum2_s[,1], coef_w_sum2_s2 = coef_w_sum2_s[,2],  
               coef_estimate_s1 = coef_estimate_s[,1], coef_estimate_s2 = coef_estimate_s[,2], converge_naive = converge_naive, 
               converge_w_sum1= converge_w_sum1,converge_w_sum2 = converge_w_sum2, converge_estimate = converge_estimate))  } 
######################################

