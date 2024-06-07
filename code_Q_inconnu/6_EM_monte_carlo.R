


library(survival)
library(rootSolve)
library(nleqslv)


setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_inconnu")
source("2_0_data_generate.R")
source("2_0_EM_matrice_Q.R")
source("2_EM_risk_function.R")
source("3_EM_naive_method.R")
source("4_EM_method_EM.R")
source("5_EM_scenarios.R")


################### nsim monte carlos
estimates_survival<- function(nsim,nA,nB,p,beta,K,prevalence){
  
  
  coef_true_s = matrix(0,nrow = nsim, ncol = p)
  
  coef_naive_s = matrix(0,nrow = nsim, ncol = p)
  converge_naive = vector()
  
  coef_estimate_s = matrix(0,nrow = nsim, ncol = p)
  converge_estimate = vector()
  
  
  for (i in 1:nsim){
    
    #data generation
    mf = Generate_data(K,nA,nB,prevalence,min_prev = 0.01)
    surv_data = mf$surv_data
    XB = mf$XB
    data_true=surv_data
    idBA=mf$idBA
    
    databaseA=mf$databaseA
    databaseB =mf$ databaseB
    
    datA=mf$datA # matching variables
    datB= mf$datB
    #####
    #record linkage 
    comp_mat = compare_binary (datA, datB, K)
    linked = Matrix_Q(comp_mat,K,datA, datB)
    Q=linked$Q
    apply(Q,1,sum)

    naive_id = linked$naive_id
    dif=linked$dif
    #####naive 
    Ndata = generate_naive_data(naive_id,XB , surv_data)
    naive_data = Ndata$data_naive
    Z = as.matrix(naive_data[,-c(1,2)])
    
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
    
    # EM method
    fit_estimate = Func_itteration(beta0,lambda0,Ts,event,XB,Q,tol= 1e-6,maxits = 300)
    coef_estimate_s[i,] = fit_estimate$beta0
    it_EM = fit_estimate$it
    
    
  }
  
  return(list( coef_true_s1 = coef_true_s[,1], coef_true_s2 = coef_true_s[,2], coef_true_s3 = coef_true_s[,3],
               coef_naive_s1 = coef_naive_s[,1], coef_naive_s2 = coef_naive_s[,2], coef_naive_s3 = coef_naive_s[,3],
               coef_estimate_s1 = coef_estimate_s[,1], coef_estimate_s2 = coef_estimate_s[,2],  coef_estimate_s3= coef_estimate_s[,3], 
               converge_estimate = converge_estimate))  } 
######################################



setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_inconnu")

for (i in (1:nrow(scenarios))){
  nsim=scenarios[i,1]
  K=scenarios[i,2]
  nA=scenarios[i,3]
  nB= scenarios[i,4]
  
  prevalence = rep(prevalence_sample, K/length(prevalence_sample))
  
  results_sample = estimates_survival(nsim,nA,nB,p,beta,K,prevalence) #monte carlos
  
  filename = paste0("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_inconnu/Results/","nsim=",nsim,
                    "_nA=",nA, "_K=",K,".Rdata")
  
  save(results_sample,file = filename)
  
}

