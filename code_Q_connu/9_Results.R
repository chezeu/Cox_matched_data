
############results


setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")
source("8_scenarios.R")

# list of files to load
#library(parallel) 

file_list <- c("scenario=1_m=1000_prob1=0.6.Rdata", 
               "scenario=2_m=1000_prob1=0.7.Rdata",  
               "scenario=3_m=1000_prob1=0.8.Rdata",
               "scenario=4_m=500_prob1=0.8.Rdata",
               "scenario=5_m=1000_prob1=0.8.Rdata",
               "scenario=6_m=2000_prob1=0.8.Rdata") 


r = nrow(scenarios) 

funct_tableau <- function(file_list,r,p){
  
  beta_true = matrix(0,r,p )
  beta_naive = matrix(0,r,p )
  beta_w_sum1 = matrix(0,r,p )
  beta_w_sum2 = matrix(0,r,p )
  beta_estimate = matrix(0,r,p )
  
  L_naive = vector()  # number of convergences
  L_w_sum1 = vector()
  L_w_sum2 = vector()
  L_estimate = vector()
  
  Sd_true = matrix(0,r,p ) # standard deviation
  Sd_naive = matrix(0,r,p )
  Sd_sum1 = matrix(0,r,p )
  Sd_sum2 = matrix(0,r,p )
  Sd_estimate = matrix(0,r,p )
  
  rmse_true = matrix(0,r,p ) # mean squart error
  rmse_naive = matrix(0,r,p )
  rmse_sum1 = matrix(0,r,p )
  rmse_sum2 = matrix(0,r,p )
  rmse_estimate = matrix(0,r,p )
  
  for (i in 1:length(file_list)) {
    
    setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu/Results")
    load(file_list[i])
    results_sample = results_sample
    
    beta_true[i,] = c (mean(results_sample$coef_true_s1), mean(results_sample$coef_true_s2))
    beta_naive[i,] = rbind(  mean(results_sample$coef_naive_s1), mean(results_sample$coef_naive_s2))
    beta_w_sum1[i,] = rbind( mean(results_sample$coef_w_sum1_s1), mean(results_sample$coef_w_sum1_s2))
    beta_w_sum2[i,] = rbind( mean(results_sample$coef_w_sum2_s1), mean(results_sample$coef_w_sum2_s2))
    beta_estimate[i,] = rbind( mean(results_sample$coef_estimate_s1 ), mean(results_sample$coef_estimate_s2 ))
    
    
    L_naive[i] = length(which(results_sample$converge_naive==0))
    L_w_sum1[i] = length(which(results_sample$converge_w_sum1==0))
    L_w_sum2[i] = length( which(results_sample$converge_w_sum2==0))
    L_estimate[i] = length(which(results_sample$converge_estimate==0))
    
    
    Sd_true[i,] = rbind( sd(results_sample$coef_true_s1), sd(results_sample$coef_true_s2))
    Sd_naive[i,] = rbind( sd( results_sample$coef_naive_s1) , sd( results_sample$coef_naive_s2) )
    Sd_sum1[i,] = rbind( sd(results_sample$coef_w_sum1_s1), sd(results_sample$coef_w_sum1_s2))
    Sd_sum2 [i,] = rbind( sd(results_sample$coef_w_sum2_s1),sd(results_sample$coef_w_sum2_s2))
    Sd_estimate[i,] = rbind( sd(results_sample$coef_estimate_s1),sd(results_sample$coef_estimate_s2))
    
    
    rmse_true[i,] = rbind (sqrt( mean((results_sample$coef_true_s1- beta[1])^2)),
                           sqrt( mean((results_sample$coef_true_s2+  beta[2])^2)))
    
    rmse_naive[i,] = rbind( sqrt( mean((results_sample$coef_naive_s1- beta[1])^2)),
                            sqrt( mean((results_sample$coef_naive_s2+  beta[2])^2)))
    
    rmse_sum1[i,] = rbind( sqrt( mean((results_sample$coef_w_sum1_s1- beta[1])^2)),
                           sqrt( mean((results_sample$coef_w_sum1_s2- beta[2])^2)))
    
    rmse_sum2[i,] = rbind( sqrt( mean((results_sample$coef_w_sum2_s1- beta[1])^2)),
                           sqrt( mean((results_sample$coef_w_sum2_s2- beta[2])^2)))
    
    rmse_estimate[i,] = rbind( sqrt( mean((results_sample$coef_estimate_s1 -  beta[1])^2)),
                               sqrt( mean((results_sample$coef_estimate_s2 -  beta[2])^2)))
    
    
  }
  
  return( list( beta_true = beta_true, beta_naive = beta_naive,beta_w_sum1=beta_w_sum1,
          beta_w_sum2=beta_w_sum2, beta_estimate=beta_estimate,
        L_naive =L_naive,L_w_sum1=L_w_sum1,L_w_sum2=L_w_sum2,L_estimate=L_estimate,
    Sd_true= Sd_true,Sd_naive=Sd_naive ,Sd_sum1=Sd_sum1,  Sd_sum2=Sd_sum2, Sd_estimate =Sd_estimate,
   rmse_true = rmse_true, rmse_naive = rmse_naive,rmse_sum1 = rmse_sum1, rmse_sum2 = rmse_sum2, rmse_estimate = rmse_estimate))
  
}
#####################

Results_scenario = as.data.frame(funct_tableau(file_list,r,p))
rownames(Results_scenario )= paste0( "Scenario", 1:nrow(Results_scenario))
