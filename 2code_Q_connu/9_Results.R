
############ results


setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")

library(xtable)
source("7_scenarios.R")
r = nrow(scenarios) # number of scenarios

# parameter beta  
beta_true = matrix(0,r,p ) 
beta_naive = matrix(0,r,p )
beta_w_sum1 = matrix(0,r,p )
beta_w_sum2 = matrix(0,r,p )
beta_estimate = matrix(0,r,p )

# number of convergences
L_naive = vector()  
L_w_sum1 = vector()
L_w_sum2 = vector()
L_estimate = vector()

# standard deviation of the parameters  
Sd_true = matrix(0,r,p ) 
Sd_naive = matrix(0,r,p )
Sd_sum1 = matrix(0,r,p )
Sd_sum2 = matrix(0,r,p )
Sd_estimate = matrix(0,r,p )

# mean square error
rmse_true = matrix(0,r,p )
rmse_naive = matrix(0,r,p )
rmse_sum1 = matrix(0,r,p )
rmse_sum2 = matrix(0,r,p )
rmse_estimate = matrix(0,r,p )

# bias 
bias_true = matrix(0,r,p )
bias_naive = matrix(0,r,p )
bias_w_sum1 = matrix(0,r,p )
bias_w_sum2 = matrix(0,r,p )
bias_estimate = matrix(0,r,p )

for (i in 1:nrow(scenarios)) {
  
  nsim=scenarios[i,1]
  m=scenarios[i,2]
  n=scenarios[i,3]
  C=vector()
  C[1]=scenarios[i,4]
  C[2]=scenarios[i,5]
  C[3]=scenarios[i,6]
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu/Results/",
              "nsim=",nsim,"_nA=",nA, "_prob1=", C[1], ".Rdata"))
  results_sample = results_sample
  
  # beta
  beta_true[i,] = c (mean(results_sample$coef_true_s1), mean(results_sample$coef_true_s2),mean(results_sample$coef_true_s3))
  beta_naive[i,] = rbind(  mean(results_sample$coef_naive_s1), mean(results_sample$coef_naive_s2), mean(results_sample$coef_naive_s3))
  beta_w_sum1[i,] = rbind( mean(results_sample$coef_w_sum1_s1), mean(results_sample$coef_w_sum1_s2), mean(results_sample$coef_w_sum1_s3))
  beta_w_sum2[i,] = rbind( mean(results_sample$coef_w_sum2_s1), mean(results_sample$coef_w_sum2_s2), mean(results_sample$coef_w_sum2_s3))
  beta_estimate[i,] = rbind( mean(results_sample$coef_estimate_s1 ), mean(results_sample$coef_estimate_s2 ), mean(results_sample$coef_estimate_s3))
  
  #bias
  bias_true[i,]  = rbind( abs( mean(results_sample$coef_true_s1)- beta[1]),
                          abs( mean(results_sample$coef_true_s2)-  beta[2]),
                          abs( mean(results_sample$coef_true_s3)-  beta[3]))
  bias_naive[i,]  = rbind( abs(mean(results_sample$coef_naive_s1)- beta[1]), 
                           abs(mean(results_sample$coef_naive_s2)-beta[2]),
                           abs(mean(results_sample$coef_naive_s3)-beta[3]))
  bias_w_sum1[i,]  = rbind(  abs(mean(results_sample$coef_w_sum1_s1)- beta[1]), 
                             abs( mean(results_sample$coef_w_sum1_s2)- beta[2]),
                             abs( mean(results_sample$coef_w_sum1_s3)- beta[3]))
  bias_w_sum2[i,]  = rbind(  abs(mean(results_sample$coef_w_sum2_s1)- beta[1]), 
                             abs( mean(results_sample$coef_w_sum2_s2)- beta[2]),
                             abs( mean(results_sample$coef_w_sum2_s3)- beta[3]))
  bias_estimate[i,]  = rbind( abs(mean(results_sample$coef_estimate_s1) -  beta[1]),
                              abs(mean(results_sample$coef_estimate_s2) -  beta[2] ),
                              abs(mean(results_sample$coef_estimate_s3) -  beta[3] ))
  
  # number of convergences
  L_naive[i] = length(which(results_sample$converge_naive==0))
  L_w_sum1[i] = length(which(results_sample$converge_w_sum1==0))
  L_w_sum2[i] = length( which(results_sample$converge_w_sum2==0))
  L_estimate[i] = length(which(results_sample$converge_estimate==0))
  
  # standard deviation
  Sd_true[i,] = rbind( sd(results_sample$coef_true_s1), sd(results_sample$coef_true_s2), sd(results_sample$coef_true_s3))
  Sd_naive[i,] = rbind( sd( results_sample$coef_naive_s1) , sd( results_sample$coef_naive_s2), sd( results_sample$coef_naive_s3) )
  Sd_sum1[i,] = rbind( sd(results_sample$coef_w_sum1_s1), sd(results_sample$coef_w_sum1_s2), sd(results_sample$coef_w_sum1_s3))
  Sd_sum2 [i,] = rbind( sd(results_sample$coef_w_sum2_s1),sd(results_sample$coef_w_sum2_s2),sd(results_sample$coef_w_sum2_s3))
  Sd_estimate[i,] = rbind( sd(results_sample$coef_estimate_s1),sd(results_sample$coef_estimate_s2),sd(results_sample$coef_estimate_s3))
  
  
  # mean square error
  rmse_true[i,] = rbind (sqrt( mean((results_sample$coef_true_s1- beta[1])^2)),
                         sqrt( mean((results_sample$coef_true_s2-  beta[2])^2)),
                         sqrt( mean((results_sample$coef_true_s3-  beta[3])^2)))
  
  rmse_naive[i,] = rbind( sqrt( mean((results_sample$coef_naive_s1- beta[1])^2)),
                          sqrt( mean((results_sample$coef_naive_s2-  beta[2])^2)),
                          sqrt( mean((results_sample$coef_naive_s3-  beta[3])^2)))
  
  rmse_sum1[i,] = rbind( sqrt( mean((results_sample$coef_w_sum1_s1- beta[1])^2)),
                         sqrt( mean((results_sample$coef_w_sum1_s2- beta[2])^2)),
                         sqrt( mean((results_sample$coef_w_sum1_s3- beta[3])^2)))
  
  rmse_sum2[i,] = rbind( sqrt( mean((results_sample$coef_w_sum2_s1- beta[1])^2)),
                         sqrt( mean((results_sample$coef_w_sum2_s2- beta[2])^2)),
                         sqrt( mean((results_sample$coef_w_sum2_s3- beta[3])^2)))
  
  rmse_estimate[i,] = rbind( sqrt( mean((results_sample$coef_estimate_s1 -  beta[1])^2)),
                             sqrt( mean((results_sample$coef_estimate_s2 -  beta[2])^2)),
                             sqrt( mean((results_sample$coef_estimate_s3 -  beta[3])^2)))
  
  
  Results_scenario = cbind( c(0,L_naive[i],L_w_sum1[i],L_w_sum2[i], L_estimate[i],L_estimate_Nd[i]),
                            c( beta_true[i,1],beta_naive[i,1],beta_w_sum1[i,1],beta_w_sum2[i,1],beta_estimate[i,1]),
                            c(bias_true[i,1],  bias_naive[i,1],bias_w_sum1[i,1],bias_w_sum2[i,1],bias_estimate[i,1]),
                            c(Sd_true[i,1], Sd_naive[i,1],Sd_sum1[i,1],Sd_sum2[i,1],Sd_estimate[i,1]),
                            c(rmse_true[i,1],rmse_naive[i,1], rmse_sum1[i,1], rmse_sum1[i,1],rmse_estimate[i,1]),
                            c( beta_true[i,2],beta_naive[i,2],beta_w_sum1[i,2],beta_w_sum2[i,2],beta_estimate[i,2]),
                            c(bias_true[i,2],  bias_naive[i,2],bias_w_sum1[i,2],bias_w_sum2[i,2],bias_estimate[i,2]),
                            c(Sd_true[i,2], Sd_naive[i,2],Sd_sum1[i,2],Sd_sum2[i,2],Sd_estimate[i,2]),
                            c(rmse_true[i,2],rmse_naive[i,2], rmse_sum1[i,2], rmse_sum1[i,2],rmse_estimate[i,2]),
                            c( beta_true[i,3],beta_naive[i,3],beta_w_sum1[i,3],beta_w_sum2[i,3],beta_estimate[i,3]),
                            c(bias_true[i,3],  bias_naive[i,3],bias_w_sum1[i,3],bias_w_sum2[i,3],bias_estimate[i,3]),
                            c(Sd_true[i,3], Sd_naive[i,3],Sd_sum1[i,3],Sd_sum2[i,3],Sd_estimate[i,3]),
                            c(rmse_true[i,3],rmse_naive[i,3], rmse_sum1[i,3], rmse_sum1[i,3],rmse_estimate[i,3]))
  
  colnames( Results_scenario) = c( "Fails","beta.1", "Bias.1", "Sd.1", "Rmse.1","beta.2", "Bias.2", "Sd.2", "Rmse.2",
                                   "beta.3", "Bias.3", "Sd.3", "Rmse.3")
  
    rownames( Results_scenario) = c("True"," Naive", "weighted average", "2-max weighted average", "EM")
    
    print(c("nA=",nA, "prob_1=", C[1] ))
    
  # print( Results_scenario)
    tb <- xtable(Results_scenario,digits=3,caption="Résultats sur les profils-lignes pour les 3 premiers axes factoriels")
   print(tb,include.rownames=TRUE,  sanitize.text.function = function(x){x},tabular.environment='longtable', 
         hline.after = c(-1,0),floating=FALSE,size="\\scriptsize")
   
  }
  
