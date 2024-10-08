
######################## ####################
##############     m fixed and Q varies            ######
### boxplots

setwd("H:/code_record")

source("7_scenarios.R")

scenarios_boxplot=scenarios[scenarios$nA==100,] 
#scenarios_boxplot=scenarios[scenarios$nA==200,] 
#scenarios_boxplot=scenarios[scenarios$nA==300,] 

varies = "Q"
fixed = "nA"

library(ggplot2)
library(gridExtra)
library(xtable)


boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios_boxplot ))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  Prob1 = scenarios_boxplot[i,4]
  Prob2 = scenarios_boxplot[i,5]
  Prob3 = scenarios_boxplot[i,6]
  
  Prob_sample = paste0("(",Prob1,",",Prob2,",",Prob3,")") 
  
  load(paste0("H:/code_record/Results/",
              "nsim=",nsim,"_nA=",nA, "_prob1=", Prob1, ".Rdata"))
  results_sample = results_sample
  
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(nA,times=nsim),
                             rep( Prob_sample , times=nsim ),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3,
                             results_sample$coef_true_s4)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(nA,times=nsim),
                              rep(  Prob_sample , times=nsim ),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3,
                              results_sample$coef_naive_s4)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(nA,times=nsim),
                               rep(  Prob_sample , times=nsim ),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3,
                               results_sample$coef_w_sum1_s4)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(nA,times=nsim),
                               rep(  Prob_sample , times=nsim ),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3,
                               results_sample$coef_w_sum2_s4)
  
  boxplot_beta0_estimate = cbind(rep(5,times=nsim),
                                 rep(nA,times=nsim),
                                 rep(  Prob_sample , times=nsim ),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3,
                                 results_sample$coef_estimate_s4)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,boxplot_beta0_w_sum1,
                        boxplot_beta0_w_sum2,boxplot_beta0_estimate  )
  
}

titlename = paste0("nB=",unique(scenarios_boxplot[,2]),
                   ", ","nA=",unique(scenarios_boxplot[,3]))

colnames(boxplot_beta0) = c("Estimator",fixed, varies, "beta0.1","beta0.2","beta0.3","beta0.4")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
           levels = 1:5, labels = c("True"," Naive", "w_sum1", "w_sum2", "EM"))


boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4])
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5])
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6])
boxplot_beta0[,7] = as.numeric(boxplot_beta0[,7])
# Boxplots for beta0

beta0_1 = ggplot(boxplot_beta0,aes(x=Q,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=0.5,lty=1,col="orange")+
  ggtitle(titlename)
beta0_1
beta0_2 = ggplot(boxplot_beta0,aes(x=Q,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=-0.5,lty=1,col="orange")+
  ggtitle(titlename)
beta0_2
beta0_3 = ggplot(boxplot_beta0,aes(x=Q,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=-0.5,lty=1,col="orange")+
  ggtitle(titlename)
beta0_3
beta0_4 = ggplot(boxplot_beta0,aes(x=Q,y= beta0.4,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[4]))+
  geom_hline(yintercept=-0.5,lty=1,col="orange")+
  ggtitle(titlename)
beta0_4
######################## ####################

##############     Q fixed and m varies            ######
### boxplots


setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")

source("7_scenarios.R")

 scenarios_boxplot=scenarios[scenarios$Prob_1==0.4,] 
#scenarios_boxplot=scenarios[scenarios$Prob_1==0.7,] 
#scenarios_boxplot=scenarios[scenarios$Prob_1==1,] 

varies="nA"
fixed = "Prob_1"

library(ggplot2)
library(gridExtra)

boxplot_beta0 = NULL
# size fix
for (i in (1:nrow(scenarios_boxplot ))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  Prob1 = scenarios_boxplot[i,4]
  Prob2 = scenarios_boxplot[i,5]
  Prob3 = scenarios_boxplot[i,6]
  
  Prob_sample = paste0("(",Prob1,",",Prob2,",",Prob3,")") 
  
  load(paste0("H:/code_record/Results/",
              "nsim=",nsim,"_nA=",nA, "_prob1=", Prob1, ".Rdata"))
  results_sample = results_sample
  
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(nA,times=nsim),
                             rep( Prob_sample , times=nsim ),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3,
                             results_sample$coef_true_s4)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(nA,times=nsim),
                              rep(  Prob_sample , times=nsim ),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3,
                              results_sample$coef_naive_s4)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(nA,times=nsim),
                               rep(  Prob_sample , times=nsim ),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3,
                               results_sample$coef_w_sum1_s4)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(nA,times=nsim),
                               rep(  Prob_sample , times=nsim ),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3,
                               results_sample$coef_w_sum2_s4)
  
  boxplot_beta0_estimate = cbind(rep(5,times=nsim),
                                 rep(nA,times=nsim),
                                 rep(  Prob_sample , times=nsim ),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3,
                                 results_sample$coef_estimate_s4)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,boxplot_beta0_w_sum1,
                        boxplot_beta0_w_sum2,boxplot_beta0_estimate )
}

titlename = paste0( "Q=", "(",unique(scenarios_boxplot[,4]), ",",
                    unique(scenarios_boxplot[,5]), ",",unique(scenarios_boxplot[,6]), ")" )

colnames(boxplot_beta0) = c("Estimator",varies,fixed,"beta0.1","beta0.2","beta0.3","beta0.4")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
      levels = 1:5, labels = c("True"," Naive", "w_sum1", "w_sum2", "EM"))

boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )
boxplot_beta0[,7] = as.numeric(boxplot_beta0[,7] )

# Boxplots for beta0

beta0_1 = ggplot(boxplot_beta0,aes(x=nA,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=0.5,lty=1,col="orange")+
  ggtitle(titlename)
beta0_1
beta0_2 = ggplot(boxplot_beta0,aes(x=nA,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=-0.5,lty=1,col="orange")+
  ggtitle(titlename)
beta0_2
beta0_3 = ggplot(boxplot_beta0,aes(x=nA,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=-0.5,lty=1,col="orange")+
  ggtitle(titlename)
beta0_3
beta0_4 = ggplot(boxplot_beta0,aes(x=nA,y= beta0.4,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[4]))+
  geom_hline(yintercept=-0.5,lty=1,col="orange")+
  ggtitle(titlename)
beta0_4
#####################################

