

### boxplots

setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")

source("8_scenarios.R")

scenarios_boxplot=scenarios # the last scenario 6

library(ggplot2)
library(gridExtra)

levelsvec = numeric(0)


file_list <- c("scenario=1_m=1000_prob1=0.6.Rdata", 
               "scenario=2_m=1000_prob1=0.7.Rdata",  
               "scenario=3_m=1000_prob1=0.8.Rdata",
               "scenario=4_m=500_prob1=0.8.Rdata",
               "scenario=5_m=1000_prob1=0.8.Rdata",
               "scenario=6_m=2000_prob1=0.8.Rdata") 

for (i in (1:nrow(scenarios_boxplot))){
  nsim = scenarios_boxplot[i,1]
  n = scenarios_boxplot[i,2]
  m = scenarios_boxplot[i,3]
  Prob1 = scenarios_boxplot[i,4]
  Prob2 = scenarios_boxplot[i,5]
  Prob3 = scenarios_boxplot[i,6]

  Prob_sample = paste0("(",Prob1,",",Prob2,",",Prob3,")") 
  
  setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu/Results")
  load(file_list[i])
  results_sample = results_sample
  
  levelsvec = c(levelsvec,Prob_sample)

  boxplot_beta0_true = cbind(rep(1,times=nsim),
                               rep(Prob_sample,times=nsim),
                               results_sample$coef_true_s1,
                               results_sample$coef_true_s2)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                      rep(Prob_sample,times=nsim),
                      results_sample$coef_naive_s1,
                      results_sample$coef_naive_s2)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                      rep(Prob_sample,times=nsim),
                      results_sample$coef_w_sum1_s1,
                      results_sample$coef_w_sum1_s2)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                                rep(Prob_sample,times=nsim),
                                results_sample$coef_w_sum2_s1,
                                results_sample$coef_w_sum2_s2)
  
  boxplot_beta0_estimate = cbind(rep(5,times=nsim),
                                rep(Prob_sample,times=nsim),
                                results_sample$coef_estimate_s1,
                                results_sample$coef_estimate_s2)
  
  boxplot_beta0 = rbind(boxplot_beta0_true,boxplot_beta0_naive,boxplot_beta0_w_sum1,
                        boxplot_beta0_w_sum2,boxplot_beta0_estimate )
}

titlename = paste0("n=",unique(scenarios_boxplot[,2]),
                   ", ","m=",unique(scenarios_boxplot[,3]))

colnames(boxplot_beta0) = c("Estimator", "Prob","beta0[1]","beta0[2]")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
      levels = 1:5, labels = c("True"," Naive", "w_sum1", "w_sum2", "EM"))

boxplot_beta0[,2] = factor(boxplot_beta0[,2],
                         levels = levelsvec,labels = levelsvec)

boxplot_beta0[,3] = as.numeric(boxplot_beta0[,3] )
boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
# Boxplots for beta0

beta0_1 = ggplot(boxplot_beta0,aes(x=Prob,y= beta0.1.,fill=Estimator))+
  geom_boxplot()+
  xlab("Q_matrix")+
  ylab(quote(hat(beta)[0]))+
  geom_hline(yintercept=0.5,lty=1,col="orange")+
  ggtitle(titlename)
beta0_1

beta0_2 = ggplot(boxplot_beta0,aes(x=Prob,y= beta0.2.,fill=Estimator))+
  geom_boxplot()+
  xlab("Q_matrix")+
  ylab(quote(hat(beta)[0]))+
  geom_hline(yintercept=-0.5,lty=1,col="orange")+
  ggtitle(titlename)
beta0_2
