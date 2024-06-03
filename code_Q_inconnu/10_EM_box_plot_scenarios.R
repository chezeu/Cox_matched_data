
######################## ####################
##############     m fixed and Q varies            ######
### boxplots

setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_inconnu")

source("7_EM_scenarios.R")

varies = "nA"
fixed="beta"

library(ggplot2)
library(gridExtra)
library(xtable)


boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios))){
  
  nsim = scenarios[i,1]
  nB = scenarios[i,2]
  nA = scenarios[i,3]
 
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_inconnu/Results/",
              "nsim=",nsim,"_nA=",nA,"_nB=",nB, ".Rdata"))
  results_sample = results_sample
  
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(nA,times=nsim),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(nA,times=nsim),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2)
  
  
  boxplot_beta0_estimate = cbind(rep(3,times=nsim),
                                 rep(nA,times=nsim),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,
         boxplot_beta0_estimate  )
  
}


titlename =paste ("error=", error)

colnames(boxplot_beta0) = c("Estimator",varies, "beta0.1","beta0.2")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:3, labels = c("True"," Naive","EM"))
boxplot_beta0[,2] = as.factor(boxplot_beta0[,2] )

boxplot_beta0[,3] = as.numeric(boxplot_beta0[,3] )
boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
# Boxplots for beta0

beta0_1 = ggplot(boxplot_beta0,aes(x=nA,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")+
  ggtitle(titlename)
beta0_1

beta0_2 = ggplot(boxplot_beta0,aes(x=nA,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")+
 ggtitle(titlename)
beta0_2

