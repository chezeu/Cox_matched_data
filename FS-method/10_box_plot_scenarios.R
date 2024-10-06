
######################## ####################
##############  
### boxplots

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/FS-method")

source("7_scenarios.R")

library(ggplot2)
library(gridExtra)
library(xtable)


  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/FS-method/Results/",
              "nsim=",nsim,"_nA=",nA,".Rdata"))
  results = results
function_plot <- function(nsim){
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(1, times=nsim),
                             results$coef_true_s1,
                             results$coef_true_s2,
                             results$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(2, times=nsim),
                              results$coef_naive_s1,
                              results$coef_naive_s2,
                              results$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(3,times=nsim),
                               results$coef_w_sum1_s1,
                               results$coef_w_sum1_s2,
                               results$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(4,times=nsim),
                               results$coef_w_sum2_s1,
                               results$coef_w_sum2_s2,
                               results$coef_w_sum2_s3)
  
  boxplot_beta0_estimate = cbind(rep(5,times=nsim),
                                 rep(5,times=nsim),
                                 results$coef_estimate_s1,
                                 results$coef_estimate_s2,
                                 results$coef_estimate_s3)
  
  boxplot_beta0 = rbind(boxplot_beta0_true,boxplot_beta0_naive,boxplot_beta0_w_sum1,
                        boxplot_beta0_w_sum2,boxplot_beta0_estimate )
  
  
return(boxplot_beta0)
}
  
boxplot_beta0 = function_plot (nsim)  


  titlename = paste0("nB=",nB, ", ","nA=",nA)
colnames(boxplot_beta0) = c("Estimator","position","beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
             levels = 1:5, labels = c("True"," Naive", "weighted average", "2-max weighted average", "EM"))


boxplot_beta0[,3] = as.numeric(boxplot_beta0[,3] )
boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
attach(boxplot_beta0)
# Boxplots for beta0

beta0_1 = ggplot(boxplot_beta0,aes(x=position,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(position)+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")+
  ggtitle(titlename)
beta0_1

beta0_2 = ggplot(boxplot_beta0,aes(x=position,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  xlab(position)+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")+
  ggtitle(titlename)
beta0_2

beta0_3 = ggplot(boxplot_beta0,aes(x=position,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  xlab(position)+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")+
  ggtitle(titlename)
beta0_3
