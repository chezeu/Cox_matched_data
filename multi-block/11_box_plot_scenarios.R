
######################## ####################
##############     m fixed and Q varies            ######
### boxplots

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/multi-block")

source("09_scenarios.R")

#scenarios_boxplot=scenarios[scenarios$nA==100,] 
#scenarios_boxplot=scenarios[scenarios$nA==500,] 
scenarios_boxplot=scenarios[scenarios$nA==100,] 

varies = "K"
fixed="nA"

library(ggplot2)
library(gridExtra)
library(xtable)

boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  K = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  nB= scenarios_boxplot[i,4]
  nA_block= c(scenarios[i,5],scenarios[i,6],scenarios[i,7] )
  nB_block= c(scenarios[i,8],scenarios[i,9],scenarios[i,10] )
  
  filename = paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/multi-block/Results/","nsim=",nsim,
                    "_nA=",nA, "_K=",K,".Rdata")
  results_sample = results_sample
  
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(K,times=nsim),
                             rep(nA, times=nsim),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(K,times=nsim),
                              rep(nA, times=nsim),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(K,times=nsim),
                               rep(nA,times=nsim),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(K,times=nsim),
                               rep(nA,times=nsim),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
  
  boxplot_beta0_estimate = cbind(rep(5,times=nsim),
                                 rep(K,times=nsim),
                                 rep(nA, times=nsim),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind( boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,
                        boxplot_beta0_w_sum1,boxplot_beta0_w_sum2, 
                        boxplot_beta0_estimate  )
  
}


titlename =paste ("nA=", nA, "nB=", nB)

colnames(boxplot_beta0) = c("Estimator", varies, fixed,"beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:5, labels = c("True"," Naive", "weighted average", "2-max weighted average", "EM"))

boxplot_beta0[,2] = as.factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = as.factor(boxplot_beta0[,3] )

boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )
# Boxplots for beta0

beta0_1 = ggplot(boxplot_beta0,aes(x=K,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")+
  ggtitle(titlename)
beta0_1

beta0_2 = ggplot(boxplot_beta0,aes(x=K,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")+
 ggtitle(titlename)
beta0_2

beta0_3 = ggplot(boxplot_beta0,aes(x=K,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")+
  ggtitle(titlename)
beta0_3

# histogram


#data_histogram = boxplot_beta0[which(boxplot_beta0$Estimator==unique(boxplot_beta0$Estimator)[3]),]
#K_hist = unique(data_histogram$K)

#for (i in 1:length(K_hist)) {
 # k= K_hist[i]
#  beta_hist =data_histogram[which(data_histogram$K==k),]
#  hist(beta_hist$beta0.1, xlab= quote(hat(beta)[1]),main = paste("nA =",unique( beta_hist$nA),",","K=",unique( beta_hist$K)))
 # hist(beta_hist$beta0.2, xlab= quote(hat(beta)[2]),main = paste("nA =",unique( beta_hist$nA),",","K=",unique( beta_hist$K)))
  #hist(beta_hist$beta0.3, xlab= quote(hat(beta)[3]),main = paste("nA =",unique( beta_hist$nA),",","K=",unique( beta_hist$K)))
#}
#########################################################
########################################

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/multi-block")

source("09_scenarios.R")

scenarios_boxplot=scenarios[scenarios$K==7,] 
#scenarios_boxplot=scenarios[scenarios$K==14,] 
#scenarios_boxplot=scenarios[scenarios$K==21,] 

varies = "nA"
fixed="K"

library(ggplot2)
library(gridExtra)
library(xtable)


boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  K = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  nB= scenarios_boxplot[i,4]
  nA_block= c(scenarios[i,5],scenarios[i,6],scenarios[i,7] )
  nB_block= c(scenarios[i,8],scenarios[i,9],scenarios[i,10] )
  
  filename = paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/multi-block/Results/","nsim=",nsim,
                    "_nA=",nA, "_K=",K,".Rdata")
  results_sample = results_sample
  
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(K,times=nsim),
                             rep(nA, times=nsim),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(K,times=nsim),
                              rep(nA, times=nsim),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(K,times=nsim),
                               rep(nA,times=nsim),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(K,times=nsim),
                               rep(nA,times=nsim),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
  
  boxplot_beta0_estimate = cbind(rep(5,times=nsim),
                                 rep(K,times=nsim),
                                 rep(nA, times=nsim),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind( boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,
                         boxplot_beta0_w_sum1,boxplot_beta0_w_sum2, 
                         boxplot_beta0_estimate  )
  
}


titlename =paste ("K=", K)

colnames(boxplot_beta0) = c("Estimator",fixed, varies, "beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:5, labels = c("True"," Naive", "weighted average", "2-max weighted average", "EM"))
boxplot_beta0[,2] = as.factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = as.factor(boxplot_beta0[,3] )

boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )
# Boxplots for beta0

beta0_1 = ggplot(boxplot_beta0,aes(x=nA,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")+
  ggtitle(titlename)
beta0_1

beta0_2 = ggplot(boxplot_beta0,aes(x=nA,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")+
  ggtitle(titlename)
beta0_2

beta0_3 = ggplot(boxplot_beta0,aes(x=nA,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")+
  ggtitle(titlename)
beta0_3
#########################################################
############# Histogram


data_histogram = boxplot_beta0[which(boxplot_beta0$Estimator==unique(boxplot_beta0$Estimator)[3]),]
nA_hist = unique(data_histogram$nA)

for (i in 1:length(nA_hist)) {
  n = nA_hist[i]
  beta_hist =data_histogram[which(data_histogram$nA==n),]
  hist(beta_hist$beta0.1, xlab= quote(hat(beta)[1]),main = paste("nA =", n,",","K=",unique(data_histogram$K)))
  hist(beta_hist$beta0.2, xlab= quote(hat(beta)[2]),main = paste("nA =", n,",","K=",unique(data_histogram$K)))
  hist(beta_hist$beta0.2, xlab= quote(hat(beta)[3]),main = paste("nA =", n,",","K=",unique(data_histogram$K)))
}
