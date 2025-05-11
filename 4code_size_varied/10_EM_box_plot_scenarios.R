
######################## ####################
##############  K varies            ######
### boxplots

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/4code_size_varied")

source("7_EM_scenarios.R")

#scenarios_boxplot=scenarios[scenarios$nB==2000,] 
#scenarios_boxplot=scenarios[scenarios$nB==4000,] 
scenarios_boxplot=scenarios[(scenarios$nB==4000)&
                              (scenarios$K!=18),] 

varies = "K"
fixed="nB"

library(ggplot2)
library(gridExtra)
library(xtable)


boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  K= scenarios_boxplot[i,4]
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/code_size_varied/Results/",
              "nsim=",nsim,"_nB=",nB,"_nA=",nA,"_K=",K,".Rdata"))
  results_sample = results_sample
  
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(K,times=nsim),
                             rep(nB, times=nsim),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(K,times=nsim),
                              rep(nB, times=nsim),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(K,times=nsim),
                               rep(nB,times=nsim),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(K,times=nsim),
                               rep(nB,times=nsim),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
  
  boxplot_beta0_estimate = cbind(rep(5,times=nsim),
                                 rep(K,times=nsim),
                                 rep(nB, times=nsim),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,
                        boxplot_beta0_w_sum1,  boxplot_beta0_w_sum2,boxplot_beta0_estimate  )
  
}


#titlename =paste ("nA=", nA, "nB=", nB)

colnames(boxplot_beta0) = c("Estimator", varies, fixed,"beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:5, labels = c("Theorical"," Naive", "weighted average", "2-max weighted average", "EM"))

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
  geom_hline(yintercept=beta[1],lty=1,col="orange")
beta0_1

beta0_2 = ggplot(boxplot_beta0,aes(x=K,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(varies)+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
beta0_2

beta0_3 = ggplot(boxplot_beta0,aes(x=K,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(varies)+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
beta0_3

#########################################################
########### nB varied ###

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/4code_size_varied")

#source("5_EM_scenarios.R")
source("7_EM_scenarios.R")
#scenarios_boxplot=scenarios[scenarios$K==9,] 
scenarios_boxplot=scenarios[scenarios$K==18,] 
#scenarios_boxplot=scenarios[scenarios$K==27,] 

varies = "nB"
fixed="K"

library(ggplot2)
library(gridExtra)
library(xtable)


boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  K= scenarios_boxplot[i,4]
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/code_size_varied/Results/",
              "nsim=",nsim,"_nB=",nB,"_nA=",nA,"_K=",K,".Rdata"))
  results_sample = results_sample
  
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(K,times=nsim),
                             rep(nB, times=nsim),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(K,times=nsim),
                              rep(nB, times=nsim),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(K,times=nsim),
                               rep(nB,times=nsim),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(K,times=nsim),
                               rep(nB,times=nsim),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
  
  boxplot_beta0_estimate = cbind(rep(5,times=nsim),
                                 rep(K,times=nsim),
                                 rep(nB, times=nsim),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,
                        boxplot_beta0_w_sum1,  boxplot_beta0_w_sum2,boxplot_beta0_estimate  )
  
}



titlename =paste ("K=", K)

colnames(boxplot_beta0) = c("Estimator",fixed, varies, "beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:5, labels = c("Theorical"," Naive", "weighted average", "2-max weighted average", "EM"))

boxplot_beta0[,2] = as.factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = as.factor(boxplot_beta0[,3] )

boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )
# Boxplots for beta0

beta0_1 = ggplot(boxplot_beta0,aes(x=nB,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")+  ggtitle(titlename)
beta0_1

beta0_2 = ggplot(boxplot_beta0,aes(x=nB,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(varies)+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")+
  ggtitle(titlename)
beta0_2

beta0_3 = ggplot(boxplot_beta0,aes(x=nB,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
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
