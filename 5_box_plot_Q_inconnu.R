
######################## ####################
##############     nA fixed and K varies            ######
### boxplots
setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/data_generation")

source("3_scenarios_Q_inconnu.R")

#scenarios_boxplot=scenarios[scenarios$nB==200,] 
#scenarios_boxplot=scenarios[scenarios$nA==1000,] 
scenarios_boxplot=scenarios[(scenarios$nB==2000)&
                              (scenarios$censor==2.948469),] 

varies = "K"
fixed="nA"

library(ggplot2)
library(gridExtra)
library(xtable)
library(ggpubr)
#theme_set(theme_bw() + theme(legend.position = "top"))
theme_set(theme(legend.position = "top"))

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results")

boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  K = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  nB= scenarios_boxplot[i,4]
  censor=scenarios_boxplot[i,5]
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results/Results_Q_inconnu/",
              "nsim=",nsim,"_nA=",nA,"_K=",K,"censor=",censor, ".Rdata"))
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
  
  boxplot_beta0_Q_sum = cbind(rep(5,times=nsim),
                              rep(K,times=nsim),
                              rep(nA,times=nsim),
                              results_sample$coef_Q_sum_s1,
                              results_sample$coef_Q_sum_s2,
                              results_sample$coef_Q_sum_s3)
  
  boxplot_beta0_estimate = cbind(rep(6,times=nsim),
                                 rep(K,times=nsim),
                                 rep(nA, times=nsim),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,
                        boxplot_beta0_w_sum1,  boxplot_beta0_w_sum2,
                        boxplot_beta0_Q_sum,boxplot_beta0_estimate  )
  
}


#titlename =paste ("nA=", nA, "nB=", nB)

colnames(boxplot_beta0) = c("Estimator", varies, fixed,"beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:6, labels = c("Theorical"," Naive", "weighted average", "2-max weighted average", "Q_weighted average","EM"))

boxplot_beta0[,2] = as.factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = as.factor(boxplot_beta0[,3] )

boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )

# Boxplots for beta0

beta0_1_K = ggplot(boxplot_beta0,aes(x=K,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")
#beta0_1

beta0_2_K = ggplot(boxplot_beta0,aes(x=K,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(varies)+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
#beta0_2

beta0_3_K = ggplot(boxplot_beta0,aes(x=K,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(varies)+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
#beta0_3

#figure <- ggarrange(beta0_1_K, beta0_2_K, beta0_3_K,
#           common.legend = TRUE,
#                    ncol = 2, nrow = 2)
#figure


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
################### censor varied 

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/data_generation")

source("3_scenarios_Q_inconnu.R")
#scenarios_boxplot=scenarios[scenarios$K==9,] 
scenarios_boxplot=scenarios[(scenarios$K==16)&
                              ( scenarios$nB==2000),] 
#scenarios_boxplot=scenarios[scenarios$K==27,] 

varies = "censor"
fixed="nA"

library(ggplot2)
library(gridExtra)
library(xtable)
setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results")

Censor_sample=c(30,40,20)# for 30%,40%,20%

boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  K = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  nB= scenarios_boxplot[i,4]
  censor=scenarios_boxplot[i,5]
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results/Results_Q_inconnu/",
              "nsim=",nsim,"_nA=",nA,"_K=",K,"censor=",censor,".Rdata"))
  results_sample = results_sample
  
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(Censor_sample[i],times=nsim),
                             rep(nA, times=nsim),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(Censor_sample[i],times=nsim),
                              rep(nA, times=nsim),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(Censor_sample[i],times=nsim),
                               rep(nA,times=nsim),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(Censor_sample[i],times=nsim),
                               rep(nA,times=nsim),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
  
  boxplot_beta0_Q_sum = cbind(rep(5,times=nsim),
                              rep(Censor_sample[i],times=nsim),
                              rep(nA,times=nsim),
                              results_sample$coef_Q_sum_s1,
                              results_sample$coef_Q_sum_s2,
                              results_sample$coef_Q_sum_s3)
  
  boxplot_beta0_estimate = cbind(rep(6,times=nsim),
                                 rep(Censor_sample[i],times=nsim),
                                 rep(nA,times=nsim),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,
                        boxplot_beta0_w_sum1,  boxplot_beta0_w_sum2,
                        boxplot_beta0_Q_sum,boxplot_beta0_estimate  )
  
}



#titlename =paste ("K=", K)

colnames(boxplot_beta0) = c("Estimator", varies,fixed, "beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:6, labels = c("Theorical"," Naive", "weighted average", "2-max weighted average", "Q_weighted average","EM"))

boxplot_beta0[,2] = as.factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = as.factor(boxplot_beta0[,3] )

boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )
# Boxplots for beta0

beta0_1_censor = ggplot(boxplot_beta0,aes(x=censor,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(nu))+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")
#+ggtitle(titlename)
#beta0_1

beta0_2_censor = ggplot(boxplot_beta0,aes(x=censor,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(nu))+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
#+  ggtitle(titlename)
#beta0_2

beta0_3_censor = ggplot(boxplot_beta0,aes(x=censor,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(nu))+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
#+  ggtitle(titlename)
#beta0_3
#figure <- ggarrange(beta0_1_censor, beta0_2_censor, beta0_3_censor,
 #                   common.legend = TRUE,
  #                  ncol = 2, nrow = 2)
#figure

#########################################################
################### nA varied 

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/data_generation")

source("3_scenarios_Q_inconnu.R")
#scenarios_boxplot=scenarios[scenarios$K==9,] 
scenarios_boxplot=scenarios[(scenarios$censor==2.948469 & scenarios$K==16)&
                           (scenarios$nB==2*scenarios$nA),] 
#scenarios_boxplot=scenarios[scenarios$K==27,] 

varies = "nA"
fixed="K"

library(ggplot2)
library(gridExtra)
library(xtable)

#Censor_sample=c(30,40,20)# for 30%,40%,20%
setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results")

boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  K = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  nB= scenarios_boxplot[i,4]
  censor=scenarios_boxplot[i,5]
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results/Results_Q_inconnu/",
              "nsim=",nsim,"_nA=",nA,"_K=",K,"censor=",censor,".Rdata"))
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
  
  boxplot_beta0_Q_sum = cbind(rep(5,times=nsim),
                              rep(K,times=nsim),
                              rep(nA,times=nsim),
                              results_sample$coef_Q_sum_s1,
                              results_sample$coef_Q_sum_s2,
                              results_sample$coef_Q_sum_s3)
  
  boxplot_beta0_estimate = cbind(rep(6,times=nsim),
                                 rep(K,times=nsim),
                                 rep(nA,times=nsim),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,
                        boxplot_beta0_w_sum1,  boxplot_beta0_w_sum2,
                        boxplot_beta0_Q_sum,boxplot_beta0_estimate  )
  
}



#titlename =paste ("K=", K)

colnames(boxplot_beta0) = c("Estimator",fixed, varies, "beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:6, labels = c("Theorical"," Naive", "weighted average", "2-max weighted average","Q_weighted average", "EM"))

boxplot_beta0[,2] = as.factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = as.factor(boxplot_beta0[,3] )

boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )
# Boxplots for beta0

beta0_1_nA = ggplot(boxplot_beta0,aes(x=nA,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")
#+  ggtitle(titlename)
#beta0_1

beta0_2_nA = ggplot(boxplot_beta0,aes(x=nA,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
#+  ggtitle(titlename)
#beta0_2

beta0_3_nA = ggplot(boxplot_beta0,aes(x=nA,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
#+  ggtitle(titlename)
#beta0_3
#beta0_3
#figure <- ggarrange(beta0_1_nA, beta0_2_nA, beta0_3_nA,
 #                   common.legend = TRUE,
  #                  ncol = 2, nrow = 2)
#figure

#########

#########################################################
################### R varied 

#scenarios_boxplot=scenarios[scenarios$K==9,] 
scenarios_boxplot=scenarios[(scenarios$K==16 & scenarios$censor!=1.553309)&
                            (scenarios$nA==1000& scenarios$censor!=6.392605), ] 
#scenarios_boxplot=scenarios[scenarios$K==27,] 

varies = "R"
fixed="nA"

library(ggplot2)
library(gridExtra)
library(xtable)

#R=c(2,4,6) rapport de taille

boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  K = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  nB= scenarios_boxplot[i,4]
  censor=scenarios_boxplot[i,5]
  R=nB/nA
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results/Results_Q_inconnu/",
              "nsim=",nsim,"_nA=",nA,"_nB=",nB,"_K=",K,"censor=",censor,".Rdata"))
  results_sample = results_sample
  
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(R,times=nsim),
                             rep(nB, times=nsim),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(R,times=nsim),
                              rep(nB, times=nsim),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(R,times=nsim),
                               rep(nB,times=nsim),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(R,times=nsim),
                               rep(nB,times=nsim),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
  
  boxplot_beta0_Q_sum = cbind(rep(5,times=nsim),
                              rep(R,times=nsim),
                              rep(nB,times=nsim),
                              results_sample$coef_Q_sum_s1,
                              results_sample$coef_Q_sum_s2,
                              results_sample$coef_Q_sum_s3)
  
  boxplot_beta0_estimate = cbind(rep(6,times=nsim),
                                 rep(R,times=nsim),
                                 rep(nB,times=nsim),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,
                        boxplot_beta0_w_sum1,  boxplot_beta0_w_sum2,
                        boxplot_beta0_Q_sum,boxplot_beta0_estimate  )
  
}



#titlename =paste ("nA=", nA)

colnames(boxplot_beta0) = c("Estimator", varies,"nB", "beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:6, labels = c("Theorical"," Naive", "weighted average", "2-max weighted average","Q_weighted average", "EM"))

boxplot_beta0[,2] = as.factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = as.factor(boxplot_beta0[,3] )

boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )
# Boxplots for beta0

beta0_1_R = ggplot(boxplot_beta0,aes(x=R,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(R))+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")
#+  ggtitle(titlename)
#beta0_1_R

beta0_2_R = ggplot(boxplot_beta0,aes(x=R,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(R))+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
#+  ggtitle(titlename)
#beta0_2_R

beta0_3_R = ggplot(boxplot_beta0,aes(x=R,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(R))+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
#+  ggtitle(titlename)
#beta0_3_R
#figure <- ggarrange(beta0_1_R, beta0_2_R, beta0_3_R,
#                    common.legend = TRUE,
#                    ncol = 2, nrow = 2)
#figure
figure <- ggarrange(beta0_1_K, beta0_2_K, beta0_3_K,
                    beta0_1_nA, beta0_2_nA, beta0_3_nA,
                    beta0_1_censor, beta0_2_censor, beta0_3_censor,
                    beta0_1_R, beta0_2_R, beta0_3_R,
                    common.legend = TRUE,
                    ncol = 3, nrow = 4)

figure
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
