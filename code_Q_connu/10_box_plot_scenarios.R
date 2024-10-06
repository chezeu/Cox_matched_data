
######################## ####################
##############     m fixed and Q varies            ######
### boxplots

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/code_Q_connu")

source("7_scenarios.R")

#scenarios_boxplot=scenarios[scenarios$nA==100,] 
#scenarios_boxplot=scenarios[scenarios$nA==500,] 
scenarios_boxplot=scenarios[scenarios$nA==1000,] 

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
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/code_Q_connu/Results/",
              "nsim=",nsim,"_nA=",nA, "_prob1=", Prob1, ".Rdata"))
  results_sample = results_sample
  
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep( Prob_sample ,times=nsim),
                             rep(nA, times=nsim),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep( Prob_sample ,times=nsim),
                              rep(nA, times=nsim),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(  Prob_sample , times=nsim ),
                               rep(nA,times=nsim),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(  Prob_sample , times=nsim ),
                               rep(nA,times=nsim),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
  
  boxplot_beta0_estimate = cbind(rep(5,times=nsim),
                                 rep(  Prob_sample , times=nsim ),
                                 rep(nA,times=nsim),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,boxplot_beta0_w_sum1,
                        boxplot_beta0_w_sum2,boxplot_beta0_estimate )
}

titlename = paste0("nB=",unique(scenarios_boxplot[,2]),
                ", ","nA=",unique(scenarios_boxplot[,3]))

colnames(boxplot_beta0) = c("Estimator",varies,fixed,  "beta0.1","beta0.2")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
             levels = 1:5, labels = c("True"," Naive", "weighted average", "2-max weighted average", "EM"))


boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
#boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )

# Boxplots for beta0

beta0_1 = ggplot(boxplot_beta0,aes(x=Q,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")+
  ggtitle(titlename)
beta0_1
beta0_2 = ggplot(boxplot_beta0,aes(x=Q,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")+
  ggtitle(titlename)
beta0_2
beta0_3 = ggplot(boxplot_beta0,aes(x=Q,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  xlab(varies)+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")+
  ggtitle(titlename)
beta0_3
# histogram


data_histogram = boxplot_beta0[which(boxplot_beta0$Estimator==unique(boxplot_beta0$Estimator)[5]),]
Q_hist = unique(data_histogram$Q)

for (i in 1:length(Q_hist)) {
  q= Q_hist[i]
  beta_hist =data_histogram[which(data_histogram$Q==q),]
  hist(beta_hist$beta0.1, xlab= quote(hat(beta)[1]),main = paste("nA =",unique( beta_hist$nA),",","Q=",unique( beta_hist$Q)))
  hist(beta_hist$beta0.2, xlab= quote(hat(beta)[2]),main = paste("nA =",unique( beta_hist$nA),",","Q=",unique( beta_hist$Q)))
}

################## ####################

##############     Q fixed and m varies            ######
### boxplots,main = paste("nA =", n,",","prob=",unique(data_histogram$Prob_1)))


setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")

source("7_scenarios.R")

scenarios_boxplot=scenarios[scenarios$Prob_1==0.4,] 
#scenarios_boxplot=scenarios[scenarios$Prob_1==0.8,] 
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
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/code_Q_connu/Results/",
              "nsim=",nsim,"_nA=",nA, "_prob1=", Prob1, ".Rdata"))
  results_sample = results_sample
  
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(nA,times=nsim),
                             rep(  Prob_sample , times=nsim ),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(nA,times=nsim),
                              rep(  Prob_sample , times=nsim ),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(nA,times=nsim),
                               rep(  Prob_sample , times=nsim ),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(nA,times=nsim),
                               rep(  Prob_sample , times=nsim ),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2)
  
  boxplot_beta0_estimate = cbind(rep(5,times=nsim),
                                 rep(nA,times=nsim),
                                 rep(  Prob_sample , times=nsim ),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2)

  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,boxplot_beta0_w_sum1,
                        boxplot_beta0_w_sum2,boxplot_beta0_estimate )
  
  
}

titlename = paste0( "Q=", "(",unique(scenarios_boxplot[,4]), ",",
                    unique(scenarios_boxplot[,5]), ",",unique(scenarios_boxplot[,6]), ")" )

colnames(boxplot_beta0) = c("Estimator",varies,fixed,  "beta0.1","beta0.2")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
       levels = 1:5,labels = c("True"," Naive", "weighted average", "2-max weighted average", "EM"))

boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
order=unique(boxplot_beta0[,2])
# Boxplots for beta0

beta0_1 = ggplot(boxplot_beta0,aes(x=nA,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")+
  ggtitle(titlename)
beta0_1 + scale_x_discrete(limits=c(order[1],order[2],order[3]))

beta0_2 = ggplot(boxplot_beta0,aes(x=nA,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")+
  ggtitle(titlename)
beta0_2+ scale_x_discrete(limits=c(order[1],order[2],order[3]))
#####################################
### Histogram


data_histogram = boxplot_beta0[which(boxplot_beta0$Estimator==unique(boxplot_beta0$Estimator)[5]),]
nA_hist = unique(data_histogram$nA)

for (i in 1:length(nA_hist)) {
  n = nA_hist[i]
  beta_hist =data_histogram[which(data_histogram$nA==n),]
  hist(beta_hist$beta0.1, xlab= quote(hat(beta)[1]),main = paste("nA =", n,",","prob=",unique(data_histogram$Prob_1)))
  hist(beta_hist$beta0.2, xlab= quote(hat(beta)[2]),main = paste("nA =", n,",","prob=",unique(data_histogram$Prob_1)))
}
