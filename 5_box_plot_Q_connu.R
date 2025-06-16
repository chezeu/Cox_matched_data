
#######  alpha varied #########################
setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Methods")

source("3_scenarios_Q_connu.R")
#scenarios_boxplot=scenarios[scenarios$nA==100,] 
#scenarios_boxplot=scenarios[scenarios$nA==500,] 
scenarios_boxplot=scenarios[(scenarios$nA==1000 & scenarios$nB==2000)&
                            (scenarios$censor==2.948469 & scenarios$v_matrix==3 ),] 

varies="alpha"
fixed= "nA"

library(ggplot2)
library(gridExtra)
library(ggpubr)
#theme_set(theme_bw() + theme(legend.position = "top"))
theme_set(theme(legend.position = "top"))

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results")
boxplot_beta0 = NULL
# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  alpha=scenarios_boxplot[i,4]
  censor=scenarios_boxplot[i,5]
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results/Results_Q_connu/","nsim=",nsim,
              "_nA=",nA,"_alpha=",alpha,"_censor=",censor,".Rdata"))
  results_sample = results_sample
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(nA,times=nsim),
                             rep( alpha , times=nsim ),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(nA,times=nsim),
                              rep(alpha, times=nsim ),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(nA,times=nsim),
                               rep(alpha,times=nsim ),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(nA,times=nsim),
                               rep(alpha,times=nsim ),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
  
  boxplot_beta0_Q_sum = cbind(rep(5,times=nsim),
                               rep(nA,times=nsim),
                               rep(alpha,times=nsim ),
                               results_sample$coef_Q_sum_s1,
                               results_sample$coef_Q_sum_s2,
                               results_sample$coef_Q_sum_s3)
  
  boxplot_beta0_estimate = cbind(rep(6,times=nsim),
                                 rep(nA,times=nsim),
                                 rep(alpha,times=nsim ),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,boxplot_beta0_w_sum1,
                        boxplot_beta0_w_sum2,boxplot_beta0_Q_sum,boxplot_beta0_estimate )
  
  
}

#f=quote(n[B])
#titlename = paste0("n_B","=",scenarios_boxplot[,2],
#   ", ","n_A=",scenarios_boxplot[,3])

colnames(boxplot_beta0) = c("Estimator",fixed,varies,"beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:6,labels = c("Theorical"," Naive", "weighted average", "Max weighted average","P_weighted average", "EM"))
boxplot_beta0[,2] = factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = factor(boxplot_beta0[,3] )

boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )

#order=unique(boxplot_beta0[,2])
#xlab= "Proportion of individual with true match"
# Boxplots for beta0

beta0_1_alpha = ggplot(boxplot_beta0,aes(x=alpha,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab( quote(alpha) )+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")
#beta0_1 
beta0_2_alpha = ggplot(boxplot_beta0,aes(x=alpha,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab( quote(alpha) )+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
#beta0_2 

beta0_3_alpha = ggplot(boxplot_beta0,aes(x=alpha,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab( quote(alpha) )+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
#beta0_3

#figure <- ggarrange(beta0_1_alpha, beta0_2_alpha, beta0_3_alpha,
 #                   common.legend = TRUE,
  #                  ncol = 2, nrow = 2)
#figure
#####################################
#



##############     alpha fixed and nA varies ######

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Methods")

source("3_scenarios_Q_connu.R")
#scenarios_boxplot=scenarios[scenarios$alpha==0.8,] 
scenarios_boxplot=scenarios[(scenarios$alpha==0.9& scenarios$v_matrix==3)&
                           (scenarios$censor==2.948469),] 
#scenarios_boxplot=scenarios[scenarios$alpha==1,] 

fixed= "alpha"
varies="nA"
 
library(ggplot2)
library(gridExtra)
library(ggpubr)
theme_set(theme(legend.position = "top"))

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results")
boxplot_beta0 = NULL
# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  alpha=scenarios_boxplot[i,4]
  censor=scenarios_boxplot[i,5]
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results/Results_Q_connu/","nsim=",nsim,
              "_nA=",nA,"_alpha=",alpha,"_censor=",censor,".Rdata"))
  results_sample = results_sample
  #size= paste0("nA=", nA)
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(alpha,times=nsim),
                             rep(nA,times=nsim ),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(alpha,times=nsim),
                              rep( nA, times=nsim ),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(alpha,times=nsim),
                               rep(nA, times=nsim ),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(alpha,times=nsim),
                               rep(nA, times=nsim ),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
  
  boxplot_beta0_Q_sum = cbind(rep(5,times=nsim),
                              rep(alpha,times=nsim),
                              rep(nA,times=nsim ),
                              results_sample$coef_Q_sum_s1,
                              results_sample$coef_Q_sum_s2,
                              results_sample$coef_Q_sum_s3)
  
  boxplot_beta0_estimate = cbind(rep(6,times=nsim),
                                 rep(alpha,times=nsim),
                                 rep( nA, times=nsim ),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)

  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,boxplot_beta0_w_sum1,
                        boxplot_beta0_w_sum2, boxplot_beta0_Q_sum, boxplot_beta0_estimate )
  
  
}
#titlename =paste ("alpha=", 0.9)

colnames(boxplot_beta0) = c("Estimator",fixed,varies, "beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
       levels = 1:6,labels = c("Theorical", "Naive", "weighted average", "Max weighted average","P_weighted average", "EM"))
boxplot_beta0[,3] = factor(boxplot_beta0[,3] )
boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )

#order=unique(boxplot_beta0[,2])
# Boxplots for beta0

beta0_1_nA = ggplot(boxplot_beta0,aes(x=nA,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")
#+  ggtitle(titlename)
#beta0_1 
#+ scale_x_discrete(limits=c(order[1],order[2],order[3]))

beta0_2_nA = ggplot(boxplot_beta0,aes(x=nA,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
#beta0_2
#+ scale_x_discrete(limits=c(order[1],order[2],order[3]))

beta0_3_nA = ggplot(boxplot_beta0,aes(x=nA,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
#beta0_3
#+ scale_x_discrete(limits=c(order[1],order[2],order[3]))

#figure <- ggarrange(beta0_1_nA, beta0_2_nA, beta0_3_nA,
 #                  common.legend = TRUE,
  #                 ncol = 2, nrow = 2)
#figure
#####################################
### Histogram


##############    nA fixed and censor varies ######


setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Methods")

source("3_scenarios_Q_connu.R")
#scenarios_boxplot=scenarios[scenarios$alpha==0.8,] 
scenarios_boxplot = scenarios[(scenarios$alpha==0.9& scenarios$v_matrix==3) &
                               (scenarios$nA==1000),] 
# order 
scenarios_boxplot= rbind(scenarios_boxplot[3,],
                         scenarios_boxplot[1,],
                         scenarios_boxplot[2,])
#scenarios_boxplot=scenarios[scenarios$alpha==1,] 

fixed= "nA"
varies="censor"

library(ggplot2)
library(gridExtra)
library(ggpubr)

setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results")
boxplot_beta0 = NULL
Censor_sample = c(20,30,40) # percentage of censor: see in scenario
# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  alpha=scenarios_boxplot[i,4]
  censor=scenarios_boxplot[i,5]
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results/Results_Q_connu/","nsim=",nsim,
              "_nA=",nA,"_alpha=",alpha,"_censor=",censor,".Rdata"))
  results_sample = results_sample
  #size= paste0("alpha=", nA)
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                       rep(Censor_sample[i],times=nsim),
                             rep(   nA , times=nsim ),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(Censor_sample[i],times=nsim),
                              rep(  nA, times=nsim ),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(Censor_sample[i],times=nsim),
                               rep( nA , times=nsim ),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(Censor_sample[i],times=nsim),
                               rep( nA, times=nsim ),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
 
 boxplot_beta0_Q_sum = cbind(rep(5,times=nsim),
                              rep(Censor_sample[i],times=nsim),
                              rep( nA, times=nsim),
                              results_sample$coef_Q_sum_s1,
                              results_sample$coef_Q_sum_s2,
                              results_sample$coef_Q_sum_s3)
  
  boxplot_beta0_estimate = cbind(rep(6,times=nsim),
                                 rep(Censor_sample[i],times=nsim),
                                 rep( nA, times=nsim ),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,boxplot_beta0_w_sum1,
                        boxplot_beta0_w_sum2, boxplot_beta0_Q_sum,boxplot_beta0_estimate )
  
  
}
#titlename =paste ("alpha=", 0.9)

colnames(boxplot_beta0) = c("Estimator",varies,fixed, "beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:6,labels = c("Theorical", "Naive", "weighted average", "Max weighted average", "P_weighted average", "EM"))
boxplot_beta0[,2] = factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = factor(boxplot_beta0[,3] )
boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )

#order=unique(boxplot_beta0[,2]) 
# Boxplots for beta0

beta0_1_censor = ggplot(boxplot_beta0,aes(x=censor,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(nu))+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")
#+  ggtitle(titlename)
#beta0_1
#+ scale_x_discrete(limits=c(order[3],order[1],order[2]))

beta0_2_censor = ggplot(boxplot_beta0,aes(x=censor,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(nu))+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
#beta0_2
#+ scale_x_discrete(limits=c(order[1],order[2],order[3]))

beta0_3_censor = ggplot(boxplot_beta0,aes(x=censor,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(nu))+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
#beta0_3
#+ scale_x_discrete(limits=c(order[1],order[2],order[3]))

#figure <- ggarrange(beta0_1_censor, beta0_2_censor, beta0_3_censor,
 #                common.legend = TRUE,
  #              ncol = 2, nrow = 2)
#figure

##############    nA fixed and v_matrix varies ######


setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Methods")

source("3_scenarios_Q_connu.R")
scenarios_boxplot = scenarios[(scenarios$alpha==0.9& scenarios$censor== 2.948469) &
                           (scenarios$nA==1000 & scenarios$nB==2000),] 

fixed= "nA"
varies="v_matrix"

library(ggplot2)
library(gridExtra)
library(ggpubr)

#setwd("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results")
boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  alpha=scenarios_boxplot[i,4]
  censor=scenarios_boxplot[i,5]
  v_matrix = scenarios_boxplot[i,6]
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox-matched-data/Results/Results_Q_connu/","nsim=",nsim,
              "_nA=",nA,"_nB=",nB,"_alpha=",alpha,"_censor=",censor,"_v_matrix=",v_matrix,".Rdata"))
  results_sample = results_sample
  #size= paste0("alpha=", nA)
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                             rep(v_matrix,times=nsim),
                             rep(   nB , times=nsim ),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(v_matrix,times=nsim),
                              rep(  nB, times=nsim ),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(v_matrix,times=nsim),
                               rep( nB , times=nsim ),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(v_matrix,times=nsim),
                               rep( nB, times=nsim ),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
  
  boxplot_beta0_Q_sum = cbind(rep(5,times=nsim),
                              rep(v_matrix,times=nsim),
                              rep( nB, times=nsim),
                              results_sample$coef_Q_sum_s1,
                              results_sample$coef_Q_sum_s2,
                              results_sample$coef_Q_sum_s3)
  
  boxplot_beta0_estimate = cbind(rep(6,times=nsim),
                                 rep(v_matrix,times=nsim),
                                 rep( nB, times=nsim ),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,boxplot_beta0_w_sum1,
                        boxplot_beta0_w_sum2, boxplot_beta0_Q_sum,boxplot_beta0_estimate )
  
  
}
#titlename =paste ("alpha=", 0.9)

colnames(boxplot_beta0) = c("Estimator",varies,"nB", "beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:6,labels = c("Theorical", "Naive", "weighted average", "Max weighted average", "P_weighted average", "EM"))
boxplot_beta0[,2] = factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = factor(boxplot_beta0[,3] )
boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )

#order=unique(boxplot_beta0[,2]) 
# Boxplots for beta0

beta0_1_R = ggplot(boxplot_beta0,aes(x=v_matrix,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(n[B]^"*"))+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")
#+  ggtitle(titlename)
#beta0_1_R
#+ scale_x_discrete(limits=c(order[3],order[1],order[2]))

beta0_2_R = ggplot(boxplot_beta0,aes(x=v_matrix,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(n[B]^"*"))+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
#beta0_2_R
#+ scale_x_discrete(limits=c(order[1],order[2],order[3]))

beta0_3_R = ggplot(boxplot_beta0,aes(x=v_matrix,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(n[B]^"*"))+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
#beta0_3_R

#figure <- ggarrange(beta0_1_R, beta0_2_R, beta0_3_R,
#                    common.legend = TRUE,
#                   ncol = 2, nrow = 2)
#figure

figure <- ggarrange(beta0_1_alpha, beta0_2_alpha, beta0_3_alpha,
                    beta0_1_nA, beta0_2_nA, beta0_3_nA,
                    beta0_1_censor, beta0_2_censor, beta0_3_censor,
                    beta0_1_R, beta0_2_R, beta0_3_R,
                    common.legend = TRUE,
                    ncol = 3, nrow = 4)

figure
#####################################


data_histogram = boxplot_beta0[which(boxplot_beta0$Estimator==unique(boxplot_beta0$Estimator)[5]),]
nA_hist = unique(data_histogram$nA)

for (i in 1:length(nA_hist)) {
  n = nA_hist[i]
  beta_hist =data_histogram[which(data_histogram$nA==n),]
  hist(beta_hist$beta0.1, xlab= quote(hat(beta)[1]),main = paste("nA =", n,",","prob=",unique(data_histogram$Prob_1)))
  hist(beta_hist$beta0.2, xlab= quote(hat(beta)[2]),main = paste("nA =", n,",","prob=",unique(data_histogram$Prob_1)))
}
