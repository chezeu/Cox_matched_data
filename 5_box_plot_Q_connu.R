
#######  alpha varied #########################
setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Methods")

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

theme_set(theme(legend.position = "top"))

setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Results")
boxplot_beta0 = NULL
# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  alpha=scenarios_boxplot[i,4]
  censor=scenarios_boxplot[i,5]
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Results/Results_Q_connu/","nsim=",nsim,
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

 
colnames(boxplot_beta0) = c("Estimator",fixed,varies,"beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:6,labels = c("Theoretical", "Naive", "Weighted average", "Max weighted average", "Weighted PS", "Full likelihood"))
boxplot_beta0[,2] = factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = factor(boxplot_beta0[,3] )

boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )
 
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
 
##############     alpha fixed and nA varies ######

setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Methods")

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

setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Results")
boxplot_beta0 = NULL
# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  alpha=scenarios_boxplot[i,4]
  censor=scenarios_boxplot[i,5]
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Results/Results_Q_connu/","nsim=",nsim,
              "_nA=",nA,"_alpha=",alpha,"_censor=",censor,".Rdata"))
  results_sample = results_sample
   
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
 
colnames(boxplot_beta0) = c("Estimator",fixed,varies, "beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
       levels = 1:6,labels = c("Theoretical", "Naive", "Weighted average", "Max weighted average", "Weighted PS", "Full likelihood"))
boxplot_beta0[,3] = factor(boxplot_beta0[,3] )
boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )
 
beta0_1_nA = ggplot(boxplot_beta0,aes(x=nA,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")
 
beta0_2_nA = ggplot(boxplot_beta0,aes(x=nA,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
 
beta0_3_nA = ggplot(boxplot_beta0,aes(x=nA,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(n[A]))+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
#beta0_3
 
##############    nA fixed and censor varies ######


setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Methods")

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

setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Results")
boxplot_beta0 = NULL
Censor_sample = c(20,30,40) # percentage of censor: see in scenario
# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  alpha=scenarios_boxplot[i,4]
   censor=scenarios_boxplot[i,5]
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Results/Results_Q_connu/","nsim=",nsim,
              "_nA=",nA,"_alpha=",alpha,"_censor=",censor,".Rdata"))
  results_sample = results_sample
   censor= Censor_sample[i]
  
  boxplot_beta0_true = cbind(rep(1,times=nsim),
                       rep(censor,times=nsim),
                             rep(   nA , times=nsim ),
                             results_sample$coef_true_s1,
                             results_sample$coef_true_s2,
                             results_sample$coef_true_s3)
  
  boxplot_beta0_naive = cbind(rep(2,times=nsim),
                              rep(censor,times=nsim),
                              rep(  nA, times=nsim ),
                              results_sample$coef_naive_s1,
                              results_sample$coef_naive_s2,
                              results_sample$coef_naive_s3)
  
  boxplot_beta0_w_sum1 = cbind(rep(3,times=nsim),
                               rep(censor,times=nsim),
                               rep( nA , times=nsim ),
                               results_sample$coef_w_sum1_s1,
                               results_sample$coef_w_sum1_s2,
                               results_sample$coef_w_sum1_s3)
  
  boxplot_beta0_w_sum2 = cbind(rep(4,times=nsim),
                               rep(censor,times=nsim),
                               rep( nA, times=nsim ),
                               results_sample$coef_w_sum2_s1,
                               results_sample$coef_w_sum2_s2,
                               results_sample$coef_w_sum2_s3)
 
 boxplot_beta0_Q_sum = cbind(rep(5,times=nsim),
                              rep(censor,times=nsim),
                              rep( nA, times=nsim),
                              results_sample$coef_Q_sum_s1,
                              results_sample$coef_Q_sum_s2,
                              results_sample$coef_Q_sum_s3)
  
  boxplot_beta0_estimate = cbind(rep(6,times=nsim),
                                 rep(censor,times=nsim),
                                 rep( nA, times=nsim ),
                                 results_sample$coef_estimate_s1,
                                 results_sample$coef_estimate_s2,
                                 results_sample$coef_estimate_s3)
  
  boxplot_beta0 = rbind(boxplot_beta0,boxplot_beta0_true,boxplot_beta0_naive,boxplot_beta0_w_sum1,
                        boxplot_beta0_w_sum2, boxplot_beta0_Q_sum,boxplot_beta0_estimate )
  
  
}

colnames(boxplot_beta0) = c("Estimator",varies,fixed, "beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:6,labels = c("Theoretical", "Naive", "Weighted average", "Max weighted average", "Weighted PS", "Full likelihood"))
boxplot_beta0[,2] = factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = factor(boxplot_beta0[,3] )
boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )
 
beta0_1_censor = ggplot(boxplot_beta0,aes(x=censor,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(nu))+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")
 
beta0_2_censor = ggplot(boxplot_beta0,aes(x=censor,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(nu))+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
#beta0_2
 
beta0_3_censor = ggplot(boxplot_beta0,aes(x=censor,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(nu))+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
 
##############    nA fixed and v_matrix varies ######


setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Methods")

source("3_scenarios_Q_connu.R")
scenarios_boxplot = scenarios[(scenarios$alpha==0.9& scenarios$censor== 2.948469) &
                           (scenarios$nA==1000 & scenarios$nB==2000),] 

fixed= "nA"
varies="v_matrix"

library(ggplot2)
library(gridExtra)
library(ggpubr)
 
boxplot_beta0 = NULL

# size fix
for (i in (1:nrow(scenarios_boxplot))){
  
  nsim = scenarios_boxplot[i,1]
  nB = scenarios_boxplot[i,2]
  nA = scenarios_boxplot[i,3]
  alpha=scenarios_boxplot[i,4]
  censor=scenarios_boxplot[i,5]
  v_matrix = scenarios_boxplot[i,6]
  
  load(paste0("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Results/Results_Q_connu/","nsim=",nsim,
              "_nA=",nA,"_nB=",nB,"_alpha=",alpha,"_censor=",censor,"_v_matrix=",v_matrix,".Rdata"))
  results_sample = results_sample
   
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
 
colnames(boxplot_beta0) = c("Estimator",varies,"nB", "beta0.1","beta0.2","beta0.3")

boxplot_beta0 = data.frame(boxplot_beta0)

boxplot_beta0[,1] = factor(boxplot_beta0[,1],
                           levels = 1:6,labels = c("Theoretical", "Naive", "Weighted average", "Max weighted average", "Weighted PS", "Full likelihood"))
boxplot_beta0[,2] = factor(boxplot_beta0[,2] )
boxplot_beta0[,3] = factor(boxplot_beta0[,3] )
boxplot_beta0[,4] = as.numeric(boxplot_beta0[,4] )
boxplot_beta0[,5] = as.numeric(boxplot_beta0[,5] )
boxplot_beta0[,6] = as.numeric(boxplot_beta0[,6] )
 
beta0_1_R = ggplot(boxplot_beta0,aes(x=v_matrix,y= beta0.1,fill=Estimator))+
  geom_boxplot()+
  xlab(quote(n[B]^"*"))+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=beta[1],lty=1,col="orange")
 
#beta0_1_R
 
beta0_2_R = ggplot(boxplot_beta0,aes(x=v_matrix,y= beta0.2,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(n[B]^"*"))+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=beta[2],lty=1,col="orange")
#beta0_2_R
 
beta0_3_R = ggplot(boxplot_beta0,aes(x=v_matrix,y= beta0.3,fill=Estimator))+
  geom_boxplot()+
  theme(legend.position = "none")+
  xlab(quote(n[B]^"*"))+
  ylab(quote(hat(beta)[3]))+
  geom_hline(yintercept=beta[3],lty=1,col="orange")
#beta0_3_R

 
figure <- ggarrange(beta0_1_alpha, beta0_2_alpha, beta0_3_alpha,
                    beta0_1_nA, beta0_2_nA, beta0_3_nA,
                    beta0_1_censor, beta0_2_censor, beta0_3_censor,
                    beta0_1_R, beta0_2_R, beta0_3_R,
                    common.legend = TRUE,
                    ncol = 3, nrow = 4)

figure

figure_revue <- figure +
  labs(title = paste0(" Fidelette_CHEZEU ", " - Figure ", 2, " - page 20")) +
  theme(plot.title = element_text(size = 12, hjust = 0.5)) 
fichier <- paste0("Fidelette_Chezeu", "_Figure", 2)
ggsave(paste0(fichier, ".pdf"), figure_revue, device = cairo_pdf,
       width = 28, height = 18, units = "cm")
browseURL(paste0(fichier, ".pdf"))
#####################################

