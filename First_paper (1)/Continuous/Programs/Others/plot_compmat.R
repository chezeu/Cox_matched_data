setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Programs")

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(EnvStats)

source("1_generate_data.R")
source("2_compare.R")


nA = 500
nB = 200
K = 1
date_max = 100

error = c(0.2,0.4,0.6)
mu = c(1,2,3)

plotlist1 = list()
plotlist2 = list()
plotlist3 = list()
plotlist4 = list()
for (i in 1:length(error)){
  for (j in 1:length(mu)){
    
    ## EXP
    #data = generate_data_2(nA, nB, K, date_max, error = error[i], mu = mu[j])
    
    ## Uniform
    data = generate_data_4_norm(nA, nB, K, date_max, error = error[i], mu = mu[j])
    datA = data$dataA
    datB = data$dataB
    
    comp_mat = compare_1(datA, datB, K)
    
    M = data.frame(value = comp_mat[comp_mat[,2]==1,1])
    U = data.frame(value = comp_mat[comp_mat[,2]==0,1])
    M$status = "Matches"
    U$status = "Non-matches"
    
    df = rbind(M,U)
    
    p1 <- df %>%
      ggplot( aes(x=value, fill=status, color = status)) +
      geom_histogram( alpha=0.3, position = 'identity', binwidth = 1)+
      ggtitle(paste("error =",error[i],"mu = ",mu[j]))+
      ylim(0,6000)
      #labs(fill="") 
    
    p2 <- ggplot(df)+
      geom_histogram( aes(x=value, y=..density.., fill=status, color = status), alpha=0.3, position = 'identity', binwidth = 1)+
      ggtitle(paste("error =",error[i],"mu = ",mu[j]))+
      ylim(0,0.85)
    
    
    plotlist1[[j+(i-1)*length(mu)]] = p1
    plotlist2[[j+(i-1)*length(mu)]] = p2
    
    s2 = var(M[,1])
    m = mean(M[,1])
    lambda = (s2+m^2)/m-1
    
    p3 <- ggplot(M)+
      geom_histogram( aes(x=value, y=..density.., fill=status, color = status), alpha=0.3, position = 'identity', binwidth = 1)+
      stat_function(geom = "line",xlim=range(M[,1]), n = max(M[,1])+1, fun = dpois, args = list(lambda = mean(M[,1])), show.legend = TRUE) +
      stat_function(geom = "line",xlim=c(0,max(M[,1])), n = max(M[,1])+1, fun = dpois, args = list(lambda = lambda), show.legend = TRUE, colour = "blue") +
      ggtitle(paste("error =",error[i],"mu = ",mu[j]))+
      ylim(0,0.85)
    p3
    
    s2 = var(U[,1])
    m = mean(U[,1])
    lambda = (s2+m^2)/m-1
    
    p3U <- ggplot(U)+
      geom_histogram( aes(x=value, y=..density.., fill=status, color = status), alpha=0.3, position = 'identity', binwidth = 1)+
      stat_function(geom = "line",xlim=range(U[,1]), n = max(U[,1])+1, fun = dpois, args = list(lambda = mean(U[,1])), show.legend = TRUE) +
      stat_function(geom = "line",xlim=c(0,max(U[,1])), n = max(U[,1])+1, fun = dpois, args = list(lambda = lambda), show.legend = TRUE, colour = "blue") +
      ggtitle(paste("error =",error[i],"mu = ",mu[j]))+
      ylim(0,0.3)
    p3U
    
    
    fitM =  egamma(M[,1]+1e-15)
    shape_M = fitM$parameters[1]
    scale_M = fitM$parameters[2]
    
    fitU =  egamma(U[,1]+1e-15)
    shape_U = fitU$parameters[1]
    scale_U = fitU$parameters[2]
    
    fitM_zig =  egamma(M[M[,1]>0,1])
    shape_M_zig = fitM_zig$parameters[1]
    scale_M_zig = fitM_zig$parameters[2]
    
    fitU_zig =  egamma(U[U[,1]>0,1])
    shape_U_zig = fitU_zig$parameters[1]
    scale_U_zig = fitU_zig$parameters[2]
    
    p4 <- ggplot(M)+
      geom_histogram( aes(x=value, y=..density.., fill=status, color = status), alpha=0.3, position = 'identity', binwidth = 1)+
      stat_function(geom = "line",n = max(M[,1])+1, fun = dgamma, args = list(shape = shape_M, scale = scale_M)) +
      stat_function(geom = "line",xlim = c(1,max(M[,1])), n = max(M[,1]), fun = dgamma, args = list(shape = shape_M_zig, scale = scale_M_zig), colour = "blue") +
      ggtitle(paste("error =",error[i],"mu = ",mu[j]))+
      ylim(0,0.85)
    p4
    
    plotlist3[[j+(i-1)*length(mu)]] = p3
    plotlist4[[j+(i-1)*length(mu)]] = p4
  }
}

fig1 = ggarrange(plotlist=plotlist1, nrow = length(error), ncol = length(mu), common.legend = TRUE)
ggsave("hist_distance2_K1_exp_exp.pdf", width = 30, height = 30, units = "cm")


fig2 = ggarrange(plotlist=plotlist2, nrow = length(error), ncol = length(mu), common.legend = TRUE)
ggsave("dens_distance2_K1_exp_exp.pdf", width = 30, height = 30, units = "cm")

fig3 = ggarrange(plotlist=plotlist3, nrow = length(error), ncol = length(mu), common.legend = TRUE)
ggsave("poisson_distance2_K1_exp_exp.pdf", width = 30, height = 30, units = "cm")

fig4 = ggarrange(plotlist=plotlist4, nrow = length(error), ncol = length(mu), common.legend = TRUE)
ggsave("gamma_distance2_K1_exp_exp.pdf", width = 30, height = 30, units = "cm")
##################################################################-

nA = 500
nB = 200
mu = 3
date_max = 100

error = c(0.1,0.3,0.5)

K = c(1,3,5,10)


plotlist1 = list()
plotlist2 = list()
for (j in 1:length(K)){
  for (i in 1:length(error)){
    
    data = generate_data_2(nA, nB, K[j], date_max, error = error[i], mu = mu)
    
    datA = data$dataA
    datB = data$dataB
    
    comp_mat = compare_1(datA, datB, K[j])
    
    M = data.frame(value = comp_mat[comp_mat[,2]==1,1])
    U = data.frame(value = comp_mat[comp_mat[,2]==0,1])
    M$status = "Matches"
    U$status = "Non-matches"
    
    df = rbind(M,U)
    
    p1 <- df %>%
      ggplot( aes(x=value, fill=status, color = status)) +
      geom_histogram( alpha=0.3, position = 'identity', binwidth = 2)+
      ggtitle(paste("error =",error[i], "K =", K[j])) +
      ylim(0,12000)
    #labs(fill="") 
    
    p2 <- ggplot(df)+
      geom_histogram( aes(x=value, y=..density.., fill=status, color = status), alpha=0.3, position = 'identity', binwidth = 2)+
      ggtitle(paste("error =",error[i], "K =", K[j]))+
      ylim(0,0.6)
    
    
    plotlist1[[i+(j-1)*length(error)]] = p1
    plotlist2[[i+(j-1)*length(error)]] = p2
    
  }
}




fig1 = ggarrange(plotlist=plotlist1, nrow = length(K), ncol = length(error), common.legend = TRUE)
ggsave("distance_errorK.pdf", width = 30, height = 40, units = "cm")


fig2 = ggarrange(plotlist=plotlist2, nrow = length(K), ncol = length(error), common.legend = TRUE)
ggsave("distance_errorK_density.pdf", width = 30, height = 40, units = "cm")




########################################
nA = 10
nB = 5
K = 3
date_max = 100

error = 0.5
mu = 3

data = generate_data_2(nA,nB,K,date_max, error, mu)
datA = data$dataA
datB = data$dataB

comp1 = compare_1(datA, datB, K)
comp1

comp2 = compare_binary(datA,datB,K)
comp2

comp_mix = cbind(comp1[,1],comp2)
M = comp_mix[comp_mix[,K+2] ==1,]
M
