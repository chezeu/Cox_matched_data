setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Programs")

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(clue)
library(matrixStats)
library(ludic)
library(ggpubr)
library(Rfast)

source("1_generate_data.R")
source("2_compare.R")
source("3_EM_gamma.R")


nA = 500
nB = 200
K = 2
date_max = 100

error = 0.2
mu = 1


data = generate_data_4(nA, nB, K, date_max, error = error, mu = mu)

datA = data$dataA
datB = data$dataB

comp_mat = compare_1(datA, datB, K)


X = comp_mat[,1]

#X = (X-min(X))/(max(X)-min(X))

M = data.frame(value = X[comp_mat[,2]==1])
U = data.frame(value = X[comp_mat[,2]==0])   
M$status = "Matches"
U$status = "Non-matches"
df = rbind(M,U)


fig1 <- ggplot(df)+
  geom_histogram( aes(x=value, fill=status, color = status), alpha=0.3, position = 'identity', binwidth = 2)+
  ggtitle(paste("Histogram"))

fig1


fig2 <- ggplot(df)+
  geom_histogram( aes(x=value, y=..density.., fill=status, color = status), alpha=0.3, position = 'identity', binwidth = 2)+
  ggtitle(paste("Density"))

fig2


##########################
X = comp_mat[,1]
X = (X-min(X))/(max(X)-min(X))

X[which(X==0)] <- 1e-6
X[which(X==1)] <- 1-(1e-6)

M = X[comp_mat[,2]==1]
U = X[comp_mat[,2]==0]

fitM = beta.mle(M)
alpha_M = fitM$param[1]
beta_M = fitM$param[2]

fitU = beta.mle(U)
alpha_U = fitU$param[1]
beta_U = fitU$param[2]

#p = nB/length(X)
p = length(M)/length(X)
probM <- dbeta(X, shape1 = alpha_M, shape2 = beta_M)
probU <- dbeta(X, shape1 = alpha_U, shape2 = beta_U)
g <- (p*probM)/(p*probM+(1-p)*probU)
g[g<0]=0
Gmat = matrix(g, nrow = nB, byrow = TRUE)
opti_cont <- solve_LSAP(Gmat, maximum = TRUE)
predict_cont = cbind(seq_along(opti_cont), opti_cont)
# Percentage of correct link
TPR_beta_true_11 = sum(predict_cont[,2] == datB[,K+1])/nB
TPR_beta_true_11

dfM_true = data.frame(x = X, y = dbeta(X, shape1 = alpha_M, shape2 = beta_M))
dfU_true = data.frame(x = X, y = dbeta(X, shape1 = alpha_U, shape2 = beta_U))
dfmix_true = data.frame(x = X, y =p*dbeta(X, shape1 = alpha_M, shape2 = beta_M)
                        + (1-p)*dbeta(X, shape1 = alpha_U, shape2 = beta_U))

fig3 = fig2 + geom_line(data = dfM_true, aes(x = x, y = y), color="red", size = 1.5) +
  geom_line(data = dfU_true, aes(x = x, y = y), color="green", size = 1.5) +
  geom_line(data = dfmix_true, aes(x = x, y = y), color="black")+
  ggtitle(paste("True fit: ","aM =", round(alpha_M,2),", bM = ", round(beta_M,2),", p=", round(p,4),
                "\n aU =", round(alpha_U,2),", bU = ", round(beta_U,2), ", TPR=", round(TPR_beta_true_11,2)))


fig3