setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Programs")

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(clue)
library(matrixStats)
library(ludic)
library(ggpubr)

source("1_generate_data.R")
source("2_compare.R")
source("3_EM_gamma.R")


nA = 500
nB = 200
K = 5
date_max = 100

error = 0.6
mu = 2

    
data = generate_data_2(nA, nB, K, date_max, error = error, mu = mu)
    
datA = data$dataA
datB = data$dataB

comp_mat = compare_1(datA, datB, K)

adding = 1e-15

X <- comp_mat[,1]
X <- X + adding

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
  ggtitle(paste("Density"))+
  ylim(0,0.6)

fig2

ggarrange(fig1, fig2, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("distance.pdf", width = 20, height = 15, units = "cm")
######################### True
X <- comp_mat[,1]
X <- X + adding

nA = nrow(datA)
nB = nrow(datB)
N = length(X)


M = X[comp_mat[,2]==1]
U = X[comp_mat[,2]==0]

fitM =  egamma(M)
shape_M = fitM$parameters[1]
scale_M = fitM$parameters[2]

fitU =  egamma(U)
shape_U = fitU$parameters[1]
scale_U = fitU$parameters[2]

p_true = length(M)/length(X)

alpha_true = c(shape_M, shape_U)
beta_true = c(scale_M, scale_U)

alpha_true
beta_true
p_true

dfM_true = data.frame(x = X, y = dgamma(X, shape = alpha_true[1], scale = beta_true[1]))
dfU_true = data.frame(x = X, y = dgamma(X, shape = alpha_true[2], scale = beta_true[2]))
dfmix_true = data.frame(x = X, y =p_true*dgamma(X, shape = alpha_true[1], scale = beta_true[1])
                                + (1-p_true)*dgamma(X, shape = alpha_true[2], scale = beta_true[2]))



probM <- dgamma(X, shape = shape_M, scale = scale_M)
probU <- dgamma(X, shape = shape_U, scale = scale_U)

g <- (p_true*probM)/(p_true*probM+(1-p_true)*probU)
g[g<0]=0
Gmat = matrix(g, nrow = nB, byrow = TRUE)
opti_cont <- solve_LSAP(Gmat, maximum = TRUE)
predict_cont = cbind(seq_along(opti_cont), opti_cont)
# Percentage of correct link
TPR_gamma_true_11 = sum(predict_cont[,2] == datB[,K+1])/nB
TPR_gamma_true_11
# 


fig3 = fig2 + geom_line(data = dfM_true, aes(x = x, y = y), color="red", size = 1.5) +
  geom_line(data = dfU_true, aes(x = x, y = y), color="green", size = 1.5) +
  geom_line(data = dfmix_true, aes(x = x, y = y), color="black")+
  ggtitle(paste("True fit: ","aM =", round(shape_M,2),", bM = ", round(scale_M,2),", p=", round(p_true,4),
                "\n aU =", round(shape_U,2),", bU = ", round(scale_U,2), ", TPR=", round(TPR_gamma_true_11,2), ", max_y=", round(max(dfM_true$y),2)))
  

fig3
######################## EM estimated
fit = EM_gamma_pfix(comp_mat,datA,datB, K, tol = 1e-6, maxits = 500)

alpha_est = fit$alpha
beta_est = fit$beta
p_est = fit$p

alpha_est
beta_est
p_est = 1/nA

dfM_est = data.frame(x = X, y = dgamma(X, shape = alpha_est[1], scale = beta_est[1]))
dfU_est = data.frame(x = X, y = dgamma(X, shape = alpha_est[2], scale = beta_est[2]))
dfmix_est = data.frame(x = X, y =p_est*dgamma(X, shape = alpha_est[1], scale = beta_est[1])
                        + (1-p_est)*dgamma(X, shape = alpha_est[2], scale = beta_est[2]))



# fig5 = fig2 + geom_point(data = dfM_est, aes(x = x, y = y), color="red" ) +
#   ylim(0,400)+
#   ggtitle(paste("Full Matches fit: ","error =",error,"mu = ",mu))
# fig5

g = fit$g
g[g<0]=0
Gmat = matrix(g, nrow = nB, byrow = TRUE)
opti_cont <- solve_LSAP(Gmat, maximum = TRUE)
predict_cont = cbind(seq_along(opti_cont), opti_cont)
# Percentage of correct link
TPR_gamma_11 = sum(predict_cont[,2] == datB[,K+1])/nB
TPR_gamma_11


fig4 = fig2 + geom_line(data = dfM_est, aes(x = x, y = y), color="red", size = 1.5) +
  geom_line(data = dfU_est, aes(x = x, y = y), color="green", size = 1.5) +
  geom_line(data = dfmix_est, aes(x = x, y = y), color="black") +
  ggtitle(paste("EM fit: ","aM =", round(alpha_est[1],2),", bM = ", round(beta_est[1],4),", p=", round(p_est,4),
                "\n aU =", round(alpha_est[2],2),", bU = ", round(beta_est[2],2),", TPR=", round(TPR_gamma_11,2), ", max_y=", round(max(dfM_est$y),2)))
fig4
# 


# ####################### noEM estimated
fit_noEM = noEM_gamma(comp_mat,datA,datB, K, tol = 1e-5, maxits = 300)

alpha_noEM = fit_noEM$alpha
beta_noEM = fit_noEM$beta
p_noEM = fit_noEM$p

alpha_noEM
beta_noEM
p_noEM


dfM_noEM = data.frame(x = X, y = dgamma(X, shape = alpha_noEM[1], scale = beta_noEM[1]))
dfU_noEM = data.frame(x = X, y = dgamma(X, shape = alpha_noEM[2], scale = beta_noEM[2]))
dfmix_noEM = data.frame(x = X, y =p_noEM*dgamma(X, shape = alpha_noEM[1], scale = beta_noEM[1])
                        + (1-p_noEM)*dgamma(X, shape = alpha_noEM[2], scale = beta_noEM[2]))



g = fit_noEM$g
g[g<0]=0
Gmat = matrix(g, nrow = nB, byrow = TRUE)
opti_cont <- solve_LSAP(Gmat, maximum = TRUE)
predict_cont = cbind(seq_along(opti_cont), opti_cont)
# Percentage of correct link
TPR_gammanoEM_11 = sum(predict_cont[,2] == datB[,K+1])/nB
TPR_gammanoEM_11

M = sort(X)[1:nB]
aa = comp_mat[comp_mat[,1]<= (max(M)+0.01),]

pM = sum(aa[,2])/nB
fig6 = fig2 + geom_line(data = dfM_noEM, aes(x = x, y = y), color="red",size = 1.5) +
  geom_line(data = dfU_noEM, aes(x = x, y = y), color="green", size = 1.5) +
  geom_line(data = dfmix_noEM, aes(x = x, y = y), color="black") +
  ggtitle(paste("noEM fit: ","aM =", round(alpha_noEM[1],2),", bM = ", round(beta_noEM[1],4),", p=", round(p_noEM,4),
                "\n aU =", round(alpha_noEM[2],2),", bU = ", round(beta_noEM[2],2),", TPR=", round(TPR_gammanoEM_11,2),", max_y=", round(max(dfM_noEM$y),2), ", pStart =",round(pM,2)))
fig6

ggarrange(fig1,fig2,fig3,fig4, fig6, nrow = 3, ncol = 2, common.legend = TRUE)

ggsave("fit_1e15.pdf", width = 20, height = 30)


