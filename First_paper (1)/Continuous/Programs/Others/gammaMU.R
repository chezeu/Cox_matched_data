library(simsalapar)
library(doParallel)
library(clue)
library(matrixStats)
library(EnvStats)
library(ggplot2)
library(ggpubr)

setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Programs")

source("1c_generate_data.R")
source("2c_compare.R")
source("3c_EM.R")
source("4c_methods.R")

nA = 500
nB = 200
K = 1
lambdaK = 0.02
error = 0.2
lambdaE = 1/2
round = TRUE
data = generate_data_exp_exp(nA = nA, nB=nB,K=K, lambdaK = lambdaK, error = error, lambdaE = lambdaE, round = round)
datA = data$dataA
datB = data$dataB

comp_mat = compare_abs(datA, datB, K)
M = comp_mat[comp_mat[,2]==1,1]
U = comp_mat[comp_mat[,2]==0,1]
M1 = M[M>0]
U1 = U[U>0]


###################### M
fitM = egamma(M1)
alphaM = fitM$parameters[1]
betaM = fitM$parameters[2]

estiM = dgamma(M1, shape = alphaM, scale = betaM)
df.M1 = data.frame(distance = M1, esti = estiM)

figM <- ggplot(data = df.M1, aes(x = distance))+
          geom_histogram(aes(y=..density..), binwidth = 1,color="darkblue", fill="lightblue")+
          stat_function(fun = dgamma, args = list(shape = alphaM, scale = betaM))

figM
