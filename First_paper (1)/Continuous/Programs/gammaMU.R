library(simsalapar)
library(doParallel)
library(clue)
library(matrixStats)
library(EnvStats)

setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Programs")

source("c1_generate_data.R")
source("c2_compare.R")
source("c3_EM.R")
source("c4_methods.R")

nA = 500
nB = 200
K = 1
lambdaK = 0.02
error = 0.2
lambdaE = 1/3
round =TRUE

data = generate_data_exp_exp(nA = nA, nB=nB,K=K, lambdaK = lambdaK, error = error, lambdaE = lambdaE, round = round)
#data = generate_data_unif_exp(500,200,1,100,0.2,1/3,round = FALSE)
datA = data$dataA
datB = data$dataB

comp_mat = compare_abs(datA, datB, K)
M = comp_mat[comp_mat[,2]==1,1]
U = comp_mat[comp_mat[,2]==0,1]

M1 = M[M>0]
U1 = U[U>0]

p0M = sum(M==0)/length(M)
p0U = sum(U==0)/length(U)

fitM = egamma(M1)
shapeM = fitM$parameters[1]
scaleM = fitM$parameters[2]

fitU = egamma(U1)
shapeU = fitU$parameters[1]
scaleU = fitU$parameters[2]

x = seq(1,max(M1)+2, by=1 )
y = (1-p0M)*dgamma(x, shape = shapeM, scale = scaleM)

hist(M, breaks=10, freq = FALSE, main = NULL, xlab = expression(gamma))
lines(x, y, col = "red" )
points(0,p0M,col = "red", pch = 19  )

x = seq(2,max(U)+2, length.out = 100)
y = (1-p0U)*dgamma(x, shape = shapeU, scale = scaleU)
hist(U1, breaks = 30, freq = FALSE, main = NULL, xlab = expression(gamma))
lines(x, y, col = "red" )
points(0,p0U, col = "red" , pch = 19 )

library(MASS)

my.mle<-fitdistr(M, densfun="poisson")
x = 0:10
y = dpois(x, lambda = my.mle$estimate)
hist(M, breaks =5, freq = FALSE, main = NULL, xlab = expression(gamma))
lines(x, y, col = "red" )

my.mle<-fitdistr(U1, densfun="poisson")
x = 0:200
y = dpois(x, lambda = my.mle$estimate)
hist(U1, breaks =30, freq = FALSE, main = NULL, xlab = expression(gamma))
lines(x, y, col = "red" )
