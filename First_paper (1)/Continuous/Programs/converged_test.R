library(simsalapar)
library(doParallel)
library(clue)
library(matrixStats)
library(klaR)
library(ludic)

setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Binary\\Programs")
nA = 500
nB = 200
K = 25
error = 0.04

source("1_generate_data.R")
source("2_compare.R")
source("3_EM.R")
source("4_methods.R")

data = generate_data_test(nA = nA, nB = nB, K = K, error = error)

datA = data$dataA
datB = data$dataB
prev = data$prev

g <- true_prob(datA, datB, K, prev, error)
g_true = g$g_true

indM = g$indM
indU = g$indU

m2.true = g$m2
m1.true = g$m1
u2.true = g$u2
u1.true = g$u1

p.true = 1/nA


comp_mat <- compare3(datA, datB, K=K)
fit = EM3_test(comp_mat, datA, datB, K, tol=1e-6, maxits = 500)


iter = fit$iter
iteration= 1:(iter +1)

loglikelihood = fit$loglikelihood
plot(iteration,loglikelihood,main=expression("Evolution of log-likelihood"))

p = fit$p
plot(iteration, p, main=expression("Evolution of p"))
points(iter+2,p.true, col = "red", pch = 19)

m2.1 = fit$m2[,1]
plot(iteration, m2.1, main=expression("Evolution of m"[2]^1))
points(iter+2,m2.true[1], col = "red", pch = 19)

m1.1 = fit$m1[,1]
plot(iteration, m1.1, main=expression("Evolution of m"[1]^1))
points(iter+2,m1.true[1], col = "red", pch = 19)

u2.1 = fit$u2[,1]
plot(iteration, u2.1, main=expression("Evolution of u"[2]^1))
points(iter+2,u2.true[1], col = "red", pch = 19)

u1.1 = fit$u1[,1]
plot(iteration, u1.1, main=expression("Evolution of u"[1]^1))
points(iter+2,u1.true[1], col = "red", pch = 19)

fit2 = EM3_test(comp_mat, datA, datB, K, tol=1e-6, maxits = 1000)
