library(simsalapar)
library(doParallel)
library(clue)
library(matrixStats)
library(klaR)
library(ludic)

setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Binary\\Programs")

source("b1_generate_data.R")
source("b2_compare.R")
source("b3_EM.R")
source("bprc4_methods.R")

prev  = 0.2
K = 40
nA = 500
nB = 200
error = 0.04


prevalence = rep(prev, K)
data = generate_data(nA = nA, nB = nB, K = K, prevalence = prevalence, error = error)

datA = data$dataA
datB = data$dataB
prev = data$prev

g <- true_prob(datA, datB, K, prev, error)
g_true = g$g_true
indM = g$indM
indU = g$indU
comp_mat4 = g$comp_mat4

tol = 1e-6
maxits = 500

comp_mat2 <- compare_binary(datA, datB, K)
fit_FS <- EM_binary(comp_mat2, datA, datB, K,tol = tol, maxits = maxits)
g_FS = fit_FS$g


comp_mat3 <- compare3(datA, datB, K=K)
fit_FS3 = EM3(comp_mat3, datA, datB, K, tol=tol, maxits = maxits)
g_FS3 = fit_FS3$g

## Using EM
fit_FS4 = EM4(comp_mat4, datA, datB, K, tol=tol, maxits = maxits)
g_FS4 = fit_FS4$g

bayes = recordLink(datB[,1:K], datA[,1:K], eps_plus =0.01, eps_minus = 0.01,use_diff = FALSE)
g_bayes = as.vector(t(bayes))

################# plots#########
library(precrec)

labels = comp_mat2[,K+1]
scores <- join_scores(g_FS, g_FS3, g_FS4, g_bayes)
# Explicitly specify model names
msmdat <- mmdata(scores, labels, modnames = c("FS", "FS3", "FS4", "Bayesian"))

# Calculate ROC and Precision-Recall curves for multiple models
mscurves <- evalmod(msmdat)

