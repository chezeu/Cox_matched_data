library(mvtnorm)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(clue)
library(matrixStats)
library(ludic)
library(ggpubr)

# Generate data

nA = 500
nB = 200
nM = nB
N = nA*nB

alpha_M = 1.4
beta_M = 4.3
alpha_U = 5.3
beta_U = 7.8
p = 1/nA


M = data.frame(value = rgamma(nM, shape = alpha_M, scale = beta_M))
U = data.frame(value =  rgamma(N-nM, shape = alpha_U, scale = beta_U))
M$status = "1"
U$status = "0"
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

X = df[,1]

#True posterior probabilities
g_true = p*dgamma(X,shape = alpha_M, scale = beta_M)/(p*dgamma(X,shape = alpha_M, scale = beta_M)+(1-p)*dgamma(X,shape = alpha_U, scale = beta_U))

################# Fitting
