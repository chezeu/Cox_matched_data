#  FINAL RESULTS TABLE
 
library(survival)
library(rootSolve)
library(nleqslv)
library(reclin2)
library(Matrix)
setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/Methods")

source("1_data_generate_Q_inconnu.R")
source("2_matrice_Q_inconnu.R")
source("3_scenarios_Q_inconnu.R")
source("4_4_method_EM.R")
source("6_variance.R")
 
i=1
nsim=scenarios[i,1]
K=scenarios[i,2]
nA=scenarios[i,3]
nB= scenarios[i,4]
censor= scenarios[i,5]
prevalence = rep(prevalence_sample,K)

#data generation
mf = Generate_data(K,nA,nB,prevalence,censor,min_prev = 0.01)

surv_data = mf$surv_data
XB = mf$XB
data_true=surv_data
idBA=mf$idBA

databaseA=mf$databaseA
databaseB =mf$ databaseB

datA=mf$datA # matching variables
datB= mf$datB
#####
#record linkage 
comp_mat = compare_binary (datA, datB, K)
linked = Matrix_Q(comp_mat,K,datA, datB)
Q=linked$Q
summary (as.vector(Q))
 
naive_id= linked$naive_id
  

Ts = data_true$Time
event = data_true$delta
XB=data.frame(XB)
 
#lambda0 =  rep(0.1,length(event))
#lambda0 [which(event == 0)] = 0

beta0 = rep(0,p)

# Estimator of Nelson-Aalen
initialiser_lambda0 <- function(Ts, event) {
  
  fit     <- survfit(Surv(Ts, event) ~ 1)
  n       <- length(Ts)
  lambda0 <- numeric(n)
  
  # For each event time, a Nelson–Aalen jump = 1/n_risk
  ev_idx  <- which(event == 1)
  Ts_ev   <- sort(Ts[ev_idx]) #order
  
  for (i in seq_along(Ts_ev)) {
    n_risk         <- sum(Ts >= Ts_ev[i])
    idx_i          <- which(Ts == Ts_ev[i] & event == 1)[1]
    lambda0[idx_i] <- 1 / n_risk
  }
  
    return(lambda0)
}

lambda0= initialiser_lambda0 (Ts, event)

######### estimation


cox2 = fit_Cox_2  (beta0,lambda0,Ts,event,XB,
                   Q,tol= 1e-6,maxits = 300 )

lambda_0 = cox2$lambda_0
  
fit_true = coxph(Surv(Time,delta)~.,data = data_true)
fit_true
  
