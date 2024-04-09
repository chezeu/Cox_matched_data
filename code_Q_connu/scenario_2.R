


library(simsalapar)
library(doParallel)


sigma=1
alpha=1
C_sample=cbind(c(0.6,0.2,0.2),c(0.7,0.15,0.15),c(0.8,0.1,0.1))
m=7
n=4
p=2


beta=c(0.5,-0.5)
lambda0<- rep(0.5,n)
beta0<- c(0.1,0.1) 

nsim=5

source("scenario_Q_connu.R")
source("data_generate.R")
source("naive_method 1.R")
source(" weighted-sum_method1.R")
source(" estimate_method1.R")
source(" simulations_Q_connu.R")
###################

estimates_survival<- function(nsim,n,m,p,sigma,alpha,C,beta, lambda0, beta0){
  
  coef_true_s<-matrix(0,nrow = nsim, ncol = p)
  
  coef_naive_s<- matrix(0,nrow = nsim, ncol = p)
  
  coef_w_sum1_s<- matrix(0,nrow = nsim, ncol = p)
  
  coef_w_sum2_s<- matrix(0,nrow = nsim, ncol = p)
  
  coef_estimate_s <- matrix(0,nrow = nsim, ncol = p)
  
  converge_estimate_s<- vector()
  
  for (i in 1:nsim){
    
    surv_data =Generate_data(m,n,beta)
    
    g_naive <- doOne2d_equat_estimat(n,m,C,beta,surv_data)
    coef_naive_s[i,] <- g_naive$coef_naive2
    
    g_sum1 <- doOne2d_w_sum1 (n,m,surv_data,C)
    coef_w_sum1_s[i,] <- g_sum1$coef_w_sum1
    
    g_sum2 <- doOne2d_w_sum2 (n,m,surv_data,C)
    coef_w_sum2_s[i,] <- g_sum2$coef_w_sum2
    
    g_estim <- doOne2d_estimate (n,m,C,beta0,lambda0,surv_data)
    coef_estimate_s [i,] <- g_estim$coef_estimate
    coef_true_s[i,] <- g_estim$coef_true
    
  }
  
  
  return(list( coef_true_s1=coef_true_s[,1],coef_true_s2=coef_true_s[,2],coef_naive_s1=coef_naive_s[,1],
               coef_naive_s2=coef_naive_s[,2], coef_w_sum1_s1=coef_w_sum1_s[,1],coef_w_sum1_s2=coef_w_sum1_s[,2],
               coef_w_sum2_s1=coef_w_sum2_s[,1], coef_w_sum2_s2=coef_w_sum2_s[,2],  
               coef_estimate_s1=coef_estimate_s[,1],coef_estimate_s2=coef_estimate_s[,2]))
  
} 

######################  Simulation parameters ############################


vlis_n <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    n = list(type="grid", value = c(500,1000,2000)),
    alpha = list(type="frozen", value = 0.85),
    val_size = list(type="frozen", value = 0.1)
  )
  return(vList)
}

vlis_alpha <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    n = list(type="frozen", value = 1000),
    alpha = list(type="grid", value = c(0.75,0.85,0.95)),
    val_size = list(type="frozen", value = 0.1)
  )
  return(vList)
}

vlis_n_alpha <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    n = list(type="grid", value = c(500,1000,2000)),
    alpha = list(type="grid", value = c(0.75,0.85,0.95)),
    val_size = list(type="frozen", value = 0.1)
  )
  return(vList)
}


vlis_sampling <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    N = list(type="frozen", value = 6666),
    alpha = list(type="frozen", value = 0.85),
    val_size = list(type="frozen", value = 0.1),
    sampling_type= list(type="grid", value = c(0,1,2,3))
  )
  return(vList)
}

vlis_informative <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    n = list(type="frozen", value = 1000),
    alpha = list(type="grid", value = c(0.75,0.85,0.95)),
    val_size = list(type="frozen", value = 0.1),
    informative = list(type="grid", value = c(1,2))
  )
  return(vList)
}



vlis_sens <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    n = list(type="grid", value = c(500,1000,2000))
  )
  return(vList)
}

runSims <- function(vList=vlis1(),doOne = doOne1, seedList){
  res <- simsalapar::doForeach(vList,  doOne = doOne,  cluster=makeCluster(8, type="PSOCK"), seed = seedList)
  return(res)
}


#setwd("C://Users//thhvo.BCOM//Documents//Git//Second_article_ver2//Program")
nsim = 1000


#res_2d_mixed_laee_and_crude <- runSims(vlis_n_alpha(nsim), doOne = doOne2d_mixed_laee_and_crude, seedList = 1:nsim)
#save(res_2d_mixed_laee_and_crude, file="0512_1000s_res_2d_mixed_laee_and_crude.RData")

#res_2d_mixed_aee_sens <- runSims(vlis_sens(nsim), doOne = doOne2d_mixed_aee_sensitivity, seedList = 1:nsim)
#save(res_2d_mixed_aee_sens, file="0517_1000s_res_2d_mixed_aee_sens.RData")

#res_2d_mixed_aee_alpha <- runSims(vlis_alpha(nsim), doOne = doOne2d_mixed_aee, seedList = (1:nsim)+1000)
#save(res_2d_mixed_aee_alpha, file="0519_10000s_res_2d_mixed_aee_alpha.RData")

#res_2d_mixed_aee_n <- runSims(vlis_n(nsim), doOne = doOne2d_mixed_aee, seedList = (1:nsim)+1000)
#save(res_2d_mixed_aee_n, file="0519_10000s_res_2d_mixed_aee_n.RData")

res_2d_mixed_aee_sampling <- runSims(vlis_sampling(nsim), doOne = doOne2d_mixed_aee_sampling, seedList = 1:nsim)
save(res_2d_mixed_aee_sampling, file="0601_1000s_res_2d_mixed_aee_sampling.RData")
