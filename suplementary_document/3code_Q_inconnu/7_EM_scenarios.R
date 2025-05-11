

############### SCENARIOS
nA_sample = c(100,500,1000) # make sure we have nA >= 100
nB_sample = 2*nA_sample
p = 3
nsim =500
beta = c(0.5,-0.5,1)
R=c(4,6)
#(1.610541 3.039472 6.341335) for nA=100 (no difference with nA=1000) 
censor_sample = c(0.2924500,0.4020902,0.5739009)# for 40%,30%,20%

K_min= 22
#K_min= round(log(nB_sample)/log(2)) # min value of K for which we have a unique identification

#K_sample=  t(matrix(rep(c(9,18,27), p), ncol=3))# MULTIPLE OF 3, vary K

K_sample = cbind((K_min-10),(K_min-6),(K_min) )# make sure we have (K_min-5) >= 2 

#prevalence_sample = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)

 prevalence_sample = 0.5

############## Table of scenarios

scenarios = NULL
for (k in 1:length(K_sample)) {
KJ =  K_sample[k]

  for (i in 1:length(nA_sample)){
    nA = nA_sample[i]
    nB= nB_sample[i]
    censor= censor_sample[2]
    scenarios = rbind(scenarios,c(nsim,KJ,nA,nB,censor))

  }
}
addition = rbind( c(nsim, K_sample[2] ,nA_sample[3], nB_sample[3],censor_sample[1]),
                 c(nsim,K_sample[2] ,nA_sample[3], nB_sample[3],censor_sample[3]),
                 c(nsim,K_sample[2] ,nA_sample[3], R[1]*nA_sample[3],censor_sample[2]),
                 c(nsim,K_sample[2] ,nA_sample[3], R[2]*nA_sample[3],censor_sample[2]))
                 
scenarios=rbind(scenarios,addition)

colnames(scenarios) = c("nsim","K","nA","nB","censor")

scenarios=data.frame(scenarios)

##################### matrix
