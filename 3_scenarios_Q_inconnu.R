

############### SCENARIOS
nA_sample = c(100,500,1000) # make sure we have nA >= 100
nB_sample = 2*nA_sample
p = 3
nsim =500
beta = c(0.5,-0.5,1 )
R=c(4,6)

#(1.610541 3.039472 6.341335) for nA=100

censor_sample=c(1.553309,2.948469,6.392605)# for 40%,30%,20%

# min value of K for which we have a unique identification

K_min= 21 # for nA=1000

K_sample = cbind((K_min-9),(K_min-5),(K_min) )# make sure we have (K_min-9) >= 2 

# probability for matching variable
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


