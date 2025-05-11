


############### SCENARIOS

nA_sample = c(100,500,1000) # must be highest than 30
nB_sample = 2*nA_sample
p = 3
nsim =500
beta = c(0.5,-0.5,1)

#(1.610541 3.039472 6.341335) for nA=100 (no difference with nA=1000) 
#(0.2890068, 0.4817086, 0.8366007)
censor_sample = c(0.2924500,0.4020902,0.5739009)# for 40%,30%,20%

############## Table of scenarios

#scenarios = NULL

#for (i in 1:length(nA_sample)){
  #nA = nA_sample[i]
 # nB= nB_sample[i]
  

#  for (j in 1: length(nA_sample)) {
  #  scenarios = rbind(scenarios,c(nsim,nB,nA))
 # }
  
#}

#colnames(scenarios) = c("nsim","nB","nA")

#scenarios=data.frame(scenarios)


##################### New scenario #####################

alpha_p = c(0.8,0.9,1) # proportion of individual in (1,0,0)


scenarios = NULL

for (i in 1:length(nA_sample)){
  nA = nA_sample[i]
  nB= nB_sample[i]
  censor= censor_sample[2]
  
  for (j in 1: length(nA_sample)) {
    scenarios = rbind(scenarios,c(nsim,nB,nA,alpha_p[j],censor))
  }
  
}

colnames(scenarios) = c("nsim","nB","nA","alpha","censor")
addition = rbind(c(nsim,nB_sample[3], nA_sample[3] ,alpha_p[2], censor_sample[1]),
                  c(nsim,nB_sample[3] , nA_sample[3] ,alpha_p[2], censor_sample[3]),
                 c(nsim,4*nA_sample[3] , nA_sample[3] ,alpha_p[2], censor_sample[2]),
                 c(nsim,6*nA_sample[3] , nA_sample[3] ,alpha_p[2], censor_sample[2]))
scenarios=rbind(scenarios,addition)
scenarios=data.frame(scenarios)


scenarios=scenarios[c(2,10,11),]
scenarios$nA=rep(400,3)
scenarios$nB=rep(800,3)