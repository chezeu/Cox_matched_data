


############### SCENARIOS

C_sample = rbind(c(0.6,0.2,0.2),c(0.7,0.2,0.1),c(1,0,0))
nA_sample = c(100,500,1500)
n0_sample = round((20*nA_sample)/100)
nB_sample = 2*nA_sample
p = 3
nsim =1000
#gamma_sample=c(0.6,0.7, 0.8)
sigma=1
alpha=1 # alpha != 1 
beta = c(0.5,-0.5)

#lambda0 = rep(0.5,n)
#beta0 = c(0.1,0.1) 
#beta = c(0.5,-0.5) 

############## Table of scenarios

scenarios = NULL

for (i in 1:length(nA_sample)){
  nA = nA_sample[i]
  n0 = n0_sample[i]
  nB= nB_sample[i]
  
  for (j in 1: nrow(C_sample)) {
    #for (j in 1: length(gamma_sample)) {
    scenarios = rbind(scenarios,c(nsim,nB,nA,C_sample[j,]))
  }
  
}

colnames(scenarios) = c("nsim","nB","nA","Prob_1","Prob_2","Prob_3")

scenarios=data.frame(scenarios)


##################### matrix

