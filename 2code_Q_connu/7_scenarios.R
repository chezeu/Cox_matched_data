


############### SCENARIOS

C_sample = rbind(c(0.4,0.3,0.3),c(0.8,0.1,0.1),c(1,0,0))
nA_sample = c(30,50,100) # n_A must be bigger than 30
nB_sample = 2*nA_sample
p = 3
nsim =5
#gamma_sample=c(0.6,0.7, 0.8)
beta = c(0.5,-0.5,1)


############## Table of scenarios

scenarios = NULL

for (i in 1:length(nA_sample)){
  nA = nA_sample[i]
  nB= nB_sample[i]
  
  for (j in 1: nrow(C_sample)) {

    scenarios = rbind(scenarios,c(nsim,nB,nA,C_sample[j,]))
  }
  
}

colnames(scenarios) = c("nsim","nB","nA","Prob_1","Prob_2","Prob_3")

scenarios=data.frame(scenarios)


##################### matrix

