


############### SCENARIOS

nA_sample = c(100,500,1000)
nB_sample = 2*nA_sample
p = 3
nsim =500
beta = c(0.5,-0.5,1)
K_sample= c(9,18,27)# MULTIPLE OF 3
prevalence_sample = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)

#prevalence = rep(0.5, 0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
#lambda0 = rep(0.5,n)
#beta0 =rep(0.1,p) 
#beta = c(0.5,-0.5,1) 

############## Table of scenarios

scenarios = NULL
for (k in 1:length(K_sample)) {
  K = K_sample[k]
  for (i in 1:length(nA_sample)){
    nA = nA_sample[i]
    nB= nB_sample[i]
    scenarios = rbind(scenarios,c(nsim,K,nA,nB))
    
  }
  
}
colnames(scenarios) = c("nsim","K","nA","nB")

scenarios=data.frame(scenarios)


##################### matrix

