


############### SCENARIOS

nB_sample = c(200,300,500)
nA_sample = round((80*nB_sample)/100)
p = 2
nsim =30
beta = c(0.5,-0.5)
K= 20# MULTIPLE OF 3
error=0.04
prevalence = rep(c(0.05,0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.7), K/10)
#lambda0 = rep(0.5,n)
#beta0 = c(0.1,0.1) 
#beta = c(0.5,-0.5) 

############## Table of scenarios

scenarios = NULL

for (i in 1:length(nA_sample)){
  nA = nA_sample[i]
  nB= nB_sample[i]
 
    scenarios = rbind(scenarios,c(nsim,nB,nA))
  

}

colnames(scenarios) = c("nsim","nB","nA")

scenarios=data.frame(scenarios)


##################### matrix

