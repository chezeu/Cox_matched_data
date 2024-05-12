


############### SCENARIOS

C_sample = rbind(c(0.6,0.2,0.2),c(0.7,0.2,0.1),c(1,0,0))
n_sample = c(100,500,900)
n0_sample = round((20*n_sample)/100)
m_sample = n0_sample + n_sample
p = 2
nsim =100
#gamma_sample=c(0.6,0.7, 0.8)
sigma=1
alpha=1 # alpha != 1 
beta = c(0.5,-0.5)

#lambda0 = rep(0.5,n)
#beta0 = c(0.1,0.1) 
#beta = c(0.5,-0.5) 

############## Table of scenarios

scenarios = NULL

for (i in 1:length(n_sample)){
  nA = n_sample[i]
  n0 = n0_sample[i]
  nB= m_sample[i]
 
  for (j in 1: nrow(C_sample)) {
  #for (j in 1: length(gamma_sample)) {
    scenarios = rbind(scenarios,c(nsim,nB,nA,C_sample[j,]))
  }

}

colnames(scenarios) = c("nsim","nB","nA","Prob_1","Prob_2","Prob_3")
#colnames(scenarios) = c("nsim","nB","nA","gamma")

scenarios=data.frame(scenarios)


##################### matrix

