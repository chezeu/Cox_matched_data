############### SCENARIOS

nA_sample = c(100,200,300)
nB_sample = 2*nA_sample
nA_block_sample=rbind(c(10,40,50), c(20,80,100), c(50,100,150))
nB_block_sample = 2*nA_block_sample
p = 3
nsim =5
beta = c(0.5,-0.5,1)
K_sample= c(7,14,21)# MULTIPLE OF 3
#prevalence_sample = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)

prevalence_sample = c(0.5, 0.5,0.5,0.5,0.5,0.5,0.5)
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
    nA_block=nA_block_sample[i,]
    nB_block=nB_block_sample[i,]
    scenarios = rbind(scenarios,c(nsim,K,nA,nB,nA_block,nB_block ))
    
  }
  
}
colnames(scenarios) = c("nsim","K","nA","nB",paste("nA_block",1:length(nA_sample), sep = ""),
                        paste("nB_block",1:length(nA_sample), sep = ""))

scenarios=data.frame(scenarios)

########### matrix

