


############### SCENARIOS

nA_sample = c(100,500,1000) # must be highest than 30
nB_sample = 2*nA_sample
p = 3
nsim =5
beta = c(0.5,-0.5,1)

#(1.610541 3.039472 6.341335) for nA=100 (no difference with nA=1000) 

censor_sample=c(1.553309,2.948469,6.392605)# for 40%,30%,20%
v_matrix=c(3,6,9) # number of possible matches

##################### New scenario #####################

alpha_p = c(0.8,0.9,1) # proportion of individual in (1,0,0)


scenarios = NULL

for (i in 1:length(nA_sample)){
  nA = nA_sample[i]
  nB= nB_sample[i]
  censor= censor_sample[2]
  matrix= v_matrix[1]
  
  for (j in 1: length(nA_sample)) {
    scenarios = rbind(scenarios,c(nsim,nB,nA,alpha_p[j],censor,matrix))
  }
  
}

colnames(scenarios) = c("nsim","nB","nA","alpha","censor","v_matrix")
addition = rbind(c(nsim,nB_sample[3], nA_sample[3] ,alpha_p[2], censor_sample[1],v_matrix[1]),
                  c(nsim,nB_sample[3], nA_sample[3] ,alpha_p[2], censor_sample[3],v_matrix[1]),
                 c(nsim,nB_sample[3], nA_sample[3] ,alpha_p[2], censor_sample[2],v_matrix[2]),
                 c(nsim,nB_sample[3], nA_sample[3] ,alpha_p[2], censor_sample[2],v_matrix[3]))
scenarios=rbind(scenarios,addition)
scenarios=data.frame(scenarios)

