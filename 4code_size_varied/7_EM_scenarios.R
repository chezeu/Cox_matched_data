


############### SCENARIOS
R=c(2,4,6)
nA = 1000 # make sure we have nA >= 100
nB_sample = c(R[1]*nA ,R[2]*nA ,R[3]*nA )

p = 3
nsim =5
beta = c(0.5,-0.5,1)

K_min= c(22,24,25)
#K_min= round(log(nB_sample)/log(2)) # min value of K for which we have a unique identification

#K_sample=  t(matrix(rep(c(9,18,27), p), ncol=3))# MULTIPLE OF 3, vary K

K_sample = cbind((K_min-8),(K_min-4),(K_min) )# make sure we have (K_min-5) >= 2 

#prevalence_sample = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)

 prevalence_sample = 0.5 # probability of bernoulli

############## Table of scenarios

scenarios = NULL
for (k in 1:length(nB_sample)) {
nB =  nB_sample[k]

  for (i in 1:nrow(K_sample)){
    K= K_sample[k,i]
    
    scenarios = rbind(scenarios,c(nsim,nB,nA,K))

  }
  
}
addition= rbind(c(nsim,nB_sample[2], nA, K_sample[1,2]),
                c(nsim,nB_sample[3], nA, K_sample[1,2]))
scenarios= rbind(scenarios,addition)
colnames(scenarios) = c("nsim","nB","nA","K")

scenarios=data.frame(scenarios)


##################### matrix
