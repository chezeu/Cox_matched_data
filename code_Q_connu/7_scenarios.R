


############### SCENARIOS

C_sample = rbind(c(0.6,0.2,0.2),c(0.7,0.2,0.1),c(0.8,0.1,0.1))
m_sample = c(100,500,800)
n_sample = round((80*m_sample)/100)
p = 2
nsim = 100 
beta = c(0.5,-0.5)

#lambda0 = rep(0.5,n)
#beta0 = c(0.1,0.1) 

############## Table of scenarios

scenarios = NULL

for (i in 1:length(m_sample)){
  m= m_sample[i]
  n= n_sample[i]
  for (j in 1: nrow(C_sample)) {
    scenarios = rbind(scenarios,c(nsim,m,n,C_sample[j,]))
  }

}

colnames(scenarios) = c("nsim","m","n","Prob_1","Prob_2","Prob_3")

scenarios=data.frame(scenarios)


#####################
