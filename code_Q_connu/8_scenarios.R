
setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")

source("7_monte_carlo.R")
############### SCENARIOS

C_sample = cbind(c(0.6,0.2,0.2),c(0.7,0.15,0.15),c(0.8,0.1,0.1))
m = c(500,1000,2000)
n = round((80*m)/100)
p = 2
nsim = 100 
beta = c(0.5,-0.5)

#lambda0 = rep(0.5,n)
#beta0 = c(0.1,0.1) 

############## Tableau

scenarios = NULL
for (i in 1:ncol(C_sample)){
  scenarios = rbind(scenarios,c(nsim,n[2],m[2],C_sample[,i]))
}
for (i in 4:6){
  scenarios = rbind(scenarios,c(nsim,n[i-3],m[i-3],C_sample[,3]))
}

colnames(scenarios) = c("nsim","n","m","Prob_1","Prob_2","Prob_3")

scenarios=data.frame(scenarios)


#####################

setwd("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu")

for (i in (1:nrow(scenarios))){
  
  nsim=scenarios[i,1]
  n=scenarios[i,2]
  m=scenarios[i,3]
  C=vector()
  C[1]=scenarios[i,4]
  C[2]=scenarios[i,5]
  C[3]=scenarios[i,6]
  
  lambda0 = rep(0.5,n)
  beta0 = c(0.1,0.1) 
  
  results_sample = estimates_survival(nsim,n,m,C,p,beta,Q,lambda0,beta0) #monte carlos
  
  filename = paste0("C:/Users/fchezeut/Documents/GitHub/Cox_matched_data/code_Q_connu/Results/","scenario=",i,
                    "_m=",m, "_prob1=", C[1], ".Rdata")
  
  save(results_sample,file = filename)
  
}
