
source("scenario_Q_connu.R")

scenarios_boxplot=scenario_Q_connu

pname = "tau"
library(ggplot2)
library(gridExtra)

levelsvec=tauvec

for (i in (1:nrow(scenarios))){
  nsim_sample=scenarios[i,1]
  n_sample=scenarios[i,2]
  m_sample=scenarios[i,3]
  p_sample=scenarios[i,4]
  sigma_sample=scenarios[i,5]
  alpha_sample=scenarios[i,6]
  C=vector()
  C[1]=scenarios[i,7]
  C[2]=scenarios[i,8]
  C[3]=scenarios[i,9]
  
  
  results_sample<-estimates_survival (nsim=nsim,n=n_sample,m=m_sample,p=p_sample,
                                      sigma=sigma_sample,alpha=alpha_sample,C=C,beta, lambda0, beta0)
  
  
  flname=paste0("./Results_sample/",flnameEst,"_",
                "nsim=",nsim,"_",
                "n=",n_sample,"_",
                "m=",m_sample,"_",
                "p=",p_sample,"_",
                "sigma=",sigma_sample,"_",
                "alpha=",alpha_sample,"_",
                "C=",C,"_",
                ".RData")
  
  load(flname)
  
  boxplotHarrel=cbind(rep(1,times=nsim),
                      rep(scenarios_boxplot[[pname]][i],times=nsim),
                      results_sample$coef_true_s1)
  
  boxplotHarrel=cbind(rep(2,times=nrep),
                      rep(scenarios_boxplot[[pname]][i],times=nsim),
                      results_sample$coef_true_s2)
  
  
} 


titlename=paste0("n=",unique(scenarios_boxplot[,2]),
                 ", ",pname)

colnames(boxplotHarrel)=c("Estimator",pname,"Harrel")
boxplotHarrel=data.frame(boxplotHarrel)

boxplotHarrel[,1]=factor(boxplotHarrel[,1],
                         levels=1:2,
                         labels=c("Cox","RF"))
boxplotHarrel[,2]=factor(boxplotHarrel[,2],
                         levels=levelsvec,
                         labels = as.character(levelsvec*100))
# Boxplots for beta
plotHarrel=ggplot(boxplotHarrel,aes(x=get(pname),y=Harrel,fill=Estimator))+
  geom_boxplot()+
  xlab("Fixed % of censoring")+
  ylab("Harrel")+
  geom_hline(yintercept=250,lty=1,col="orange")+
  ggtitle(titlename)

