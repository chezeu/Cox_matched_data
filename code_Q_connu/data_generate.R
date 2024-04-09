
## generer les bases de donnEes######################################################
Generate_data <- function(m,n,beta){
  # set.seed(42)
  ##covarites data
  X1 <- rnorm(m,0,1)
  X2 <- rbinom(m,size = 1, prob = 0.7)
  X <- as.matrix(cbind(X1, X2))
  
  U <- runif(m,0,1) # Time to event variable
  Tt = sigma*(-log(U)/exp((X%*%beta)))^(1/alpha)  #lambda=1
  # Constant right censored
  #c = quantile(Tt[,1], probs = rho)
  c=1.777907
  
  Time <- pmin(Tt,c)
  delta <- as.numeric(Tt<=c)
  surv_data <- data.frame(Time,delta,X)
  colnames(surv_data)=c('Time','delta',paste("X", 1:p, sep = ""))
  
  # true_data <- surv_data[1:n,]
  
  # database_A <- true_data[,1:2]  #base de donnee de A
  # XB <- surv_data[,3:dim(surv_data)[2]] #base de donnee de covariables
  return(surv_data=surv_data )
}
