
## compute the matrix sigma = (I-J)^-1

compute_sigma<- function(beta0, lambda0, Ts, event, XB,  prob ) {
  
  XB  <- as.matrix(XB)
  nA  <- length(Ts)
  nB   <- nrow(XB)
  p   <- ncol(XB)
  
  eXbeta = as.vector(exp(XB %*% beta0)) # nB
  
    #   pi_ij * exp(beta^top x_j)  (nA x nB)
  W <- prob * matrix(eXbeta, nrow = nA, ncol = nB, byrow = TRUE)
  
  I_hat =matrix(0,p,p)
  for (i in 1:nA) {
    
    Lambda_i = sum(lambda0[Ts <= Ts[i]])
    
    for (j in 1:nB) {
      
      mat_X = XB[j,]%*% t(XB[j,]) # pxp
      I_hat = I_hat +  (W[i, j] * Lambda_i * mat_X)
    }
     } 
  I_hat = I_hat/nA
  
  J_hat <- matrix(0, p, p)
  
  for (i in seq_len(nA)) {
    
    if (event[i] == 0) next   # null value if no event occurs
    
    lambda0_i =  lambda0[i] 
    
    # indicator 1_{T_i <= T_k} for k  (n_A x 1)
    at_risk <-  Ts[i] <= Ts   # (n_A x 1)
    
    # sum_k sum_l pi_kl * x_l * exp(beta^top x_l) * 1_{T_i <= T_k}
    
    w_at_risk <- colSums(W[at_risk, , drop = FALSE])  # (nB)
    mu_i      <- as.vector(w_at_risk %*% XB) / nA     # (p)
    
     # contribution : delta_i * (mu_i %o% mu_i)  (p x p)
    J_hat <- J_hat + (event[i] * outer(mu_i, mu_i)* (lambda0_i)^2)
  }
  J_hat =nA*J_hat
  
  #### sigma
  sigma1= I_hat - J_hat 
   
  # --- Check before reversing---
  eig <- eigen(sigma1, symmetric = TRUE)$values
  cat("eigen value min/max de (HESS) :", 
      round(min(eig), 6), "/", round(max(eig), 4), "\n")
  
  if (min(eig) < 1e-8) {
    warning("Matrix quasi-singular : Check variables")
    # Regularisation  
    sigma1 <- sigma1 + 1e-6 * diag(p)
  }
  
  sigma = solve(sigma1) #invers
  
  return( sigma = sigma)
  
}
    ########### tABLE OF RESULTS

EM_results_table <- function(beta0,lambda0,Ts,event,XB,Q) {
  
  p       <- ncol(as.matrix(XB))
  nA       <- length(Ts)
  
  lambda2 <- Funct_lambda2(lambda0, Ts)
  prob    <- Functio_prob(beta0, lambda2, Ts, event, XB, Q)
  
  var <-   compute_sigma (beta0, lambda0, Ts, event, XB,  prob )
  var  = var /nA
  
  var_beta = diag(var)
  
  se=  sqrt(var_beta)
  # z-score and p-valeur
  
  z_score =(beta0 / se)
  
  # H_0: beta_0=0 ,the variable has no significant effect
  p_val  <- 2 * (1-pnorm(abs(z_score))) 
  
  tab <- data.frame(
    coef = round(beta0, 4),
    "exp(coef)" = round(exp(beta0), 4),
    "se(coef)"    = round(se, 4),
    z  = round(z_score, 3),
    "Pr(>|z|)"=  p_val ,
    check.names = FALSE
  )
  
  rownames(tab) <- colnames(as.matrix(XB))
  print(tab)
  return(tab)
}

############  estimation parameters with z score
fit_Cox_2 <- function(beta0,lambda0,Ts,event,XB,
                      Q,tol= 1e-6,maxits = 300 ){ 
  
  fit_estimate = Func_itteration(beta0,lambda0,Ts,event,XB,
                                 Q,tol= 1e-6,maxits = 300)
  beta_0   <- fit_estimate$beta0
  lambda_0 <- fit_estimate$lambda0
  
  tab_final <- EM_results_table(
    beta_0,lambda_0,Ts,event,XB,Q)
    
  cat(
    "The number of iterations is:", fit_estimate$it, "\n",
    "n =", length(Ts), "\n",
    "Number of events =", sum(event), "\n"
  )
 
  return(list( tab_final=tab_final, lambda_0=lambda_0))
}
 