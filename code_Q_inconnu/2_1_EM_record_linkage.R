

library(matrixStats)

###################################"" Record Linkage#######################

compare_binary <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid( XB.k,XA.k)
    
    gamma.k = as.numeric(temp[,1]==temp[,2])
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  
  
  return(comp_mat)
}


#T_id <- expand.grid(datB[,K+1],datA[,K+1] )


EM_binary <- function( comp_mat,datA, datB, K, e =0.02, tol = 1e-5, maxits = 2000){
  # Starting point
  nA =  nrow(datA)
  N <- nrow(comp_mat)
  p = nA/N
  
  prev = (colMeans(datB[,1:K]))
  
  u = ((1-e)*(1-prev) + e*prev)*(1-prev) + prev*((1-e)*prev+e*(1-prev))
  m = rep(1 - e,K)
  #u=(colSums(datA[-(K+1)])/ nA )*prev  + ((nA-colSums(datA[-(K+1)]))/ nA )*(1-prev)
  #m= ( colSums(datA[-(K+1)])/colSums(datB[-(K+1)]))*prev + ((nA -colSums(datA[-(K+1)]) )/(nB -colSums(datB[-(K+1)]) ))*(1-prev)
  #u= rep(0.5,K)
  #m=rep(0.5,K)

  
  # initializations
  comp_mat  <- comp_mat[,1:K]
  
  g = rep(0,N) # probability of being in Match  for each pair l
  it = 0
  converge = FALSE
  
  
  while ((!converge) & (it < maxits)){ 
    
    p.old = p
    m.old = m
    u.old = u
    ### E
    # Compute expectation
    
    m_mat = matrix(rep(m,N), nrow =N, byrow = TRUE)
    u_mat = matrix(rep(u,N), nrow =N, byrow = TRUE)
    
    
    probM = rowProds(m_mat^(comp_mat)*(1-m_mat)^(1-comp_mat))
    probU = rowProds(u_mat^(comp_mat)*(1-u_mat)^(1-comp_mat))
    
    
    g = p*probM/(p*probM+(1-p)*probU)
    
    ### Maximization
    g_mat = matrix(rep(g,K),ncol = K)
    
    p = sum(g)/N
    m = colSums(g_mat*comp_mat)/sum(g)
    u = colSums((1-g_mat)*comp_mat)/sum(1-g)
    
    if (length(which(m > 0.99999)) > 0) {
      m[which(m > 0.99999)] <- 0.99999
    }
    if (length(which(m < 1e-05)) > 0) {
      m[which(m < 1e-05)] <- 1e-05
    }
    if (length(which(u > 0.99999)) > 0) {
      u[which(u > 0.99999)] <- 0.99999
    }
    if (length(which(u < 1e-05)) > 0) {
      u[which(u < 1e-05)] <- 1e-05
    }
    
    
    it = it + 1
    
    converge <- (abs(p.old - p)/p.old < tol) && 
      all(abs(m.old - m)/m.old < tol) && 
      all(abs(u.old - u)/u.old < tol)
    
    if (it == maxits) {
      cat("WARNING! NOT CONVERGENT FOR LINKAGE!", "\n")
      converge = FALSE
    }
  }
  
  comp_mat.1 <- as.data.frame(comp_mat)
  comp_mat.1 $g <- g
  
  
  nA <- nrow(datA)
  nB <-  nrow(datB)
  matrice_norm <- t( matrix(g,ncol=nA,nrow=nB))
  
  naive_id = apply( matrice_norm,1, which.max)
  #chaque ligne correspond aux probabilites pour chaque elt de la base B
  
  # h <- which(  matrice_norm<=10e-9) 
  # matrice_norm[h] <- 0
  
  vect_sum <- apply(matrice_norm,1, sum)
  vect_sum <-matrix( rep(vect_sum,nB), ncol = nB,nrow = nA)
  
  ##matrice de prob normalisee
  
  Q <- matrice_norm/vect_sum
  apply(Q,1,sum)
  
  return(list(comp_mat.1=comp_mat.1,Q=Q, it = it, naive_id= naive_id))
}

