
##################################### 1 block ###########################
generate_1block <- function(K,nA, nB, p,prevalence, beta,block){
  
  blockA = rep(block, nA)
  blockB = rep(block, nB)
  
  # Survival data (Database B)
  ft<- Generate_data(K,nA,nB,prevalence, min_prev = 0.01)
  surv_data<-ft$surv_data
  XB <- ft$XB #Only covariates
  
  ###### 
  datA<-ft$datA  # matching variables
  datB<-ft$datB
  datA<- cbind(datA, block= blockA)
  datB<- cbind(datB, block= blockB)
  XB <- as.data.frame(cbind(XB, block = blockB))
  
   return(list( surv_data = surv_data, datA = datA, XB = XB, datB = datB))
}

# nA_block== vector size of blocks in the database A
# q = number of block
multi_block <- function( K, nA_block, nB_block, p, prevalence, beta){
  
  q = length(nA_block)
  # for one block (q=1)
  nA = nA_block[1]
  nB = nB_block[1]
 
  
  data_block = generate_1block (K,nA, nB, p,prevalence, beta,block=1)
  datA = data_block$datA
  datB = data_block$datB
  XB <- data_block$XB
  surv_data = data_block$surv_data
  
   # q >= 2  
  for (i in 2:q){
    nA = nA_block[i]
    nB = nB_block[i]
    
    data_block = generate_1block (K,nA, nB, p,prevalence, beta,block=i)
    
    datA = rbind(datA,data_block$datA)
    datB = rbind(datB,data_block$datB)
    XB = rbind(XB,data_block$XB)
    surv_data = rbind(surv_data, data_block$surv_data)
  } 
  
  return(list(datA = datA,  XB = XB, datB = datB,surv_data = surv_data))
  
}


###############  Record linkage and matrix gamma###################################

compare_binary <- function(datA, datB, K, nA_block){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid( XB.k,XA.k)
    
    gamma.k = as.numeric(temp[,1]==temp[,2])
    
    return(gamma.k)
  }
  
  #S=NULL
  comp_mat =NULL
  
  for (i in 1:length(nA_block)) {
 
  comp_mat_block = sapply(1:K, FUN = compare1, datA = datA[which(datA$block==i),], 
                          datB = datB[which(datB$block==i),], K =K)
  vect_gamma= vector()
  for (j in 1: nrow(comp_mat_block)) {
    if(sum(comp_mat_block [j,])==K){vect_gamma[j]=1 
    }else {vect_gamma[j]=0 }
    
  }
  
  #S = c(S,which(vect_gamma==1) )
  comp_mat_block2= cbind(comp_mat_block,vect_gamma,rep(i,nrow(comp_mat_block)))
  
  comp_mat=rbind(comp_mat,comp_mat_block2)
  }
  
 colnames(comp_mat)=c(paste("C_R", 1:K, sep=""), "vect_gamma", "block" )
  
  return(comp_mat)
}
##################

############### matrix Q with blocks

Matrix_Q = function(comp_mat,K,datA,datB,nA_block,nB_block){
  
  Q= matrix(0, nrow=nrow(datA),ncol = nrow(datB))
  comp_mat= data.frame(comp_mat) 
  vect_gamma= comp_mat$vect_gamma # comparison pairs
  
  l=1
  m=1
  for (q in 1:length(nA_block)) {
    nA=nA_block[q]
    nB=nB_block[q]
    mat_block = t( matrix(comp_mat$vect_gamma[which(comp_mat$block==q)],
                          ncol=nA,nrow=nB))
    
    Q[m:(nA+m-1),l:(nB+l-1)] = mat_block
    m= nA + m
    l = nB + l
    
  }
  
  #identification = apply(Q,1, FUN= function(x){which(x==1)}) # common variables
  number = apply(Q,1, FUN= function(x){length(which(x==1))}) # number of common individual
  
  Q= Q/number
  
  naive_id = apply(Q,1, which.max)
  
  return(list(Q=Q, naive_id=naive_id))
  
}
