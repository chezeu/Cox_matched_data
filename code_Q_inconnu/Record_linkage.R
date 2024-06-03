library(RecordLinkage)
library(plyr)
library(dplyr)

library(reclin2)

linkage_data <- function(databaseB, databaseA) {

  
  pp <- pair_minsim(databaseA, databaseB, on= c(paste("R", 1:K, sep = "")), minsim = 0)
  pp <- compare_pairs(pp, on=  c(paste("R", 1:K, sep = "")), inplace = TRUE)
  
  ## classification
  model <- problink_em ( ~R1+R2+R3+R4+R5+R6,
                         data=pp)
  pp <- predict(model, pp, type ="mpost", add = TRUE, binary = TRUE)
  
  matrix_norm <- pp$mpost
  nA<- nrow(databaseA)
  nB <-  nrow(databaseB)
  matrix_norm <- t( matrix(matrix_norm,ncol=nA,nrow=nB))
  #chaque ligne correspond aux probabilites pour chaque elt de la base B
  naive_id = apply(matrix_norm, 1, FUN= function(x){which.max(x)})
  idA=databaseA$id
  vect_sum <- apply(matrix_norm,1, sum)
  
  Q<-sapply(1:nA, FUN=function(i){
    matrix_norm[i,]/vect_sum[i]
  })  
  Q= t(Q)
  
  
  ##creat a linked data
  #linked_data <- link(pp)
  return(list( Q=Q,naive_id =naive_id ,pp=pp))
}
