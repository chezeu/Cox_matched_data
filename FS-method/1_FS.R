

########### record linkage#####################################


linkage_function<- function(datasetA,datasetB,var_block,matching_variables){
  
  # pairs with blocking variable
  p <- pair_blocking (datasetA, datasetB,var_block , FALSE)
  #p <- print(p)
  
  #  comparison matrix of pairs
  p <- compare_pairs ( p, on= matching_variables, 
                       comparators = list(fname_c2 = jaro_winkler(),
                  lname_c1 = jaro_winkler(), lname_c2  = jaro_winkler() ) )
  #p <- print(p)
  
  ## classification
  p0 = 1/nrow(datasetB)
  
  model <- problink_em ( ~ fname_c1+fname_c2 +lname_c1 +lname_c2 +by +bm+bd , data=p,
                         mprobs0 = list(0.9),uprobs0 = list(0.02),p0 = p0 ,tol = 1e-05,
                         mprob_max = 0.999, uprob_min = 1e-04)
  
  #print(model)       
  
  p_compar <- predict(model, p, type ="mpost", add = TRUE, binary = TRUE)
  prob_match <- p_compar[,c(1:2,10)]
  
  return(prob_match=prob_match)
}

#prob_match = linkage_function (datasetA,datasetB,var_block)

####################### Generate data with blocks#################

function_block <- function( prob_match, datasetA, datasetB){
  
     # generate the blocks A
p1 = unique( prob_match$.x)
new_dataA = datasetA[p1,]

blockA = unique(new_dataA$fname_c1) # number of blocks A with variable fname_c1
new_dataA$idA = 1:nrow(new_dataA)

     # generate the blocks B
VEC_idB = vector()
for (i in 1: length (blockA)) {
  ni = which(datasetB$fname_c1 == blockA[i])
  VEC_idB = c(VEC_idB,ni)
}
datB2= datasetB[VEC_idB,] 
n=as.numeric(VEC_idB)
rest_B = datasetB[-n,]
new_dataB = rbind(datB2, rest_B)
new_dataB$idB = 1:nrow(new_dataB) # new data B

return( list(blockA =blockA ,new_dataA=new_dataA, new_dataB=new_dataB))
}

 #data_block = function_block ( prob_match, datasetA, datasetB)
 #blockA =data_block$blockA
 #new_dataA= data_block$new_dataA
 #new_dataB = data_block$new_dataB
 
# ######BLOCK MATRIX################################
#### Q is with the new databases
 
Matrix_block <- function(blockA,datasetA,datasetB,new_dataA,new_dataB,prob_match){
 
  n_A= nrow(datasetA)
  n_B= nrow(datasetB) 
  Q = matrix(0,nrow = n_A,ncol = n_B) # general matrix
  p=1
  q=1
for (i in 1:length(blockA) ) {
  
  name_A = blockA[i]
  n =  which(datasetA$fname_c1==name_A)  
  m = which(datasetB$fname_c1==name_A)  
  idA = new_dataA[p:(p+length(n)-1),]$idA
  idB = new_dataB[q:(q+length(m)-1),]$idB
  
  if (length(n)==1){
  mprob = prob_match$mpost[ which(prob_match$.x == n) ]
  Q_I = matrix(mprob,length(n),length(m)) # block matrix
  }
  
  if (length(n)!=1){
    Q_I = matrix(0,nrow = length(n), ncol = length(m) )
    for (j in 1: length(n)) {
      mprob = prob_match$mpost[ which(prob_match$.x == n[j]) ]
      Q_I[j,] = mprob
    }
  }
  
  Q[idA, idB] = Q_I
  p=p+length(n)
  q=q+length(m)
}

  #normalize Q 
  number = apply(Q,1, sum) # number of common individual
  Q= Q/number
  
return( Q=Q)
}

#Q= Matrix_block (blockA,datasetA,datasetB,new_dataA,new_dataB,prob_match)
################ 

#naive data
#naive_id = apply( Q,1, which.max)
#data_naive=NULL
#for (i in 1:nB) {
#  a = naive_id[i]
#  data = new_dataB[which(new_dataB$idB == a),]
#  data_naive = rbind(data_naive, data)
#}


