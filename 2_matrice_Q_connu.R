
# matrix with variation of the number of possible match
genere_matrix = function(nA,nB,alpha){

Q1 = matrix(0, nA, nB)

t1= round((alpha*nA)) #individuals with (1,0,0,0)
t2= nA - t1 #individuals with different possible matche
v1= c(0.6,0.2,0.2) # 3 possible matches
v2= c(0.3,0.2,0.2,rep(0.1,3)) # 6 possible matches
#v3= c(0.2,rep(0.1,8) ) # 9 possible matches
v3= rep(1/9, 9) # 9 possible matches
nbre= c(3,6,9)

if(v_matrix==nbre[1]){ 
  for (i in 1:t1) {
    vec = sample(1:nB, 1)
    Q1[i,vec] = 1
  }
  if(t1!=nA){ 
    for (i in (t1+1):nA ){
      vec = sample(1:nB, nbre[1])
      Q1[i,vec] = v1
    }
  }
}

if(v_matrix== nbre[2]){ 
  for (i in 1:t1) {
    vec = sample(1:nB, 1)
    Q1[i,vec] = 1
  }
  if(t1!=nA){ 
    for (i in (t1+1):nA ){
      vec = sample(1:nB, nbre[2])
      Q1[i,vec] = v2
    }
  }
}

if(v_matrix==nbre[3]){ 
  for (i in 1:t1) {
    vec = sample(1:nB, 1)
    Q1[i,vec] = 1
  }
  if(t1!=nA){ 
    for (i in (t1+1):nA ){
      vec = sample(1:nB, nbre[3])
      Q1[i,vec] = v3
    }
  }
}
return(Q1)
}
