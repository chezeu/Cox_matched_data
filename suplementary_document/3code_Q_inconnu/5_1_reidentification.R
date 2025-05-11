
nA=1000
nB=2*nA
K_min= round(log(nB)/log(2)) # min value of K for which we have a unique identification

K=1:K_min
N=round(nB/(2^K)) #Pairs with common variables K
base= cbind(K,N)
base= base[K_min:(K_min-8), ]
plot(base[,1] ,base[,2],type = "l",col="blue",
     ylab = " Idenfication degree N",
     xlab = " Number of matching variables K")
abline(h=1, col = "red")

plot(K,1/N,lty = 1,type = "l", 
     ylab = " Probability", xlab = " Number of matching variables K")# proba d'obtenir un K-vecteur de B (size n_B)




############
############## NEW version ##########################
R=2
#R=rep(3,3)
nA_sample = c(100,500,1000) # make sure we have nA >= 100
nB_sample = R*nA_sample
# nB=2000

f= function(K,nB){
  N= round((nB^2)*(1/(2^K))*(1-(1/(2^K)))^(nB-1))-1
  return(N)
}
#f(7.65)
# f(22)

#solution dans l'intervalle [11, 50] où f est strictement 
# décroissante, f(11) et f(50) sont de signes crontraires
# methode de dichotomie sur l'intervalle [11, 50]
a=11
b=50

solution <- function(a,b,nB){
  error=5
  while (error >= 10e-5) { 
    c= round((a+b)/2)
    if (f(a,nB) * f(c,nB) > 0) { a = c}
    if (f(a,nB) * f(c,nB) < 0) {b=c}
    if (f(c,nB)== 0){a=b=c}
    error = abs(b-a)
  }
  return (a)
}

S=NULL
for (j in 1:3) { 
  nB =nB_sample[j]
  s= solution(a,b,nB)
  S= c(S,s)
}
S  # 16 20 22
f(S[3],nB_sample[3])
#f(7.65)=0 mais on ne peut pas d'arrondir cette valeur

################ plot for the first case
S=22
K=15:S
nB=2000
N= f(K,nB)+1 #unicité
base= cbind(K,N)
plot(base[,1] ,base[,2],type = "l",col="blue",
     ylab = " Degree of idenfication N",
     xlab = " Number of matching variables K")
abline(h=1, col = "red")
