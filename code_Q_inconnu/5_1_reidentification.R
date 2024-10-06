
nA=1000
nB=2*nA
K_min= log(nB)/log(2) # min value of K for which we have a unique identification

K=c(seq(2,K_min,1), K_min+1)

N=nB/(2^K)
plot(K,N,type = "l",col="blue",
     ylab = " Degree of identification N",
     xlab = " Number of matching variables K")
abline(h=1, col = "red")

 plot(K,1/N,lty = 1,type = "l", 
      ylab = " Probability", xlab = " Number of matching variables K")# proba d'obtenir un K-vecteur de B (size n_B)
###################################################
 
 nA=1000
 nB=2*nA
 K=10

identification<- function(k) { nB/(2^k) }
 N= sapply(1:K, identification)

 plot(1:K,N,type = "l",col="blue",
      ylab = " Degree of identification N",
      xlab = " Number of matching variables K")
 abline(h=1, col = "red")
 
 K_min= log(nB)/log(2)
 