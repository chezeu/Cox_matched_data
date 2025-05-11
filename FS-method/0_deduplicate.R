
library(RecordLinkage)
library(plyr)
#library(dplyr)
#library(reclin2)

data("RLdata500")
data("RLdata10000")

s1 <- 1:100
s2 <- 1:100
datasetB <- rbind( RLdata10000[s2,],RLdata500[s1,])
datasetA <- RLdata500[s1,]


################################ enlever les doublures

pairs1=compare.dedup(datasetB,identity=c(1:nrow(datasetB)),
                     blockfld=list(c(5,6),c(6,7),c(5,7)))#prendre comme variable de blockline list
summary(pairs1)
pairs1=emWeights(pairs1) #calculer le poids de record linkage par EM
pairs1=emClassify(pairs1, threshold.upper=20, threshold.lower=-5)
summary(pairs1)
possibles <- getPairs(pairs1, show="possible")
linke=getPairs(pairs1,show="links", single.rows=TRUE)


pairs2=compare.dedup(datasetA,identity=c(1:100),
                     blockfld=c(1,5))#prendre comme variable de blockline list
summary(pairs2)
pairs2=emWeights(pairs2) #calculer le poids de record linkage par EM
pairs2=emClassify(pairs2, threshold.upper=20, threshold.lower=-5)
summary(pairs2)
possibles <- getPairs(pairs2, show="possible")
linke=getPairs(pairs2,show="links", single.rows=TRUE)


