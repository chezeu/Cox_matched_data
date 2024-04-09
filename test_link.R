library(RecordLinkage)
library(plyr)
library(dplyr)

data("RLdata500")
pairs=compare.dedup(RLdata500,identity=identity.RLdata500,
                    blockfld=list(c(5,6),c(6,7),c(5,7)))#prendre comme variable de blockline list
summary(pairs)
pairs=emWeights(pairs) #calculer le poids de record linkage par EM
hist(pairs$Wdata, plot=FALSE)
getPairs(pairs,30,20)
pairs=emClassify(pairs, threshold.upper=24, threshold.lower=-7)
summary(pairs)
possibles <- getPairs(pairs, show="possible")
links=getPairs(pairs,show="links", single.rows=TRUE)
link_ids <- links[, c("id1", "id2")]
link_ids

data <- RLdata500[-link_ids$id1,]
paire=compare.dedup(data,identity=identity.RLdata500[-link_ids$id1],
                    blockfld=list(c(5,6),c(6,7),c(5,7)))#prendre comme variable de blockline list
summary(paire)
paire=emWeights(paire) #calculer le poids de record linkage par EM

paire=emClassify(paire, threshold.upper=24, threshold.lower=-7)
summary(paire)
possibles <- getPairs(paire, show="possible")
linke=getPairs(paire,show="links", single.rows=TRUE)



#########################record linkage#####################################

showClass("RLBigData")
showClass("RLBigDataDedup")
showClass("RLBigDataLinkage")

data("RLdata500")
data("RLdata10000")


s1 <- 471:500
s2 <- sample(1:10000, 300)
identity2 <- c(identity.RLdata500[s1], rep(NaN, length(s2)))
dataset <- rbind(RLdata500[s1,], RLdata10000[s2,])
rpairs2 <- RLBigDataLinkage(RLdata500, dataset, 
                            identity1 = identity.RLdata500,
                            identity2 = identity2, phonetic = 1:4, 
                            exclude = "lname_c2")

match <- rpairs2@pairs$is_match
ismatch <- which(match[,3]==1)
result <- rpairs2@pairs[ismatch,]

rpairs <- epiWeights(rpairs2)
getPairs(rpairs, max.weight=0.5, filter.match="match")

####################################


data("RLdata500")
data("RLdata10000")

library(reclin2)
s1 <- 496:500
s2 <- sample(1:10000, 15)
dataset1 <- rbind( RLdata10000[s2,],RLdata500[s1,])
dataset1$id1 <- c(identity.RLdata10000[s2], identity.RLdata500[s1])

dataset2 <- RLdata500[476:500,]
dataset2$id2 <- identity.RLdata500[476:500]

attach(dataset1)
attach(dataset2)
#dataset1$id1 <- identity1
#enlever les doublures

pairs1=compare.dedup(dataset1,identity=id1,
                     blockfld=c(1,5))#prendre comme variable de blockline list
summary(pairs1)
pairs1=emWeights(pairs1) #calculer le poids de record linkage par EM
pairs1=emClassify(pairs1, threshold.upper=24, threshold.lower=-7)
summary(pairs1)
possibles <- getPairs(pairs1, show="possible")
linke=getPairs(pairs1,show="links", single.rows=TRUE)


pairs2=compare.dedup(dataset2,identity=id2,
                     blockfld=c(1,5))#prendre comme variable de blockline list
summary(pairs2)
pairs2=emWeights(pairs2) #calculer le poids de record linkage par EM
pairs2=emClassify(pairs2, threshold.upper=24, threshold.lower=-7)
summary(pairs2)
possibles <- getPairs(pairs2, show="possible")
linke=getPairs(pairs2,show="links", single.rows=TRUE)


#blocking
p <- pair_blocking (dataset1, dataset2,"fname_c1" , FALSE)

p <- print(p)

# show the pairs of comparison en laissant la variable "bd"
p <- compare_pairs ( p, on= c("by","fname_c2","lname_c1","lname_c2","bm"), inplace = TRUE, 
                     comparators = list(fname_c2 = jaro_winkler(),
                                        lname_c1 = jaro_winkler(), lname_c2  = jaro_winkler() ))
p <- print(p)

## classification

model <- problink_em ( ~ by +fname_c2 +lname_c1 +lname_c2 +bm , data=p)

# print(model)       
#summary(model)

p <- predict(model, p, type ="mpost", add = TRUE, binary = TRUE)
prob_match <- p[,c(1:2,8)]

# Select pairs with a mpost > 0.5
p <- select_threshold(p, "selected", "mpost", 0.5, inplace = TRUE)



p <- add_from_x(p,"id1", id_x="id1")
p <- add_from_y(p,"id2", id_y="id2")

#p$true <- p$id_x==p$id_y
#table(as.data.frame(p[c("true","mpost")]))

# Select pairs with a mpost > 0.5 and force one-to-one linkage
p <- select_n_to_m(p,"selected2", "mpost", 0.5)

##creat a linked data
linked_data <- link(p)

##################considerons le cas ou toute la base A est contenue dans la base B#######
library(RecordLinkage)
library(plyr)
library(dplyr)
library(reclin2)

data("RLdata10000")
s1 <- sample(1:30,15)
s2 <- 1:30
dataset1 <- RLdata10000[s1,]
dataset1 <- dataset1[,-7]# il manque la variable "bd"
dataset2 <- RLdata10000[s2,]
dataset1$id1 <- identity.RLdata10000[s1]
dataset2$id2 <- identity.RLdata10000[s2]
attach(dataset1)
attach(dataset2)

pairs1=compare.dedup(dataset1,identity=id1,
                     blockfld=c(1,5))#prendre comme variable de blockline list

pairs2=compare.dedup(dataset2,identity=id2,
                     blockfld=c(1,5))#prendre comme variable de blockline list


#blocking
p <- pair_blocking (dataset1, dataset2,"fname_c1" , FALSE)
p <- print(p)

# show the pairs of comparison en laissant la variable "bd"
p <- compare_pairs ( p, on= c("by","fname_c2","lname_c1","lname_c2","bm"), inplace = TRUE, 
                     comparators = list(fname_c2 = jaro_winkler(),
                                        lname_c1 = jaro_winkler(), lname_c2  = jaro_winkler() ))
p <- print(p)

## classification using em method

model <- problink_em ( ~ by +fname_c2 +lname_c1 +lname_c2 +bm , data=p)

p <- predict(model, p, type ="mpost", add = TRUE, binary = TRUE)

# Select pairs with a mpost > 0.5
p <- select_threshold(p, "selected", "mpost", 0.5, inplace = TRUE)

p
p <- add_from_x(p,"id1", id_x="id1")
p <- add_from_y(p,"id2", id_y="id2")
p$true <- p$id1==p$id2
prob_match <- which(p$selected==TRUE)
p2 <- p[prob_match,]

# Select pairs with a mpost > 0.5 and force one-to-one linkage
#p <- select_n_to_m(p,"selected2", "mpost", 0.5)

##creat a linked data
linked_data <- link(p2)

id_pred_data <- linked_data$.y
pred_data <- dataset2[id_pred_data,]
Z_pred <- pred_data$bd


