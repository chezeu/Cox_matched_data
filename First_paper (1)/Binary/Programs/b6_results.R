library(simsalapar)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)

##############################################################

#setwd("C:\\Users\\Admin\\Dropbox\\R_program\\First_paper\\Binary\\Results")
#setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Binary\\Results")
setwd("D:/chezeu article/code-vanessa/First_paper (1)/Binary/Results")

load("0720_1000s_res1_binary.RData")
val <- getArray(res1)


dimnames(val)
names(dimnames(val))=c("methods","K","error","n.sim")

dimnames(val)[[1]]=list("TPR FS","PPV FS", "converge FS",
                        "TPR FS3","PPV FS3","converge FS3",
                        "TPR FS4","PPV FS4", "converge FS4",
                        "TPR Bayes","PPV Bayes","converge Bayes",
                        "time FS", "time FS3", "time FS4", "time Bayes")
df <- array2df(val)


df$methods = factor(df$methods)

df$prettyK = factor(df$K, c("30","40", "50"), c("K == 30", "K==40", "K==50"))
df$prettyerror = factor(df$error, c("0.02", "0.04","0.06"), c("e==0.02", "e==0.04", "e==0.06"))

df$Methods = factor(df$methods, c("TPR FS","PPV FS",  "converge FS", "TPR FS3","PPV FS3", "converge FS3",
                                  "TPR FS4","PPV FS4", "converge FS4", "TPR Bayes","PPV Bayes", "converge Bayes",
                                  "time FS", "time FS3", "time FS4", "time Bayes"), 
                    c("FS","FS","FS","FS3","FS3","FS3",
                      "FS4","FS4","FS4","Bayesian","Bayesian","Bayesian",
                      "FS", "FS3", "FS4", "Bayes"))

df$Criteria = factor(df$methods, c("TPR FS","PPV FS", "TPR FS3","PPV FS3",
                                   "TPR FS4","PPV FS4", "TPR Bayes","PPV Bayes"), 
                     c("TPR","PPV","TPR","PPV",
                       "TPR","PPV","TPR","PPV"))

df_TPRPPV = subset(df, df$methods == "TPR FS"|df$methods == "PPV FS"|
                     df$methods == "TPR FS3"|df$methods == "PPV FS3"|
                     df$methods == "TPR FS4"|df$methods == "PPV FS4"|
                     df$methods == "TPR Bayes"|df$methods == "PPV Bayes")

mean_TPRPPV <- df_TPRPPV %>% group_by(prettyK, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV

df_converge = subset(df, df$methods == "converge FS"|
                       df$methods == "converge FS3"|
                       df$methods == "converge FS4"|
                       df$methods == "converge Bayes")

sum_converge <- df_converge %>% group_by(n.sim,error,K) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum_converge$all_converge <- (sum_converge$sum==4)

df_TPRPPV$all_converge <- rep(sum_converge$all_converge, each = 8)
df_TPRPPV_converge = df_TPRPPV[df_TPRPPV$all_converge == TRUE,]
df_TPRPPV_NOconverge = df_TPRPPV[df_TPRPPV$all_converge == FALSE,]

mean_TPRPPV_converge <- df_TPRPPV_converge %>% group_by(prettyK, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_converge

mean_TPRPPV_NOconverge <- df_TPRPPV_NOconverge %>% group_by(prettyK, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_NOconverge

mean_converge <- df_converge %>% group_by(prettyK, prettyerror, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_converge


df_time <- subset(df, df$methods == "time FS"|df$methods == "time FS3"|df$methods == "time FS4"|df$methods == "time Bayes")
mean_time <- df_time%>% group_by(prettyK, error, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_time

df_time$converge <- df_converge$value
df_time_converge <- df_time[df_time$converge == 1,]
mean_time_converge <- df_time_converge %>% group_by(prettyK, error,  Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_time_converge
##############################################################
##############################################################

#setwd("C:\\Users\\Admin\\Dropbox\\R_program\\First_paper\\Binary\\Results")
#setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Binary\\Results")
setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Binary\\Results")

load("0308_1000s_res1_binary.RData")
val <- getArray(res1)


dimnames(val)
names(dimnames(val))=c("methods","K","error","n.sim")

dimnames(val)[[1]]=list("TPR FS","PPV FS", "converge FS",
                        "TPR FS3","PPV FS3","converge FS3",
                        "TPR FS4","PPV FS4", "converge FS4",
                        "TPR Bayes","PPV Bayes","converge Bayes")
df <- array2df(val)


df$methods = factor(df$methods)

df$prettyK = factor(df$K, c("30","40", "50"), c("K == 30", "K==40", "K==50"))
df$prettyerror = factor(df$error, c("0.02", "0.04","0.06"), c("e==0.02", "e==0.04", "e==0.06"))

df$Methods = factor(df$methods, c("TPR FS","PPV FS",  "converge FS", "TPR FS3","PPV FS3", "converge FS3",
                                        "TPR FS4","PPV FS4", "converge FS4", "TPR Bayes","PPV Bayes", "converge Bayes"), 
                                  c("FS","FS","FS","FS3","FS3","FS3",
                                    "FS4","FS4","FS4","Bayesian","Bayesian","Bayesian"))

df$Criteria = factor(df$methods, c("TPR FS","PPV FS", "TPR FS3","PPV FS3",
                                  "TPR FS4","PPV FS4", "TPR Bayes","PPV Bayes"), 
                          c("TPR","PPV","TPR","PPV",
                            "TPR","PPV","TPR","PPV"))

df_TPRPPV = subset(df, df$methods == "TPR FS"|df$methods == "PPV FS"|
                     df$methods == "TPR FS3"|df$methods == "PPV FS3"|
                     df$methods == "TPR FS4"|df$methods == "PPV FS4"|
                     df$methods == "TPR Bayes"|df$methods == "PPV Bayes")

mean_TPRPPV <- df_TPRPPV %>% group_by(prettyK, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV

df_converge = subset(df, df$methods == "converge FS"|
                  df$methods == "converge FS3"|
                  df$methods == "converge FS4"|
                  df$methods == "converge Bayes")

sum_converge <- df_converge %>% group_by(n.sim,error,K) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum_converge$all_converge <- (sum_converge$sum==4)

df_TPRPPV$all_converge <- rep(sum_converge$all_converge, each = 8)
df_TPRPPV_converge = df_TPRPPV[df_TPRPPV$all_converge == TRUE,]
df_TPRPPV_NOconverge = df_TPRPPV[df_TPRPPV$all_converge == FALSE,]

mean_TPRPPV_converge <- df_TPRPPV_converge %>% group_by(prettyK, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_converge

mean_TPRPPV_NOconverge <- df_TPRPPV_NOconverge %>% group_by(prettyK, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_NOconverge

mean_converge <- df_converge %>% group_by(prettyK, prettyerror, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_converge

library(xtable)
print(xtable(mean_converge, type = "latex"), file = "binary_converge.tex")

save(df, file = "0308_1000s_res1_df.RData")
save(df_TPRPPV,file = "0308_1000s_res1_df_TPRPPV.RData")
save(df_TPRPPV_converge,file = "0308_1000s_res1_df_TPRPPV_converge.RData")
save(mean_TPRPPV,file = "0308_1000s_res1_mean_TPRPPV.RData")
save(mean_TPRPPV_converge,file = "0308_1000s_res1_mean_TPRPPV_converge.RData")
save(df_converge,file = "0308_1000s_res1_df_converge.RData")
save(mean_converge, file = "0308_1000s_res1_mean_converge.RData")

######################################################
####
######################################################
library(simsalapar)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)

setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Binary\\Results")
#setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Binary\\Results")

load("0324_1000s_res2_binary.RData")
val <- getArray(res2)


dimnames(val)
names(dimnames(val))=c("methods","prev","n.sim")

dimnames(val)[[1]]=list("TPR FS","PPV FS", "converge FS",
                        "TPR FS3","PPV FS3","converge FS3",
                        "TPR FS4","PPV FS4", "converge FS4",
                        "TPR Bayes","PPV Bayes","converge Bayes")
df <- array2df(val)


df$methods = factor(df$methods)

df$prettyprev = factor(df$prev, c("0.1","0.2", "0.3"), c("p^k == 0.1", "p^k == 0.2", "p^k == 0.3"))

df$Methods = factor(df$methods, c("TPR FS","PPV FS",  "converge FS", "TPR FS3","PPV FS3", "converge FS3",
                                  "TPR FS4","PPV FS4", "converge FS4", "TPR Bayes","PPV Bayes", "converge Bayes"), 
                    c("FS","FS","FS","FS3","FS3","FS3",
                      "FS4","FS4","FS4","Bayesian","Bayesian","Bayesian"))

df$Criteria = factor(df$methods, c("TPR FS","PPV FS", "TPR FS3","PPV FS3",
                                   "TPR FS4","PPV FS4", "TPR Bayes","PPV Bayes"), 
                     c("TPR","PPV","TPR","PPV",
                       "TPR","PPV","TPR","PPV"))

df_TPRPPV = subset(df, df$methods == "TPR FS"|df$methods == "PPV FS"|
                     df$methods == "TPR FS3"|df$methods == "PPV FS3"|
                     df$methods == "TPR FS4"|df$methods == "PPV FS4"|
                     df$methods == "TPR Bayes"|df$methods == "PPV Bayes")

mean_TPRPPV <- df_TPRPPV %>% group_by(prev,Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV

df_converge = subset(df, df$methods == "converge FS"|
                       df$methods == "converge FS3"|
                       df$methods == "converge FS4"|
                       df$methods == "converge Bayes")

mean_converge <- df_converge %>% group_by(prettyprev, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_converge

sum_converge <- df_converge %>% group_by(n.sim,prettyprev) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum_converge$all_converge <- (sum_converge$sum==4)

df_TPRPPV$all_converge <- rep(sum_converge$all_converge, each = 8)

df_TPRPPV_converge = df_TPRPPV[df_TPRPPV$all_converge == TRUE,]
df_TPRPPV_NOconverge = df_TPRPPV[df_TPRPPV$all_converge == FALSE,]

mean_TPRPPV_converge <- df_TPRPPV_converge %>% group_by(prev, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_converge



save(df, file = "0324_1000s_res2_df.RData")
save(df_TPRPPV_converge,file = "0324_1000s_res2_df_TPRPPV_converge.RData")
save(mean_TPRPPV_converge,file = "0324_1000s_res2_mean_TPRPPV_converge.RData")
save(df_converge,file = "0324_1000s_res2_df_converge.RData")
save(mean_converge, file = "0324_1000s_res2_mean_converge.RData")

#####################################################################
##
####################################################################
library(simsalapar)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)

setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Binary\\Results")
#setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Binary\\Results")

load("0324_1000s_res3_binary.RData")
val <- getArray(res3)


dimnames(val)
names(dimnames(val))=c("methods","nA","n.sim")

dimnames(val)[[1]]=list("TPR FS","PPV FS", "converge FS",
                        "TPR FS3","PPV FS3","converge FS3",
                        "TPR FS4","PPV FS4", "converge FS4",
                        "TPR Bayes","PPV Bayes","converge Bayes")
df <- array2df(val)


df$methods = factor(df$methods)

df$prettynA = factor(df$nA, c("400","800", "1200"), c("n[A] == 400", "n[A] == 800", "n[A] == 1200"))

df$Methods = factor(df$methods, c("TPR FS","PPV FS",  "converge FS", "TPR FS3","PPV FS3", "converge FS3",
                                  "TPR FS4","PPV FS4", "converge FS4", "TPR Bayes","PPV Bayes", "converge Bayes"), 
                    c("FS","FS","FS","FS3","FS3","FS3",
                      "FS4","FS4","FS4","Bayesian","Bayesian","Bayesian"))

df$Criteria = factor(df$methods, c("TPR FS","PPV FS", "TPR FS3","PPV FS3",
                                   "TPR FS4","PPV FS4", "TPR Bayes","PPV Bayes"), 
                     c("TPR","PPV","TPR","PPV",
                       "TPR","PPV","TPR","PPV"))

df_TPRPPV = subset(df, df$methods == "TPR FS"|df$methods == "PPV FS"|
                     df$methods == "TPR FS3"|df$methods == "PPV FS3"|
                     df$methods == "TPR FS4"|df$methods == "PPV FS4"|
                     df$methods == "TPR Bayes"|df$methods == "PPV Bayes")

mean_TPRPPV <- df_TPRPPV %>% group_by(nA,Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV

df_converge = subset(df, df$methods == "converge FS"|
                       df$methods == "converge FS3"|
                       df$methods == "converge FS4"|
                       df$methods == "converge Bayes")


mean_converge <- df_converge %>% group_by(prettynA, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_converge

sum_converge <- df_converge %>% group_by(n.sim,nA) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum_converge$all_converge <- (sum_converge$sum==4)

df_TPRPPV$all_converge <- rep(sum_converge$all_converge, each = 8)

df_TPRPPV_converge = df_TPRPPV[df_TPRPPV$all_converge == TRUE,]
df_TPRPPV_NOconverge = df_TPRPPV[df_TPRPPV$all_converge == FALSE,]

mean_TPRPPV_converge <- df_TPRPPV_converge %>% group_by(nA, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_converge

save(df, file = "0324_1000s_res3_df.RData")
save(df_TPRPPV_converge,file = "0324_1000s_res3_df_TPRPPV_converge.RData")
save(mean_TPRPPV_converge,file = "0324_1000s_res3_mean_TPRPPV_converge.RData")
save(df_converge,file = "0324_1000s_res3_df_converge.RData")
save(mean_converge, file = "0324_1000s_res3_mean_converge.RData")

