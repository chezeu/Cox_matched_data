library(simsalapar)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)
##############################################################

#setwd("C:\\Users\\Admin\\Dropbox\\R_program\\First_paper\\Continuous\\Results")
#setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")
setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

#########################################################################
load("0623_100s_res_exp_exp1.RData")
res <- res_exp_exp1
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods","error","lambdaE","n.sim")

dimnames(val)[[1]]=list("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                        "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                        "TPR-5 FS", "PPV-5 FS","converge FS" ,
                        "time gamma", "time FS3", "time FS")
df <- array2df(val)

df$methods = factor(df$methods)

df$prettyerror = factor(df$error, c("0.1", "0.2", "0.3"), c("e^k==0.1", "e^k==0.2", "e^k==0.3"))
df$prettylambdaE = factor(df$lambdaE, c("0.5000000","0.3333333", "0.2500000"), c("lambda[e]^k == 1/2", "lambda[e]^k == 1/3", "lambda[e]^k == 1/4"))
df$Methods = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                                  "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                                  "TPR-5 FS", "PPV-5 FS","converge FS",
                                  "time gamma", "time FS3", "time FS"), 
                    c("FS-HGa","FS-HGa","FS-HGa",
                      "FS3","FS3","FS3",
                      "FS","FS","FS",
                      "FS-HGa", "FS3", "FS"))
df$Criteria = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma",
                                   "TPR-5 FS3", "PPV-5 FS3",
                                   "TPR-5 FS", "PPV-5 FS"), 
                     c("TPR","PPV", "TPR","PPV","TPR","PPV"))


df_TPRPPV = subset(df, df$methods == "TPR-5 gamma"|df$methods == "PPV-5 gamma"|
                     df$methods == "TPR-5 FS3"|df$methods == "PPV-5 FS3"|
                     df$methods == "TPR-5 FS"|df$methods == "PPV-5 FS")

mean_TPRPPV <- df_TPRPPV %>% group_by(prettylambdaE, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV

df_time = subset(df, df$methods == "time gamma"|df$methods == "time FS3"|df$methods == "time FS")

mean_time <- df_time %>% group_by(prettyerror, prettylambdaE, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_time


df_converge = subset(df, df$methods == "converge gamma"|
                       df$methods == "converge FS3"|
                       df$methods == "converge FS")

mean_converge <- df_converge %>% group_by(prettyerror, prettylambdaE,Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)


sum_converge <- df_converge %>% group_by(n.sim, lambdaE, error) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum_converge$all_converge <- (sum_converge$sum==3)

df_TPRPPV$all_converge <- rep(sum_converge$all_converge, each = 6)
df_TPRPPV_converge = df_TPRPPV[df_TPRPPV$all_converge == TRUE,]
df_TPRPPV_NOconverge = df_TPRPPV[df_TPRPPV$all_converge == FALSE,]

mean_TPRPPV_converge <- df_TPRPPV_converge %>% group_by(prettylambdaE, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_converge

mean_TPRPPV_NOconverge <- df_TPRPPV_NOconverge %>% group_by(prettylambdaE, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_NOconverge


df_time$all_converge <-  rep(sum_converge$all_converge, each = 3)
df_time_converge <- df_time[df_time$all_converge == TRUE,]
mean_time_converge <- df_time_converge %>% group_by(prettyerror, prettylambdaE, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_time_converge

#library(xtable)
#print(xtable(mean_converge, type = "latex",digits=c(0,0,0,0,0,3)),include.rownames = FALSE, file = "cont_converge.tex")











#########################################################################
load("0406_1000s_res_exp_exp1.RData")
res <- res_exp_exp1
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods","error","lambdaE","n.sim")

dimnames(val)[[1]]=list("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                        "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                        "TPR-5 FS", "PPV-5 FS","converge FS" )
df <- array2df(val)

df$methods = factor(df$methods)

df$prettyerror = factor(df$error, c("0.1", "0.2", "0.3"), c("e^k==0.1", "e^k==0.2", "e^k==0.3"))
df$prettylambdaE = factor(df$lambdaE, c("0.5000000","0.3333333", "0.2500000"), c("lambda[e]^k == 1/2", "lambda[e]^k == 1/3", "lambda[e]^k == 1/4"))
df$Methods = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                                      "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                                      "TPR-5 FS", "PPV-5 FS","converge FS"), 
                              c("FS-HGa","FS-HGa","FS-HGa",
                                "FS3","FS3","FS3",
                                "FS","FS","FS"))
df$Criteria = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma",
                                  "TPR-5 FS3", "PPV-5 FS3",
                                  "TPR-5 FS", "PPV-5 FS"), 
                    c("TPR","PPV", "TPR","PPV","TPR","PPV"))


df_TPRPPV = subset(df, df$methods == "TPR-5 gamma"|df$methods == "PPV-5 gamma"|
                   df$methods == "TPR-5 FS3"|df$methods == "PPV-5 FS3"|
                  df$methods == "TPR-5 FS"|df$methods == "PPV-5 FS")

mean_TPRPPV <- df_TPRPPV %>% group_by(prettylambdaE, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV

df_converge = subset(df, df$methods == "converge gamma"|
                     df$methods == "converge FS3"|
                     df$methods == "converge FS")

mean_converge <- df_converge %>% group_by(prettyerror, prettylambdaE,Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)


sum_converge <- df_converge %>% group_by(n.sim, lambdaE, error) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum_converge$all_converge <- (sum_converge$sum==3)

df_TPRPPV$all_converge <- rep(sum_converge$all_converge, each = 6)
df_TPRPPV_converge = df_TPRPPV[df_TPRPPV$all_converge == TRUE,]
df_TPRPPV_NOconverge = df_TPRPPV[df_TPRPPV$all_converge == FALSE,]

mean_TPRPPV_converge <- df_TPRPPV_converge %>% group_by(prettylambdaE, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_converge

mean_TPRPPV_NOconverge <- df_TPRPPV_NOconverge %>% group_by(prettylambdaE, error, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_NOconverge


#library(xtable)
#print(xtable(mean_converge, type = "latex",digits=c(0,0,0,0,0,3)),include.rownames = FALSE, file = "cont_converge.tex")

save(df, file = "0406_1000s_res_exp_exp1_df.RData")
save(df_TPRPPV,file = "0406_1000s_res_exp_exp1_df_TPRPPV.RData")
save(df_TPRPPV_converge,file = "0406_1000s_res_exp_exp1_df_TPRPPV_converge.RData")
save(mean_TPRPPV,file = "0406_1000s_res_exp_exp1_mean_TPRPPV.RData")
save(mean_TPRPPV_converge,file = "0406_1000s_res_exp_exp1_mean_TPRPPV_converge.RData")
save(df_converge,file = "0406_1000s_res_exp_exp1_df_converge.RData")
save(mean_converge, file = "0406_1000s_res_exp_exp1_mean_converge.RData")



##############################################################

library(simsalapar)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)

#setwd("C:\\Users\\Admin\\Dropbox\\R_program\\First_paper\\Continuous\\Results")
setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

load("0421_1000s_res_exp_exp2.RData")
res <- res_exp_exp2
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods","lambdaK","n.sim")

dimnames(val)[[1]]=list("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                        "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                        "TPR-5 FS", "PPV-5 FS","converge FS" )
df <- array2df(val)

df$methods = factor(df$methods)


df$Methods = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                                  "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                                  "TPR-5 FS", "PPV-5 FS","converge FS"), 
                    c("FS-HGa","FS-HGa","FS-HGa",
                      "FS3","FS3","FS3",
                      "FS","FS","FS"))
df$Criteria = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma",
                                   "TPR-5 FS3", "PPV-5 FS3",
                                   "TPR-5 FS", "PPV-5 FS"), 
                     c("TPR","PPV", "TPR","PPV","TPR","PPV"))


df_TPRPPV = subset(df, df$methods == "TPR-5 gamma"|df$methods == "PPV-5 gamma"|
                     df$methods == "TPR-5 FS3"|df$methods == "PPV-5 FS3"|
                     df$methods == "TPR-5 FS"|df$methods == "PPV-5 FS")

df_converge = subset(df, df$methods == "converge gamma"|
                       df$methods == "converge FS3"|
                       df$methods == "converge FS")



mean_converge <- df_converge %>% group_by(lambdaK,Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)

sum_converge <- df_converge %>% group_by(n.sim, lambdaK) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum_converge$all_converge <- (sum_converge$sum==3)

df_TPRPPV$all_converge <- rep(sum_converge$all_converge, each = 6)
df_TPRPPV_converge = df_TPRPPV[df_TPRPPV$all_converge == TRUE,]
df_TPRPPV_NOconverge = df_TPRPPV[df_TPRPPV$all_converge == FALSE,]

mean_TPRPPV_converge <- df_TPRPPV_converge %>% group_by(lambdaK, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_converge


save(df, file = "0421_1000s_res_exp_exp2_df.RData")
save(df_TPRPPV_converge,file = "0421_1000s_res_exp_exp2_df_TPRPPV_converge.RData")
save(mean_TPRPPV_converge,file = "0421_1000s_res_exp_exp2_mean_TPRPPV_converge.RData")
save(df_converge,file = "0421_1000s_res_exp_exp2_df_converge.RData")
save(mean_converge, file = "0421_1000s_res_exp_exp2_mean_converge.RData")

##########################################################
##############################################################

library(simsalapar)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)

#setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Continuous\\Results")
setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

load("0422_1000s_res_exp_exp3.RData")
res <- res_exp_exp3
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods","nA","n.sim")

dimnames(val)[[1]]=list("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                        "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                        "TPR-5 FS", "PPV-5 FS","converge FS" )
df <- array2df(val)

df$methods = factor(df$methods)


df$Methods = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                                  "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                                  "TPR-5 FS", "PPV-5 FS","converge FS"), 
                    c("FS-HGa","FS-HGa","FS-HGa",
                      "FS3","FS3","FS3",
                      "FS","FS","FS"))
df$Criteria = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma",
                                   "TPR-5 FS3", "PPV-5 FS3",
                                   "TPR-5 FS", "PPV-5 FS"), 
                     c("TPR","PPV", "TPR","PPV","TPR","PPV"))


df_TPRPPV = subset(df, df$methods == "TPR-5 gamma"|df$methods == "PPV-5 gamma"|
                     df$methods == "TPR-5 FS3"|df$methods == "PPV-5 FS3"|
                     df$methods == "TPR-5 FS"|df$methods == "PPV-5 FS")

df_converge = subset(df, df$methods == "converge gamma"|
                       df$methods == "converge FS3"|
                       df$methods == "converge FS")



mean_converge <- df_converge %>% group_by(nA,Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)

sum_converge <- df_converge %>% group_by(n.sim, nA) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum_converge$all_converge <- (sum_converge$sum==3)

df_TPRPPV$all_converge <- rep(sum_converge$all_converge, each = 6)
df_TPRPPV_converge = df_TPRPPV[df_TPRPPV$all_converge == TRUE,]
df_TPRPPV_NOconverge = df_TPRPPV[df_TPRPPV$all_converge == FALSE,]

mean_TPRPPV_converge <- df_TPRPPV_converge %>% group_by(nA, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_converge

save(df, file = "0422_1000s_res_exp_exp3_df.RData")
save(df_TPRPPV_converge,file = "0422_1000s_res_exp_exp3_df_TPRPPV_converge.RData")
save(mean_TPRPPV_converge,file = "0422_1000s_res_exp_exp3_mean_TPRPPV_converge.RData")
save(df_converge,file = "0422_1000s_res_exp_exp3_df_converge.RData")
save(mean_converge, file = "0408_1000s_res_exp_exp3_mean_converge.RData")

##########################################################
##############################################################

library(simsalapar)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)

#setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Continuous\\Results")
setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

load("0420_1000s_res_exp_exp4.RData")
res <- res_exp_exp4
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods","K","lambdaK","n.sim")

dimnames(val)[[1]]=list("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                        "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                        "TPR-5 FS", "PPV-5 FS","converge FS" )
df <- array2df(val)

### 
df1 <- df[as.numeric(df$lambdaK) == 1 , ]
df2 <- df[as.numeric(df$lambdaK) == 2 , ]


####
df  = df1
df$methods = factor(df$methods)


df$Methods = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                                  "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                                  "TPR-5 FS", "PPV-5 FS","converge FS"), 
                    c("FS-HGa","FS-HGa","FS-HGa",
                      "FS3","FS3","FS3",
                      "FS","FS","FS"))
df$Criteria = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma",
                                   "TPR-5 FS3", "PPV-5 FS3",
                                   "TPR-5 FS", "PPV-5 FS"), 
                     c("TPR","PPV", "TPR","PPV","TPR","PPV"))


df_TPRPPV = subset(df, df$methods == "TPR-5 gamma"|df$methods == "PPV-5 gamma"|
                     df$methods == "TPR-5 FS3"|df$methods == "PPV-5 FS3"|
                     df$methods == "TPR-5 FS"|df$methods == "PPV-5 FS")

df_converge = subset(df, df$methods == "converge gamma"|
                       df$methods == "converge FS3"|
                       df$methods == "converge FS")



mean_converge <- df_converge %>% group_by(K,Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)

sum_converge <- df_converge %>% group_by(n.sim, K) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum_converge$all_converge <- (sum_converge$sum==3)

df_TPRPPV$all_converge <- rep(sum_converge$all_converge, each = 6)
df_TPRPPV_converge = df_TPRPPV[df_TPRPPV$all_converge == TRUE,]
df_TPRPPV_NOconverge = df_TPRPPV[df_TPRPPV$all_converge == FALSE,]

mean_TPRPPV_converge <- df_TPRPPV_converge %>% group_by(K, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_converge


save(df, file = "0420_1000s_res_exp_exp4_df.RData")
save(df_TPRPPV,file = "0420_1000s_res_exp_exp4_df_TPRPPV.RData")
save(df_TPRPPV_converge,file = "0420_1000s_res_exp_exp4_df_TPRPPV_converge.RData")
save(mean_TPRPPV_converge,file = "0420_1000s_res_exp_exp4_mean_TPRPPV_converge.RData")
save(df_converge,file = "0420_1000s_rres_exp_exp4_df_converge.RData")
save(mean_converge, file = "0420_1000s_rres_exp_exp4_mean_converge.RData")

##########################################################
##############################################################

#setwd("C:\\Users\\Admin\\Dropbox\\R_program\\First_paper\\Continuous\\Results")
setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

load("0410_1000s_res_unif_exp.RData")
res <- res_unif_exp
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods","error","lambdaE","n.sim")

dimnames(val)[[1]]=list("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                        "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                        "TPR-5 FS", "PPV-5 FS","converge FS" )
df <- array2df(val)

df$methods = factor(df$methods)

df$prettyerror = factor(df$error, c("0.1", "0.2", "0.3"), c("e==0.1", "e==0.2", "e==0.3"))
df$prettylambdaE = factor(df$lambdaE, c("0.5000000","0.3333333","0.2500000"), c("lambda[e] == 1/2", "lambda[e] == 1/3", "lambda[e] == 1/4"))
df$Methods = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                                  "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                                  "TPR-5 FS", "PPV-5 FS","converge FS"), 
                    c("FS-HGa","FS-HGa","FS-HGa",
                      "FS3","FS3","FS3",
                      "FS","FS","FS"))
df$Criteria = factor(df$methods, c("TPR-5 gamma", "PPV-5 gamma",
                                   "TPR-5 FS3", "PPV-5 FS3",
                                   "TPR-5 FS", "PPV-5 FS"), 
                     c("TPR","PPV", "TPR","PPV","TPR","PPV"))


df_TPRPPV = subset(df, df$methods == "TPR-5 gamma"|df$methods == "PPV-5 gamma"|
                     df$methods == "TPR-5 FS3"|df$methods == "PPV-5 FS3"|
                     df$methods == "TPR-5 FS"|df$methods == "PPV-5 FS")

df_converge = subset(df, df$methods == "converge gamma"|
                       df$methods == "converge FS3"|
                       df$methods == "converge FS")



mean_converge <- df_converge %>% group_by(prettyerror, prettylambdaE,Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)

sum_converge <- df_converge %>% group_by(n.sim, error, lambdaE) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum_converge$all_converge <- (sum_converge$sum==3)

df_TPRPPV$all_converge <- rep(sum_converge$all_converge, each = 6)
df_TPRPPV_converge = df_TPRPPV[df_TPRPPV$all_converge == TRUE,]


mean_TPRPPV_converge <- df_TPRPPV_converge %>% group_by(error, prettylambdaE, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_converge


#save(df, file = "df.RData")
save(df_TPRPPV_converge,file = "0410_1000s_res_unif_exp_df_TPRPPV_converge.RData")
save(mean_TPRPPV_converge,file = "0410_1000s_res_unif_exp_mean_TPRPPV_converge.RData")
#save(df_converge,file = "df_converge.RData")
#save(mean_converge, file = "mean_converge.RData")

##########################################################
##############################################################

#setwd("C:\\Users\\Admin\\Dropbox\\R_program\\First_paper\\Continuous\\Results")
setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

load("0329_1000s_res_exp_norm.RData")
res <- res_exp_norm
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods","error","meanE","n.sim")

dimnames(val)[[1]]=list("TPR-5 gamma2", "PPV-5 gamma2","converge gamma2" ,
                        "TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                        "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                        "TPR-5 FS", "PPV-5 FS","converge FS" )
df <- array2df(val)

df$methods = factor(df$methods)

df$prettyerror = factor(df$error, c("0.1", "0.2", "0.3"), c("e==0.1", "e==0.2", "e==0.3"))
df$prettymeanE = factor(df$meanE, c("1", "2","3"), c("mu[e] == 1", "mu[e] == 2", "mu[e] == 3"))
df$Methods = factor(df$methods, c("TPR-5 gamma2", "PPV-5 gamma2","converge gamma2" ,
                                  "TPR-5 gamma", "PPV-5 gamma","converge gamma" ,
                                  "TPR-5 FS3", "PPV-5 FS3","converge FS3" ,
                                  "TPR-5 FS", "PPV-5 FS","converge FS"), 
                    c("FS-HGa2","FS-HGa2","FS-HGa2",
                      "FS-HGa","FS-HGa","FS-HGa",
                      "FS3","FS3","FS3",
                      "FS","FS","FS"))
df$Criteria = factor(df$methods, c("TPR-5 gamma2", "PPV-5 gamma2",
                                   "TPR-5 gamma", "PPV-5 gamma",
                                   "TPR-5 FS3", "PPV-5 FS3",
                                   "TPR-5 FS", "PPV-5 FS"), 
                     c("TPR","PPV","TPR","PPV", "TPR","PPV","TPR","PPV"))


df_TPRPPV = subset(df, df$methods == "TPR-5 gamma2"|df$methods == "PPV-5 gamma2"|
                     df$methods == "TPR-5 gamma"|df$methods == "PPV-5 gamma"|
                     df$methods == "TPR-5 FS3"|df$methods == "PPV-5 FS3"|
                     df$methods == "TPR-5 FS"|df$methods == "PPV-5 FS")

df_converge = subset(df, df$methods == "converge gamma2"|
                       df$methods == "converge gamma"|
                       df$methods == "converge FS3"|
                       df$methods == "converge FS")


mean_converge <- df_converge %>% group_by(prettyerror, prettymeanE,Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)

sum_converge <- df_converge %>% group_by(n.sim, prettyerror, prettymeanE) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum_converge$all_converge <- (sum_converge$sum==4)

df_TPRPPV$all_converge <- rep(sum_converge$all_converge, each = 8)
df_TPRPPV_converge = df_TPRPPV[df_TPRPPV$all_converge == TRUE,]


mean_TPRPPV_converge <- df_TPRPPV_converge %>% group_by(error, prettymeanE, Criteria, Methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_TPRPPV_converge


#save(df, file = "df.RData")
save(df_TPRPPV_converge,file = "0329_1000s_res_exp_norm_df_TPRPPV_converge.RData")
save(mean_TPRPPV_converge,file = "0329_1000s_res_exp_norm_mean_TPRPPV_converge.RData")
#save(df_converge,file = "df_converge.RData")
#save(mean_converge, file = "mean_converge.RData")


