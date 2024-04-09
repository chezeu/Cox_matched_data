library(simsalapar)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)

setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Mixed\\Results")

load("0624_1000s_resv1.RData")
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods", "error3","lambda3","n.sim")

dimnames(val)[[1]]=list("TPR-5 mixed", "PPV-5 mixed", "converge mixed", 
                        "TPR-5 binary", "PPV-5 binary", "converge binary", 
                        "time mixed", "time binary")
df <- array2df(val)

df$methods = factor(df$methods)
df$error3 =  factor(df$error3)
df$lambda3 = factor(df$lambda3)

df_converge = subset(df, df$methods == "converge mixed"|df$methods == "converge binary")
mean_converge <- df_converge %>% group_by(error3,lambda3,methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_converge

df_time = subset(df, df$methods == "time mixed"|df$methods == "time binary")
mean_time <- df_time %>% group_by(error3,lambda3,methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE), 
  Sd = sd(value, na.rm = TRUE)
)
mean_time


df_TPRPPV = subset(df, df$methods == "TPR-5 mixed"|df$methods == "PPV-5 mixed"|
                       df$methods == "TPR-5 binary"|df$methods == "PPV-5 binary")
mean_TPRPPV <- df_TPRPPV %>% group_by(error3,lambda3,methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE),
  Sd = sd(value, na.rm = TRUE)
)
mean_TPRPPV


##########################################################################
setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Mixed\\Results")

load("0708_1000s_resv1.RData")
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods", "error3","n.sim")

dimnames(val)[[1]]=list("TPR-5 mixed", "PPV-5 mixed", "converge mixed", 
                        "TPR-5 binary", "PPV-5 binary", "converge binary", 
                        "TPR binary2", "PPV binary2",
                        "time mixed", "time binary", "time binary2")
df <- array2df(val)

df$methods = factor(df$methods)
df$error3 =  factor(df$error3)

df_converge = subset(df, df$methods == "converge mixed"|df$methods == "converge binary")
mean_converge <- df_converge %>% group_by(error3,methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_converge

df_time = subset(df, df$methods == "time mixed"|df$methods == "time binary"|df$methods == "time binary2")
mean_time <- df_time %>% group_by(error3,methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE), 
  Sd = sd(value, na.rm = TRUE)
)
mean_time


df_TPRPPV = subset(df, df$methods == "TPR-5 mixed"|df$methods == "PPV-5 mixed"|
                     df$methods == "TPR-5 binary"|df$methods == "PPV-5 binary"|
                     df$methods == "TPR binary2"|df$methods == "PPV binary2")
mean_TPRPPV <- df_TPRPPV %>% group_by(error3,methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE),
  Sd = sd(value, na.rm = TRUE)
)
mean_TPRPPV


##########################################################################
setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Mixed\\Results")

load("1028_1000s_resv1.RData")
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods", "error3","n.sim")

dimnames(val)[[1]]=list("TPR-5 mixed", "PPV-5 mixed", "converge mixed", 
                        "TPR-5 binary", "PPV-5 binary", "converge binary", 
                        "TPR binary2", "PPV binary2",
                        "TPR binary3", "PPV binary3",
                        "time mixed", "time binary", "time binary2", "time binary3")
df <- array2df(val)

df$methods = factor(df$methods)
df$error3 =  factor(df$error3)

df_converge = subset(df, df$methods == "converge mixed"|df$methods == "converge binary")
mean_converge <- df_converge %>% group_by(error3,methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE)
)
mean_converge

df_time = subset(df, df$methods == "time mixed"|df$methods == "time binary"|df$methods == "time binary2"|df$methods == "time binary3")
mean_time <- df_time %>% group_by(error3,methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE), 
  Sd = sd(value, na.rm = TRUE)
)
mean_time


df_TPRPPV = subset(df, df$methods == "TPR-5 mixed"|df$methods == "PPV-5 mixed"|
                     df$methods == "TPR-5 binary"|df$methods == "PPV-5 binary"|
                     df$methods == "TPR binary2"|df$methods == "PPV binary2"|
                     df$methods == "TPR binary3"|df$methods == "PPV binary3")
mean_TPRPPV <- df_TPRPPV %>% group_by(error3,methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE),
  Sd = sd(value, na.rm = TRUE)
)
mean_TPRPPV
