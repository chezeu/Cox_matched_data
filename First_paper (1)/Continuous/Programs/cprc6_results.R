library(simsalapar)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)
##############################################################

setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Continuous\\Results")
#setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

load("0406_1000s_res1_exp_exp_prc.RData")
res <- res1_exp_exp_prc
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods","lambdaE", "n.sim")

dimnames(val)[1]=list(c(c("HGa-converged",paste("HGa-TPR-", 1:20, sep = "")), c(paste("HGa-PPV-", 1:20, sep = "")),
                        c(paste("HGa_obs-TPR-", 1:20, sep = "")), c(paste("HGa_obs-PPV-", 1:20, sep = "")),
                        c("FS3-converged",paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")),
                        c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")), 
                        c("FS-converged",paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")),
                        c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = ""))))
df <- array2df(val)

# factor error
df1 <- df[as.numeric(df$lambdaE) == 1 , ]
df2 <- df[as.numeric(df$lambdaE) == 2 , ]
#### 
df = df2

df$methods = factor(df$methods)

df.TPRPPV = subset(df, df$methods %in% c(c(paste("HGa-TPR-", 1:20, sep = "")), c(paste("HGa-PPV-", 1:20, sep = "")),
                                         c(paste("HGa_obs-TPR-", 1:20, sep = "")), c(paste("HGa_obs-PPV-", 1:20, sep = "")),
                                         c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")),
                                         c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")), 
                                         c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")),
                                         c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = ""))))

df.converge = subset(df, df$methods == "HGa-converged"|
                       df$methods == "FS3-converged"|
                       df$methods == "FS-converged")

sum.converge <- df.converge %>% group_by(n.sim, lambdaE) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum.converge$all_converge <- (sum.converge$sum==3)

df.TPRPPV$all_converge <- rep(sum.converge$all_converge, each = 240)
df.TPRPPV_converge <-  df.TPRPPV[df.TPRPPV$all_converge == TRUE,]
df.TPRPPV_NOconverge <-  df.TPRPPV[df.TPRPPV$all_converge == FALSE,]



df.mean <- df.TPRPPV_converge %>% group_by(lambdaE, methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE),
  Sd = sd(value, na.rm =TRUE)
)

df.mean

df.mean$fullmethods = factor(df.mean$methods,c(c(paste("HGa-TPR-", 1:20, sep = "")), c(paste("HGa-PPV-", 1:20, sep = "")),
                                               c(paste("HGa_obs-TPR-", 1:20, sep = "")), c(paste("HGa_obs-PPV-", 1:20, sep = "")),
                                               c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")),
                                               c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")), 
                                               c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")),
                                               c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = ""))),
                             c( rep("HGa-est",20),  rep("HGa-est",20),
                                rep("HGa-obs",20),  rep("HGa-obs",20),
                                rep("FS3-est",20),  rep("FS3-est",20),
                                rep("FS3-obs",20),  rep("FS3-obs",20),
                                rep("FS-est",20),  rep("FS-est",20),
                                rep("FS-obs",20), rep("FS-obs",20)))

df.mean$Methods = factor(df.mean$methods,c(c(paste("HGa-TPR-", 1:20, sep = "")), c(paste("HGa-PPV-", 1:20, sep = "")),
                                           c(paste("HGa_obs-TPR-", 1:20, sep = "")), c(paste("HGa_obs-PPV-", 1:20, sep = "")),
                                           c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")),
                                           c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")), 
                                           c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")),
                                           c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = ""))),
                         c(rep("FS-HGa",80),  
                           rep("FS3",80),
                           rep("FS",80)))

df.mean$Criteria = factor(df.mean$methods,c(c(paste("HGa-TPR-", 1:20, sep = "")), c(paste("HGa-PPV-", 1:20, sep = "")),
                                            c(paste("HGa_obs-TPR-", 1:20, sep = "")), c(paste("HGa_obs-PPV-", 1:20, sep = "")),
                                            c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")),
                                            c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")), 
                                            c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")),
                                            c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = ""))),
                          c(rep("TPR",20), rep("PPV",20),
                            rep("TPR",20), rep("PPV",20),
                            rep("TPR",20), rep("PPV",20),
                            rep("TPR",20), rep("PPV",20),
                            rep("TPR",20), rep("PPV",20),
                            rep("TPR",20), rep("PPV",20)))


df.mean$Status = factor(df.mean$methods,c(c(paste("HGa-TPR-", 1:20, sep = "")), c(paste("HGa-PPV-", 1:20, sep = "")),
                                          c(paste("HGa_obs-TPR-", 1:20, sep = "")), c(paste("HGa_obs-PPV-", 1:20, sep = "")),
                                          c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")),
                                          c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")), 
                                          c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")),
                                          c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = ""))),
                        c(rep("estimated",40),
                          rep("observed",40),
                          rep("estimated",40),
                          rep("observed",40),
                          rep("estimated",40),
                          rep("observed",40)))

x = c(0,
      1e-100,1e-20,1e-5,1e-2,0.05,
      0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
      0.95,0.975,0.995,0.999,
      1)
df.mean$thresholds = rep(x,12)

df.mean_TPR_est = subset(df.mean, df.mean$methods %in% c(paste("HGa-TPR-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS3-TPR-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS-TPR-", 1:20, sep = "")))

df.mean_TPR_obs = subset(df.mean, df.mean$methods %in% c(paste("HGa_obs-TPR-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS3_obs-TPR-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS_obs-TPR-", 1:20, sep = "")))


df.mean_PPV_est = subset(df.mean, df.mean$methods %in% c(paste("HGa-PPV-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS3-PPV-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS-PPV-", 1:20, sep = "")))

df.mean_PPV_obs = subset(df.mean, df.mean$methods %in% c(paste("HGa_obs-PPV-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS3_obs-PPV-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS_obs-PPV-", 1:20, sep = "")))

df.mean_fscore <- df.mean_TPR
df.mean_fscore$Mean <- (df.mean_TPR$Mean + df.mean_PPV$Mean)/2
df.mean_fscore$Criteria <- "f-score"
df.mean_fscore_est <- df.mean_fscore[df.mean_fscore$Status == "estimated",-9]
############################### Plots
library(ggplot2)
library(ggpubr)



df.mean_TPR = as.data.frame(rbind(df.mean_TPR_est, df.mean_TPR_obs))
df.mean_PPV = as.data.frame(rbind(df.mean_PPV_est, df.mean_PPV_obs))
df.mean_TPRPPV_est <- as.data.frame(rbind(df.mean_TPR_est, df.mean_PPV_est))

df.mean_fscore <- df.mean_TPR
df.mean_fscore$Mean <- (df.mean_TPR$Mean + df.mean_PPV$Mean)/2

df_prc = df.mean_TPR
df_prc$TPR = df.mean_TPR$Mean
df_prc$PPV = df.mean_PPV$Mean
# ############################## obs TPR, PPV
# tpr_obs <- ggplot(data = as.data.frame(df.mean_TPR), aes(x = thresholds, y = Mean, color =Methods, shape = Criteria )) +
#   geom_point(size = 2)+
#   geom_line(aes(linetype = Status)) +
#   xlab("thresholds") +
#   ylab("") +
#   guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE), shape = guide_legend(nrow = 2, byrow = FALSE)) +
#   scale_linetype_manual( values = c("solid", "dashed")) +
#   #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
#   scale_shape_manual(values = c(17)) +
#   theme_bw() +
#   theme(legend.position="top")
# 
# ppv_obs <- ggplot(data = as.data.frame(df.mean_PPV), aes(x = thresholds, y = Mean, color =Methods, shape = Criteria )) +
#   geom_point(size = 2)+
#   geom_line(aes(linetype = Status)) +
#   xlab("thresholds") +
#   ylab("") +
#   guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE), shape = guide_legend(nrow = 2, byrow = FALSE)) +
#   scale_linetype_manual( values = c("solid", "dashed")) +
#   #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
#   scale_shape_manual("", values = c(19) ) +
#   theme_bw() +
#   theme(legend.position="top")
# 
# 
# #ggarrange(tpr_obs1, ppv_obs1, ncol = 2, nrow = 1)
# #ggsave("0222_tpr_ppv_obs.pdf", width = 8, height = 4)
# 
# tpr_ppv <- ggplot(data = as.data.frame(df.mean_TPRPPV_est), aes(x = thresholds, y = Mean, color = Methods, shape = Criteria )) +
#   geom_point(size = 2)+
#   geom_line(aes(linetype = Status)) +
#   xlab("thresholds") +
#   ylab("") +
#   guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE), shape = guide_legend(nrow = 2, byrow = FALSE)) +
#   scale_linetype_manual( values = c("solid", "dashed")) +
#   #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
#   scale_shape_manual(values = c(17,19)) +
#   theme_bw() +
#   theme(legend.position="right")
# 
# tpr_ppv_obs <- ggplot(data = as.data.frame(df.mean), aes(x = thresholds, y = Mean, color = Methods, shape = Criteria )) +
#   geom_point(size = 2)+
#   geom_line(aes(linetype = Status)) +
#   xlab("thresholds") +
#   ylab("") +
#   guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE), shape = guide_legend(nrow = 2, byrow = FALSE)) +
#   scale_linetype_manual( values = c("solid", "dashed")) +
#   #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
#   scale_shape_manual(values = c(17,19))+
#   theme_bw() +
#   theme(legend.position="right")
# 
# ####################### grid 
# 
# library(gridExtra)
# library(lemon)
# legend <- g_legend(tpr_ppv_obs + theme(legend.position='right'))
# 
# ggarrange(tpr_obs+theme(legend.position='hidden'),  ppv_obs+theme(legend.position='hidden'),
#           tpr_ppv+theme(legend.position='hidden'), legend,  labels = c("A", "B", "C"), ncol = 2, nrow = 2)
# 
# ggsave("0305_cont_curves.pdf", width = 8, height = 8)
# 
# ################ PRC plot ##########
# df_prc = df.mean_TPR
# df_prc$TPR = df.mean_TPR$Mean
# df_prc$PPV = df.mean_PPV$Mean
# 
# prc <- ggplot(data = df_prc, aes(x=TPR, y=PPV,  color = Methods)) +
#   #geom_point(size = 2)+
#   geom_line(aes(linetype = Status)) +
#   xlab("TPR") +
#   ylab("PPV") +
#   guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE)) +
#   scale_linetype_manual( values = c("solid", "dashed")) +
#   ylim(0,1)+
#   xlim(0,1)+
#   theme_bw()+
#   theme(legend.position="top")
# prc
# ggsave("0305_cont_prc.pdf", width = 5, height = 5)




################Another version
############################## obs TPR, PPV
tpr_obs <- ggplot(data = as.data.frame(df.mean_TPR), aes(x = thresholds, y = Mean, color =Methods)) +
  #geom_point(size = 2)+
  geom_line(aes(linetype = Status)) +
  xlab("thresholds") +
  ylab("TPR") +
  guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE)) +
  scale_linetype_manual( values = c("dashed", "solid")) +
  #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
  scale_shape_manual(values = c(17)) +
  theme_bw() +
  theme(legend.position="top")
tpr_obs
ppv_obs <- ggplot(data = as.data.frame(df.mean_PPV), aes(x = thresholds, y = Mean, color =Methods)) +
  #geom_point(size = 2)+
  geom_line(aes(linetype = Status)) +
  xlab("thresholds") +
  ylab("PPV") +
  guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE)) +
  scale_linetype_manual( values = c("dashed", "solid")) +
  #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
  scale_shape_manual("", values = c(19) ) +
  theme_bw() +
  theme(legend.position="top")
ppv_obs

fscore_obs <- ggplot(data = as.data.frame(df.mean_fscore), aes(x = thresholds, y = Mean, color =Methods)) +
  #geom_point(size = 2)+
  geom_line(aes(linetype = Status)) +
  xlab("thresholds") +
  ylab("f-score") +
  guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE)) +
  scale_linetype_manual( values = c("dashed", "solid")) +
  #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
  scale_shape_manual("", values = c(19) ) +
  theme_bw() +
  theme(legend.position="top")
fscore_obs


prc <- ggplot(data = df_prc, aes(x=TPR, y=PPV,  color = Methods)) +
  #geom_point(size = 2)+
  geom_line(aes(linetype = Status)) +
  xlab("TPR") +
  ylab("PPV") +
  guides(color = guide_legend(nrow = 3, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE)) +
  scale_linetype_manual( values = c("dashed", "solid")) +
  ylim(0,1)+
  xlim(0,1)+
  theme_bw()+
  theme(legend.position="right")
prc
ggsave("0701_1000s_cont_prc.pdf", width = 6, height = 4)
ggsave("0701_1000s_cont_prc.png", width = 6, height = 4)

library(gridExtra)
library(lemon)
legend <- g_legend(ppv_obs + theme(legend.position='right'))

ggarrange(tpr_obs+theme(legend.position='hidden'),  ppv_obs+theme(legend.position='hidden'),
          fscore_obs+theme(legend.position='hidden'), legend,  labels = c("A", "B", "C"), ncol = 2, nrow = 2)

ggsave("0701_1000s_cont_curves.pdf", width = 8, height = 8)



################Another Status for estimated only
tpr_obs <- ggplot(data = as.data.frame(df.mean_TPR_est), aes(x = thresholds, y = Mean, color =Methods)) +
  #geom_point(size = 2)+
  geom_line() +
  xlab("thresholds") +
  ylab("TPR") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE)) +
  #scale_linetype_manual( values = c("solid", "dashed")) +
  #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
  #scale_shape_manual(values = c(17)) +
  theme_bw() +
  theme(legend.position="top")
tpr_obs
ppv_obs <- ggplot(data = as.data.frame(df.mean_PPV_est), aes(x = thresholds, y = Mean, color =Methods)) +
  #geom_point(size = 2)+
  geom_line() +
  xlab("thresholds") +
  ylab("PPV") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE)) +
  #scale_linetype_manual( values = c("solid", "dashed")) +
  #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
  #scale_shape_manual("", values = c(19) ) +
  theme_bw() +
  theme(legend.position="top")
ppv_obs

fscore_obs <- ggplot(data = as.data.frame(df.mean_fscore_est), aes(x = thresholds, y = Mean, color =Methods)) +
  #geom_point(size = 2)+
  geom_line() +
  xlab("thresholds") +
  ylab("f-score") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE)) +
  #scale_linetype_manual( values = c("solid", "dashed")) +
  #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
  scale_shape_manual("", values = c(19) ) +
  theme_bw() +
  theme(legend.position="top")
fscore_obs

df_prc_est = df_prc[df_prc$Status=="estimated",-9]
prc <- ggplot(data = df_prc_est, aes(x=TPR, y=PPV,  color = Methods)) +
  #geom_point(size = 2)+
  geom_line() +
  xlab("TPR") +
  ylab("PPV") +
  guides(color = guide_legend(nrow = 1, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE)) +
  #scale_linetype_manual( values = c("solid", "dashed")) +
  ylim(0,1)+
  xlim(0,1)+
  theme_bw()+
  theme(legend.position="top")
prc


ggarrange(tpr_obs,  ppv_obs,
          fscore_obs , prc,  labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, common.legend = TRUE)

ggsave("0406_1000s_cont_curves_slides_estimated.pdf", width = 8, height = 8)










#################################### with HGa2 method ######################
library(tidyr)
library(tidyverse)
library(remotes)
##############################################################

#setwd("C:\\Users\\Admin\\Dropbox\\R_program\\First_paper\\Continuous\\Results")
setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

load("2502_res_exp_norm.RData")
res <- res_exp_norm
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods","error",  "n.sim")

dimnames(val)[1]=list(c("HGa2-TPR-oneone", "HGa2-PPV-oneone",c(paste("HGa2-TPR-", 1:20, sep = "")), c(paste("HGa2-PPV-", 1:20, sep = "")),
                        "HGa2_obs-TPR-oneone", "HGa2_obs-PPV-oneone",c(paste("HGa2_obs-TPR-", 1:20, sep = "")), c(paste("HGa2_obs-PPV-", 1:20, sep = "")),
                        "HGa-TPR-oneone", "HGa-PPV-oneone",c(paste("HGa-TPR-", 1:20, sep = "")), c(paste("HGa-PPV-", 1:20, sep = "")),
                        "HGa_obs-TPR-oneone", "HGa_obs-PPV-oneone",c(paste("HGa_obs-TPR-", 1:20, sep = "")), c(paste("HGa_obs-PPV-", 1:20, sep = "")),
                        "FS3-TPR-oneone", "FS3-PPV-oneone",c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")),
                        "FS3_obs-TPR-oneone", "FS3_obs-PPV-oneone",c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")), 
                        "FS-TPR-oneone", "FS-PPV-oneone",c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")),
                        "FS_obs-TPR-oneone", "FS_obs-PPV-oneone",c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = ""))))
df <- array2df(val)

# factor error
df1 <- df[as.numeric(df$error) == 1 , ]
df2 <- df[as.numeric(df$error) == 2 , ]
#### 
df = df1

df$methods = factor(df$methods)

df.mean <- df %>% group_by(error, methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE),
  Sd = sd(value, na.rm =TRUE)
)

df.mean

df.mean$fullmethods = factor(df.mean$methods,c("HGa2-TPR-oneone", "HGa2-PPV-oneone",c(paste("HGa2-TPR-", 1:20, sep = "")), c(paste("HGa2-PPV-", 1:20, sep = "")),
                                               "HGa2_obs-TPR-oneone", "HGa2_obs-PPV-oneone",c(paste("HGa2_obs-TPR-", 1:20, sep = "")), c(paste("HGa2_obs-PPV-", 1:20, sep = "")),
                                               "HGa-TPR-oneone", "HGa-PPV-oneone",c(paste("HGa-TPR-", 1:20, sep = "")), c(paste("HGa-PPV-", 1:20, sep = "")),
                                               "HGa_obs-TPR-oneone", "HGa_obs-PPV-oneone",c(paste("HGa_obs-TPR-", 1:20, sep = "")), c(paste("HGa_obs-PPV-", 1:20, sep = "")),
                                               "FS3-TPR-oneone", "FS3-PPV-oneone",c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")),
                                               "FS3_obs-TPR-oneone", "FS3_obs-PPV-oneone",c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")), 
                                               "FS-TPR-oneone", "FS-PPV-oneone",c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")),
                                               "FS_obs-TPR-oneone", "FS_obs-PPV-oneone",c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = ""))),
                             c("HGa2-est-20","HGa2-est-20", rep("HGa2-est",20),  rep("HGa2-est",20),
                               "HGa2-obs-20","HGa2-obs-20", rep("HGa2-obs",20),  rep("HGa2-obs",20),
                               "HGa-est-20","HGa-est-20", rep("HGa-est",20),  rep("HGa-est",20),
                               "HGa-obs-20","HGa-obs-20", rep("HGa-obs",20),  rep("HGa-obs",20),
                               "FS3-est-20","FS3-est-20", rep("FS3-est",20),  rep("FS3-est",20),
                               "FS3-obs-20","FS3-obs-20", rep("FS3-obs",20),  rep("FS3-obs",20),
                               "FS-est-20","FS-est-20", rep("FS-est",20),  rep("FS-est",20),
                               "FS-obs-20","FS-obs-20", rep("FS-obs",20), rep("FS-obs",20)))

df.mean$Methods = factor(df.mean$methods,c("HGa2-TPR-oneone", "HGa2-PPV-oneone",c(paste("HGa2-TPR-", 1:20, sep = "")), c(paste("HGa2-PPV-", 1:20, sep = "")),
                                           "HGa2_obs-TPR-oneone", "HGa2_obs-PPV-oneone",c(paste("HGa2_obs-TPR-", 1:20, sep = "")), c(paste("HGa2_obs-PPV-", 1:20, sep = "")),
                                           "HGa-TPR-oneone", "HGa-PPV-oneone",c(paste("HGa-TPR-", 1:20, sep = "")), c(paste("HGa-PPV-", 1:20, sep = "")),
                                           "HGa_obs-TPR-oneone", "HGa_obs-PPV-oneone",c(paste("HGa_obs-TPR-", 1:20, sep = "")), c(paste("HGa_obs-PPV-", 1:20, sep = "")),
                                           "FS3-TPR-oneone", "FS3-PPV-oneone",c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")),
                                           "FS3_obs-TPR-oneone", "FS3_obs-PPV-oneone",c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")), 
                                           "FS-TPR-oneone", "FS-PPV-oneone",c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")),
                                           "FS_obs-TPR-oneone", "FS_obs-PPV-oneone",c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = ""))),
                         c(rep("HGa2",84),rep("HGa",84),  
                           rep("FS3",84),
                           rep("FS",84)))

df.mean$Criteria = factor(df.mean$methods,c("HGa2-TPR-oneone", "HGa2-PPV-oneone",c(paste("HGa2-TPR-", 1:20, sep = "")), c(paste("HGa2-PPV-", 1:20, sep = "")),
                                            "HGa2_obs-TPR-oneone", "HGa2_obs-PPV-oneone",c(paste("HGa2_obs-TPR-", 1:20, sep = "")), c(paste("HGa2_obs-PPV-", 1:20, sep = "")),
                                            "HGa-TPR-oneone", "HGa-PPV-oneone",c(paste("HGa-TPR-", 1:20, sep = "")), c(paste("HGa-PPV-", 1:20, sep = "")),
                                            "HGa_obs-TPR-oneone", "HGa_obs-PPV-oneone",c(paste("HGa_obs-TPR-", 1:20, sep = "")), c(paste("HGa_obs-PPV-", 1:20, sep = "")),
                                            "FS3-TPR-oneone", "FS3-PPV-oneone",c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")),
                                            "FS3_obs-TPR-oneone", "FS3_obs-PPV-oneone",c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")), 
                                            "FS-TPR-oneone", "FS-PPV-oneone",c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")),
                                            "FS_obs-TPR-oneone", "FS_obs-PPV-oneone",c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = ""))),
                          c("TPR","PPV", rep("TPR",20), rep("PPV",20),
                            "TPR","PPV", rep("TPR",20), rep("PPV",20),
                            "TPR","PPV", rep("TPR",20), rep("PPV",20),
                            "TPR","PPV", rep("TPR",20), rep("PPV",20),
                            "TPR","PPV", rep("TPR",20), rep("PPV",20),
                            "TPR","PPV", rep("TPR",20), rep("PPV",20),
                            "TPR","PPV", rep("TPR",20), rep("PPV",20),
                            "TPR","PPV", rep("TPR",20), rep("PPV",20)))


df.mean$Status = factor(df.mean$methods,c("HGa2-TPR-oneone", "HGa2-PPV-oneone",c(paste("HGa2-TPR-", 1:20, sep = "")), c(paste("HGa2-PPV-", 1:20, sep = "")),
                                          "HGa2_obs-TPR-oneone", "HGa2_obs-PPV-oneone",c(paste("HGa2_obs-TPR-", 1:20, sep = "")), c(paste("HGa2_obs-PPV-", 1:20, sep = "")),
                                          "HGa-TPR-oneone", "HGa-PPV-oneone",c(paste("HGa-TPR-", 1:20, sep = "")), c(paste("HGa-PPV-", 1:20, sep = "")),
                                          "HGa_obs-TPR-oneone", "HGa_obs-PPV-oneone",c(paste("HGa_obs-TPR-", 1:20, sep = "")), c(paste("HGa_obs-PPV-", 1:20, sep = "")),
                                          "FS3-TPR-oneone", "FS3-PPV-oneone",c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")),
                                          "FS3_obs-TPR-oneone", "FS3_obs-PPV-oneone",c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")), 
                                          "FS-TPR-oneone", "FS-PPV-oneone",c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")),
                                          "FS_obs-TPR-oneone", "FS_obs-PPV-oneone",c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = ""))),
                        c(rep("estimated",42),
                          rep("observed",42),
                          rep("estimated",42),
                          rep("observed",42),
                          rep("estimated",42),
                          rep("observed",42),
                          rep("estimated",42),
                          rep("observed",42)))

x = c(0,1e-100,1e-50,1e-30,1e-20,1e-10,1e-5,1e-2,
                     0.2,0.3,0.4,0.5,0.6,0.7,0.8,
                     0.9,0.95,0.99,0.995,1)
df.mean_TPR_est = subset(df.mean, df.mean$methods %in% c(paste("HGa2-TPR-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("HGa-TPR-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS3-TPR-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS-TPR-", 1:20, sep = "")))
df.mean_TPR_est$thresholds = rep(x,4)

df.mean_TPR_obs = subset(df.mean, df.mean$methods %in% c(paste("HGa2_obs-TPR-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("HGa_obs-TPR-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS3_obs-TPR-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS_obs-TPR-", 1:20, sep = "")))
df.mean_TPR_obs$thresholds = rep(x,4)


df.mean_PPV_est = subset(df.mean, df.mean$methods %in% c(paste("HGa2-PPV-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("HGa-PPV-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS3-PPV-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS-PPV-", 1:20, sep = "")))
df.mean_PPV_est$thresholds = rep(x,4)

df.mean_PPV_obs = subset(df.mean, df.mean$methods %in% c(paste("HGa2_obs-PPV-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("HGa_obs-PPV-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS3_obs-PPV-", 1:20, sep = ""))|
                           df.mean$methods %in% c(paste("FS_obs-PPV-", 1:20, sep = "")))
df.mean_PPV_obs$thresholds = rep(x,4)


############################### Plots
library(ggplot2)
library(ggpubr)
library(rstatix)


df.mean_TPR = as.data.frame(rbind(df.mean_TPR_est, df.mean_TPR_obs))
df.mean_PPV = as.data.frame(rbind(df.mean_PPV_est, df.mean_PPV_obs))

############################## threshold - TPR
tpr_obs1 <- ggplot(data = as.data.frame(df.mean_TPR), aes(x = thresholds, y = Mean, color =Methods, shape = Methods )) +
  geom_point(size = 2)+
  geom_line(aes(linetype = Status)) +
  xlab("thresholds") +
  ylab("TPR") +
  guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE), shape = guide_legend(nrow = 2, byrow = FALSE)) +
  #scale_linetype_manual("", values = ltbis, breaks = names(ltbis), drop = FALSE) +
  #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
  #scale_shape_manual("", values = shapsbis, breaks = names(shapsbis), drop = FALSE) +
  theme_bw() +
  theme(legend.position="top")
tpr_obs1

tpr_obs2 <- ggplot(data = as.data.frame(df.mean_TPR), aes(x = thresholds, y = 1- Mean, color =Methods, shape = Methods )) +
  geom_point(size = 2)+
  geom_line(aes(linetype = Status)) +
  xlab("thresholds") +
  ylab("1-TPR") +
  guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE), shape = guide_legend(nrow = 2, byrow = FALSE)) +
  #scale_linetype_manual("", values = ltbis, breaks = names(ltbis), drop = FALSE) +
  #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
  #scale_shape_manual("", values = shapsbis, breaks = names(shapsbis), drop = FALSE) +
  theme_bw() +
  theme(legend.position="top")
tpr_obs2

ggarrange(tpr_obs1, tpr_obs2, ncol = 2, nrow = 1, common.legend = TRUE)
ggsave("0225_threshold_tpr.pdf", width = 8, height = 4)

############################## threshold - PPV
ppv_obs1 <- ggplot(data = as.data.frame(df.mean_PPV), aes(x = thresholds, y = Mean, color =Methods, shape = Methods )) +
  geom_point(size = 2)+
  geom_line(aes(linetype = Status)) +
  xlab("thresholds") +
  ylab("PPV") +
  guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE), shape = guide_legend(nrow = 2, byrow = FALSE)) +
  #scale_linetype_manual("", values = ltbis, breaks = names(ltbis), drop = FALSE) +
  #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
  #scale_shape_manual("", values = shapsbis, breaks = names(shapsbis), drop = FALSE) +
  theme_bw() +
  theme(legend.position="top")
ppv_obs1

ggsave("0225_threshold_ppv.pdf", width = 8, height = 8)
############################ threshold - TPR - PPV
df.mean_TPRPPV_est <- as.data.frame(rbind(df.mean_TPR_est, df.mean_PPV_est))
tpr_ppv <- ggplot(data = as.data.frame(df.mean_TPRPPV_est), aes(x = thresholds, y = Mean, color = Methods, shape = Methods )) +
  geom_point(size = 2)+
  geom_line(aes(linetype = Criteria)) +
  xlab("thresholds") +
  ylab("") +
  guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE), shape = guide_legend(nrow = 2, byrow = FALSE)) +
  #scale_linetype_manual("", values = ltbis, breaks = names(ltbis), drop = FALSE) +
  #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
  #scale_shape_manual("", values = shapsbis, breaks = names(shapsbis), drop = FALSE) +
  theme_bw() +
  theme(legend.position="top")
tpr_ppv

ggsave("0225_tpr_ppv_02.pdf", width = 8, height = 8)

################ PRC plot ##########

df_prc = df.mean_TPR
df_prc$TPR = df.mean_TPR$Mean
df_prc$PPV = df.mean_PPV$Mean

prc <- ggplot(data = df_prc, aes(x=TPR, y=PPV,  color = Methods, shape = Methods)) +
  geom_point(size = 2)+
  geom_line(aes(linetype = Status)) +
  xlab("TPR") +
  ylab("PPV") +
  guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE), shape = guide_legend(nrow = 2, byrow = FALSE)) +
  ylim(0,1)+
  theme_bw()+
  theme(legend.position="top")
prc
ggsave("0225_prc.pdf", width = 8, height = 8)

