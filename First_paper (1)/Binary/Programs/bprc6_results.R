library(simsalapar)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)
##############################################################
setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Binary\\Results")
#setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Binary\\Results")

load("0402_1000s_res1_binary_prc.RData")
res <- res1
val <- getArray(res)


dimnames(val)
names(dimnames(val))=c("methods", "K", "n.sim")

dimnames(val)[1]=list(c(c("FS-converged", paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")), 
                        c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = "")),
                        c("FS3-converged",paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")), 
                        c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")),
                        c("FS4-converged",paste("FS4-TPR-", 1:20, sep = "")), c(paste("FS4-PPV-", 1:20, sep = "")),
                        c(paste("FS4_obs-TPR-", 1:20, sep = "")), c(paste("FS4_obs-PPV-", 1:20, sep = "")),
                        c("Bayesian-converged",paste("Bayesian-TPR-", 1:20, sep = "")), c(paste("Bayesian-PPV-", 1:20, sep = ""))))
df <- array2df(val)

df1 <- df[as.numeric(df$K) == 1, ]
df2 <- df[as.numeric(df$K) == 2, ]

#### 
df = df1

df$methods = factor(df$methods)

df.TPRPPV = subset(df, df$methods %in% c(c( paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")), 
                                         c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = "")),
                                         c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")), 
                                         c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")),
                                         c(paste("FS4-TPR-", 1:20, sep = "")), c(paste("FS4-PPV-", 1:20, sep = "")),
                                         c(paste("FS4_obs-TPR-", 1:20, sep = "")), c(paste("FS4_obs-PPV-", 1:20, sep = "")),
                                         c(paste("Bayesian-TPR-", 1:20, sep = "")), c(paste("Bayesian-PPV-", 1:20, sep = ""))))


df.converge = subset(df, df$methods == "FS-converged"|
                       df$methods == "FS3-converged"|
                       df$methods == "FS4-converged"|
                       df$methods == "Bayesian-converged")

sum.converge <- df.converge %>% group_by(n.sim,K) %>% dplyr::summarize(
  B = n(),
  sum = sum(value, na.rm = TRUE)
)
sum.converge$all_converge <- (sum.converge$sum==4)

df.TPRPPV$all_converge <- rep(sum.converge$all_converge, each = 280)
df.TPRPPV_converge <-  df.TPRPPV[df.TPRPPV$all_converge == TRUE,]
df.TPRPPV_NOconverge <-  df.TPRPPV[df.TPRPPV$all_converge == FALSE,]

df.mean <- df.TPRPPV_converge %>% group_by(K, methods) %>% dplyr::summarize(
  B = n(),
  Mean = mean(value, na.rm = TRUE),
  Sd = sd(value, na.rm =TRUE)
)

df.mean



df.mean$fullmethods = factor(df.mean$methods, c(c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")), 
                                        c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = "")),
                                        c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")), 
                                       c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")),
                                       c(paste("FS4-TPR-", 1:20, sep = "")), c(paste("FS4-PPV-", 1:20, sep = "")),
                                       c(paste("FS4_obs-TPR-", 1:20, sep = "")), c(paste("FS4_obs-PPV-", 1:20, sep = "")),
                                     c(paste("Bayesian-TPR-", 1:20, sep = "")), c(paste("Bayesian-PPV-", 1:20, sep = ""))), 
                           c(rep("FS-est",20),  rep("FS-est",20),
                             rep("FS-obs",20),  rep("FS-obs",20),
                             rep("FS3-est",20),  rep("FS3-est",20),
                             rep("FS3-obs",20),  rep("FS3-obs",20),
                             rep("FS4-est",20),  rep("FS4-est",20),
                             rep("FS4-obs",20), rep("FS4-obs",20),
                             rep("Bayesian-est",20),  rep("Bayesian-est",20)))
df.mean$Methods = factor(df.mean$methods, c(c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")), 
                                            c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = "")),
                                            c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")), 
                                            c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")),
                                            c(paste("FS4-TPR-", 1:20, sep = "")), c(paste("FS4-PPV-", 1:20, sep = "")),
                                            c(paste("FS4_obs-TPR-", 1:20, sep = "")), c(paste("FS4_obs-PPV-", 1:20, sep = "")),
                                            c(paste("Bayesian-TPR-", 1:20, sep = "")), c(paste("Bayesian-PPV-", 1:20, sep = ""))), 
                         c(rep("FS",20),  rep("FS",20),
                            rep("FS",20),  rep("FS",20),
                            rep("FS3",20),  rep("FS3",20),
                            rep("FS3",20),  rep("FS3",20),
                           rep("FS4",20),  rep("FS4",20),
                            rep("FS4",20), rep("FS4",20),
                           rep("Bayesian",20),  rep("Bayesian",20)))

df.mean$Criteria = factor(df.mean$methods, c(c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")), 
                                             c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = "")),
                                             c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")), 
                                             c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")),
                                             c(paste("FS4-TPR-", 1:20, sep = "")), c(paste("FS4-PPV-", 1:20, sep = "")),
                                             c(paste("FS4_obs-TPR-", 1:20, sep = "")), c(paste("FS4_obs-PPV-", 1:20, sep = "")),
                                             c(paste("Bayesian-TPR-", 1:20, sep = "")), c(paste("Bayesian-PPV-", 1:20, sep = ""))), 
                          c( rep("TPR",20), rep("PPV",20),
                            rep("TPR",20), rep("PPV",20),
                             rep("TPR",20), rep("PPV",20),
                            rep("TPR",20), rep("PPV",20),
                             rep("TPR",20), rep("PPV",20),
                            rep("TPR",20), rep("PPV",20),
                             rep("TPR",20), rep("PPV",20)))

df.mean$Status= factor(df.mean$methods, c(c(paste("FS-TPR-", 1:20, sep = "")), c(paste("FS-PPV-", 1:20, sep = "")), 
                                          c(paste("FS_obs-TPR-", 1:20, sep = "")), c(paste("FS_obs-PPV-", 1:20, sep = "")),
                                          c(paste("FS3-TPR-", 1:20, sep = "")), c(paste("FS3-PPV-", 1:20, sep = "")), 
                                          c(paste("FS3_obs-TPR-", 1:20, sep = "")), c(paste("FS3_obs-PPV-", 1:20, sep = "")),
                                          c(paste("FS4-TPR-", 1:20, sep = "")), c(paste("FS4-PPV-", 1:20, sep = "")),
                                          c(paste("FS4_obs-TPR-", 1:20, sep = "")), c(paste("FS4_obs-PPV-", 1:20, sep = "")),
                                          c(paste("Bayesian-TPR-", 1:20, sep = "")), c(paste("Bayesian-PPV-", 1:20, sep = ""))), 
                          c( rep("estimated",40),
                             rep("observed",40),
                            rep("estimated",40),
                            rep("observed",40),
                            rep("estimated",40),
                            rep("observed",40),
                            rep("estimated",40)))

x = c(0,
      1e-100,1e-20,1e-5,1e-2,0.05,
      0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
      0.95,0.975,0.995,0.999,
      1)

df.mean$thresholds = rep(x,14)


df.mean_TPR_est = subset(df.mean, df.mean$methods %in% c(paste("FS-TPR-", 1:20, sep = ""))|
                       df.mean$methods %in% c(paste("FS3-TPR-", 1:20, sep = ""))|
                       df.mean$methods %in% c(paste("FS4-TPR-", 1:20, sep = ""))|
                       df.mean$methods %in% c(paste("Bayesian-TPR-", 1:20, sep = "")))


df.mean_TPR_obs = subset(df.mean, df.mean$methods %in% c(paste("FS_obs-TPR-", 1:20, sep = ""))|
                       df.mean$methods %in% c(paste("FS3_obs-TPR-", 1:20, sep = ""))|
                       df.mean$methods %in% c(paste("FS4_obs-TPR-", 1:20, sep = "")))


df.mean_PPV_est = subset(df.mean, df.mean$methods %in% c(paste("FS-PPV-", 1:20, sep = ""))|
                       df.mean$methods %in% c(paste("FS3-PPV-", 1:20, sep = ""))|
                       df.mean$methods %in% c(paste("FS4-PPV-", 1:20, sep = ""))|
                       df.mean$methods %in% c(paste("Bayesian-PPV-", 1:20, sep = "")))


df.mean_PPV_obs = subset(df.mean, df.mean$methods %in% c(paste("FS_obs-PPV-", 1:20, sep = ""))|
                       df.mean$methods %in% c(paste("FS3_obs-PPV-", 1:20, sep = ""))|
                       df.mean$methods %in% c(paste("FS4_obs-PPV-", 1:20, sep = "")))



############################### Plots
library(ggplot2)
library(ggpubr)
library(rstatix)


df.mean_TPR = as.data.frame(rbind(df.mean_TPR_est, df.mean_TPR_obs))
df.mean_PPV = as.data.frame(rbind(df.mean_PPV_est, df.mean_PPV_obs))
df.mean_TPRPPV_est <- as.data.frame(rbind(df.mean_TPR_est, df.mean_PPV_est))

df.mean_fscore <- df.mean_TPR
df.mean_fscore$Mean <- (df.mean_TPR$Mean + df.mean_PPV$Mean)/2
df.mean_fscore$Criteria <- "f-score"
df.mean_fscore_est <- df.mean_fscore[df.mean_fscore$Status == "estimated",-9]

df_prc = df.mean_TPR
df_prc$TPR = df.mean_TPR$Mean
df_prc$PPV = df.mean_PPV$Mean
df_prc_est = df_prc[df_prc$Status=="estimated",-9]
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
#             tpr_ppv+theme(legend.position='hidden'), legend,  labels = c("A", "B", "C"), ncol = 2, nrow = 2)
# 
# ggsave("0303_binary_curves.pdf", width = 8, height = 8)
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
# ggsave("0303_binary_prc.pdf", width = 5, height = 5)





################Another Status
############################## obs TPR, PPV
tpr_obs <- ggplot(data = as.data.frame(df.mean_TPR), aes(x = thresholds, y = Mean, color =Methods)) +
  #geom_point(size = 2)+
  geom_line(aes(linetype = Status)) +
  xlab("thresholds") +
  ylab("TPR") +
  guides(color = guide_legend(nrow = 2, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE)) +
  scale_linetype_manual( values = c("dashed", "solid")) +
  #scale_color_manual("", values = colsbis, breaks = names(colsbis), drop = FALSE) +
  #scale_shape_manual(values = c(17)) +
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
  #scale_shape_manual("", values = c(19) ) +
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
  guides(color = guide_legend(nrow = 4, byrow = FALSE), linetype = guide_legend(nrow = 2, byrow = FALSE)) +
  scale_linetype_manual( values = c("dashed", "solid")) +
  ylim(0,1)+
  xlim(0,1)+
  theme_bw()+
  theme(legend.position="right")
prc
ggsave("0701_1000s_binary_prc_converge.pdf", width = 6, height = 4)
ggsave("0701_1000s_binary_prc_converge.png", width = 6, height = 4)

library(gridExtra)
library(lemon)
legend <- g_legend(ppv_obs+ theme(legend.position='right'))

ggarrange(tpr_obs+theme(legend.position='hidden'),  ppv_obs+theme(legend.position='hidden'),
          fscore_obs+theme(legend.position='hidden'), legend,  labels = c("A", "B", "C"), ncol = 2, nrow = 2)

ggsave("0701_1000s_binary_curves_converge.pdf", width = 8, height = 8)







################Another Status for estimated only
############################## obs TPR, PPV
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

ggsave("0319_1000s_binary_curves_slides_estimated.pdf", width = 8, height = 8)












