library(ggplot2)
library(ggpubr)
#library(rstatix)

#setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")
setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

#########################################
### res_exp_exp1
########################################
load("0406_1000s_res_exp_exp1_df_TPRPPV.Rdata")

TPRPPV <- ggplot(data = as.data.frame(df_TPRPPV), aes(x=error, y=value, fill=Methods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid( prettylambdaE ~  Criteria , labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0.5,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")
TPRPPV
ggsave("0406_1000s_cont_boxplots_exp_exp1.pdf", width = 8, height = 8)

load("0406_1000s_res_exp_exp1_mean_TPRPPV.Rdata")

mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV), aes(x=as.numeric(levels(error))[error], y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid( prettylambdaE ~  Criteria , labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0.5,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")+
  scale_x_continuous(breaks=c(0.1,0.2,0.3))
mean_fig
ggsave("0406_1000s_cont_means_exp_exp1.pdf", width = 8, height = 8)


load("0406_1000s_res_exp_exp1_df_TPRPPV_converge.Rdata")
TPRPPV <- ggplot(data = as.data.frame(df_TPRPPV_converge), aes(x=error, y=value, fill=Methods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid( prettylambdaE ~  Criteria , labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0.5,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")
TPRPPV
ggsave("0406_1000s_cont_boxplots_exp_exp1_converge.pdf", width = 8, height = 8)

load("0406_1000s_res_exp_exp1_mean_TPRPPV_converge.Rdata")
mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV_converge), aes(x=as.numeric(levels(error))[error], y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid( prettylambdaE ~  Criteria , labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0.5,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")+
  scale_x_continuous(breaks=c(0.1,0.2,0.3))
mean_fig
ggsave("0406_1000s_cont_means_exp_exp1_converge.pdf", width = 8, height = 8)
ggsave("0406_1000s_cont_means_exp_exp1_converge.png", width = 8, height = 8)


#########################################
### res_exp_exp2
########################################
library(ggplot2)
library(ggpubr)
#library(rstatix)

setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

load("0421_1000s_res_exp_exp2_df_TPRPPV_converge.Rdata")
TPRPPV <- ggplot(data = as.data.frame(df_TPRPPV_converge), aes(x=lambdaK, y=value, fill=Methods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid( cols = vars(Criteria) , labeller = label_parsed) +
  xlab(expression(lambda^k)) +
  ylab("") +
  ylim(0.5,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="right")
TPRPPV
ggsave("0421_1000s_cont_boxplots_exp_exp2_converge.pdf", width = 6, height = 3)

load("0421_1000s_res_exp_exp2_mean_TPRPPV_converge.Rdata")
mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV_converge), aes(x=as.numeric(levels(lambdaK))[lambdaK], y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid( cols = vars(Criteria) , labeller = label_parsed) +
  xlab(expression(lambda^k)) +
  ylab("") +
  ylim(0.5,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")+
  scale_x_continuous(breaks=c(0.01,0.02,0.03))
mean_fig
ggsave("0421_1000s_cont_means_exp_exp2_converge.pdf", width = 8, height = 8)

#########################################
### res_exp_exp3
########################################
library(ggplot2)
library(ggpubr)
#library(rstatix)

setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

load("0422_1000s_res_exp_exp3_df_TPRPPV_converge.Rdata")
TPRPPV <- ggplot(data = as.data.frame(df_TPRPPV_converge), aes(x=nA, y=value, fill=Methods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid( cols = vars(Criteria) , labeller = label_parsed) +
  xlab(expression(n[A])) +
  ylab("") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="right")
TPRPPV
ggsave("0422_1000s_cont_boxplots_exp_exp3_converge.pdf", width = 6, height = 3)


load("0422_1000s_res_exp_exp3_mean_TPRPPV_converge.Rdata")
mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV_converge), aes(x=as.numeric(levels(nA))[nA], y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid( cols = vars(Criteria) , labeller = label_parsed) +
  xlab(expression(n[A])) +
  ylab("") +
  ylim(0.5,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")+
  scale_x_continuous(breaks=c(400,800,1200))
mean_fig
ggsave("0422_1000s_cont_means_exp_exp3_converge.pdf", width = 8, height = 8)

#########################################
### res_exp_exp4
########################################
library(ggplot2)
library(ggpubr)
#library(rstatix)

#setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

load("0420_1000s_res_exp_exp4_df_TPRPPV_converge.Rdata")
TPRPPV <- ggplot(data = as.data.frame(df_TPRPPV_converge), aes(x=K, y=value, fill=Methods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid( cols = vars(Criteria) , labeller = label_parsed) +
  xlab("K") +
  ylab("") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="right")
TPRPPV
ggsave("0420_1000s_cont_boxplots_exp_exp4_converge.pdf", width = 6, height = 3)

load("0420_1000s_res_exp_exp4_mean_TPRPPV_converge.Rdata")
mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV_converge), aes(x=as.numeric(levels(K))[K], y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid( cols = vars(Criteria) , labeller = label_parsed) +
  xlab("K") +
  ylab("") +
  ylim(0.5,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")+
  scale_x_continuous(breaks=c(2,3,4))
mean_fig
ggsave("0420_1000s_cont_means_exp_exp4_converge.pdf", width = 8, height = 8)
########################################################
### Robustness Uniform-exp
#######################################################
load("0410_1000s_res_unif_exp_df_TPRPPV_converge.RData")
TPRPPV <- ggplot(data = as.data.frame(df_TPRPPV_converge), aes(x=error, y=value, fill=Methods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid( prettylambdaE ~  Criteria , labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")
TPRPPV
ggsave("0410_1000s_cont_boxplots_unif_exp.pdf", width = 8, height = 8)

load("0410_1000s_res_unif_exp_mean_TPRPPV_converge.RData")
mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV_converge), aes(x=as.numeric(levels(error))[error], y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid( prettylambdaE ~  Criteria, labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0.5,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")+
  scale_x_continuous(breaks=c(0.1,0.2,0.3)) 
mean_fig
ggsave("0410_1000s_cont_means_unif_exp.pdf", width = 8, height = 8)

########################################################
### Robustness Exp-norm
#######################################################
load("0329_1000s_res_exp_norm_df_TPRPPV_converge.RData")
TPRPPV <- ggplot(data = as.data.frame(df_TPRPPV_converge), aes(x=error, y=value, fill=Methods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid( prettymeanE ~  Criteria , labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0.6,1)+
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")
TPRPPV
ggsave("0329_1000s_cont_boxplots_exp_norm.pdf", width = 8, height = 8)

load("0329_1000s_res_exp_norm_mean_TPRPPV_converge.RData")
mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV_converge), aes(x=as.numeric(levels(error))[error], y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid( prettymeanE ~  Criteria, labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0.5,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")+
  scale_x_continuous(breaks=c(0.1,0.2,0.3)) 
mean_fig
ggsave("0329_1000s_cont_means_exp_norm.pdf", width = 8, height = 8)
