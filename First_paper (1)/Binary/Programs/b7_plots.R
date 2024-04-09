library(ggplot2)
library(ggpubr)



#setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Binary\\Results")
setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Binary\\Results")
###############################################################
load("0308_1000s_res1_df_TPRPPV.Rdata")
TPRPPV <- ggplot(data = as.data.frame(df_TPRPPV), aes(x=error, y=value, fill=Methods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid( prettyK ~  Criteria , labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")
TPRPPV
ggsave("0308_1000s_binary_res1_boxplots.pdf", width = 8, height = 8)


load("0308_1000s_res1_mean_TPRPPV.Rdata")
mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV), aes(x=as.numeric(levels(error))[error], y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid( prettyK ~  Criteria , labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")+
  scale_x_continuous(breaks=c(0.02,0.04,0.06)) 
mean_fig
ggsave("0308_1000s_binary_res1_means.pdf", width = 8, height = 8)

load("0308_1000s_res1_df_TPRPPV_converge.Rdata")
TPRPPV <- ggplot(data = as.data.frame(df_TPRPPV_converge), aes(x=error, y=value, fill=Methods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid( prettyK ~  Criteria , labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")
TPRPPV

ggsave("0308_1000s_binary_res1_boxplots_converge.pdf", width = 8, height = 8)


load("0308_1000s_res1_mean_TPRPPV_converge.Rdata")
mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV_converge), aes(x=as.numeric(levels(error))[error], y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid( prettyK ~  Criteria , labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0.4,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top", legend.text=element_text(size=12), strip.text.x = element_text(size = 12))+
  scale_x_continuous(breaks=c(0.02,0.04,0.06)) 
mean_fig

ggsave("0308_1000s_binary_res1_means_converge.pdf", width = 8, height = 8)
ggsave("0308_1000s_binary_res1_means_converge.png", width = 8, height = 8)
###################### Res2
load("0324_1000s_res2_df_TPRPPV_converge.Rdata")

TPRPPV <- ggplot(data = as.data.frame(df_TPRPPV_converge), aes(x=prev, y=value, fill=Methods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid(  cols = vars(Criteria) , labeller = label_parsed) +
  xlab(expression(p^k)) +
  ylab("") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="right")
TPRPPV

ggsave("0324_1000s_binary_res2_boxplots_converge.pdf", width = 6, height = 3)

load("0324_1000s_res2_mean_TPRPPV_converge.Rdata")

mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV_converge), aes(x=as.numeric(levels(prev))[prev], y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid(  cols = vars(Criteria) , labeller = label_parsed) +
  xlab(expression(p^k)) +
  ylab("") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")+
  scale_x_continuous(breaks=c(0.1,0.2,0.3)) 
mean_fig

ggsave("0324_1000s_binary_res2_means_converge.pdf", width = 8, height = 8)

###################### Res3
load("0324_1000s_res3_df_TPRPPV_converge.Rdata")

TPRPPV <- ggplot(data = as.data.frame(df_TPRPPV_converge), aes(x=nA, y=value, fill=Methods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid(  cols = vars(Criteria) , labeller = label_parsed) +
  xlab(expression(n[A])) +
  ylab("") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="right")
TPRPPV

ggsave("0324_1000s_binary_res3_boxplots_converge.pdf", width = 6, height = 3)

load("0324_1000s_res3_mean_TPRPPV_converge.Rdata")

mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV), aes(x=as.numeric(levels(nA))[nA], y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid(  cols = vars(Criteria) , labeller = label_parsed) +
  xlab(expression(n[A])) +
  ylab("") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")+
  scale_x_continuous(breaks=c(400,800,1200)) 
mean_fig

ggsave("0324_1000s_binary_res3_means_converge.pdf", width = 8, height = 8)
