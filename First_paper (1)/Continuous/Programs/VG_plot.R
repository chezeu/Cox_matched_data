library(ggplot2)
library(ggpubr)

setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Results")

load("0406_1000s_res_exp_exp1_mean_TPRPPV_converge.Rdata")
mean_fig <- ggplot(data = as.data.frame(mean_TPRPPV_converge), aes(x=error, y=Mean, color=Methods, group = Methods)) +
  geom_point() +
  geom_line() +
  facet_grid( prettylambdaE ~  Criteria , labeller = label_parsed) +
  xlab(expression(e^k)) +
  ylab("") +
  ylim(0.5,1)+
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_bw() +
  theme(legend.position="top")
mean_fig
