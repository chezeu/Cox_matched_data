df3_TPR11 = subset(df, df$methods == "TPR-11 FS3" & df$K == 20 & df$error == 0.04)
df3_converge = subset(df, df$methods == "converge FS3" & df$K == 20 & df$error == 0.04)
df3_TPR11$converge = df3_converge$value

df3_TPR11_converge = df3_TPR11[df3_TPR11$converge==1,]

df3_TPR11_nonconverge = df3_TPR11[df3_TPR11$converge==0,]

df3_PPV11 = subset(df, df$methods == "PPV-11 FS3" & df$K == 20 & df$error == 0.04)
df3_converge = subset(df, df$methods == "converge FS3" & df$K == 20 & df$error == 0.04)
df3_PPV11$converge = df3_converge$value

df3_PPV11_converge = df3_PPV11[df3_PPV11$converge==1,]

df3_PPV11_nonconverge = df3_PPV11[df3_PPV11$converge==0,]

library(ggplot2)
library(ggpubr)
library(rstatix)
setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Binary\\Results")

TPR <- ggplot(data = as.data.frame(df3_TPR11_converge), aes(x=prettyerror, y=value, fill=prettymethods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid( prettyK ~  prettyrules , labeller = label_parsed) +
  xlab("") +
  ylab("TPR") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")

PPV <- ggplot(data = as.data.frame(df3_PPV11_converge), aes(x=prettyerror, y=value, fill=prettymethods)) +
  geom_boxplot(outlier.size=0) +
  # stat_summary(geom = "crossbar", width=0.3, fatten=0, color="white", position = position_dodge(width = 0.75), fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  facet_grid( prettyK ~  prettyrules , labeller = label_parsed) +
  xlab("e") +
  ylab("PPV") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() + #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")

ggarrange(TPR, PPV, ncol = 2, nrow = 1, common.legend = TRUE)
ggsave("0802_FS3_converge.pdf", width = 10, height = 8, unit = "cm")

TPR <- ggplot(data = as.data.frame(df3_TPR11_nonconverge), aes(x=prettyerror, y=value, fill=prettymethods)) +
  geom_boxplot(outlier.size=0) +
  facet_grid( prettyK ~  prettyrules , labeller = label_parsed) +
  xlab("") +
  ylab("TPR") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() +# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")

PPV <- ggplot(data = as.data.frame(df3_PPV11_nonconverge), aes(x=prettyerror, y=value, fill=prettymethods)) +
  geom_boxplot(outlier.size=0) +
  # stat_summary(geom = "crossbar", width=0.3, fatten=0, color="white", position = position_dodge(width = 0.75), fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  facet_grid( prettyK ~  prettyrules , labeller = label_parsed) +
  xlab("e") +
  ylab("PPV") +
  ylim(0,1)+
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_fill_discrete("Methods") +
  theme_bw() + #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top")

ggarrange(TPR, PPV, ncol = 2, nrow = 1, common.legend = TRUE)
ggsave("0802_FS3_nonconverge.pdf", width = 10, height = 8, unit = "cm")