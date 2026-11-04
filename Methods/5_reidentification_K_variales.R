
################ plot  N(K)
library(ggplot2)

S  <- 21 #max value of k
K  <- 14:S
nB <- 2000
nA <- 1000
N<- 1+ (nB - 1)*(2^(-K))

base <- data.frame(K = K, N = N)  # data.frame  

p_a <- ggplot(base, aes(K, N)) +
  geom_line(aes(linetype = "Degree of identification"), linewidth = 0.5) +
  geom_point(size = 1.2) +
  geom_hline(aes(yintercept = 1, linetype = "N(K) = 1"), linewidth = 0.5) +
  scale_linetype_manual(values = c("Degree of identification" = "solid",
                                   "N(K) = 1"                 = "dashed")) +
  scale_x_continuous(breaks = K) +              
  labs(x = "Number of matching variables K",
       y = "Degree of identification N(K)") +
  theme_bw(base_size = 8) +
  theme(panel.grid        = element_blank(),
        legend.title      = element_blank(),
        legend.position   = c(0.72, 0.85),
        legend.background = element_rect(colour = "black", linewidth = 0.3))

p_a

# export the figure in .pdf
figure=p_a

figure_revue <- figure +
  labs(title = paste0(" Fidelette_CHEZEU ", " - Figure ", 1, " - page 16")) +
  theme(plot.title = element_text(size = 16, hjust = 0.5)) 
fichier <- paste0("Fidelette_Chezeu", "_Figure", 1)
ggsave(paste0(fichier, ".pdf"), figure_revue, device = cairo_pdf,
       width = 20, height = 10, units = "cm")
browseURL(paste0(fichier, ".pdf"))
 
