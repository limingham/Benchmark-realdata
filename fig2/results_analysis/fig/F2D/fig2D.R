
source("~/Benchmark-realdata-main/fig2/results_analysis/fig/celltype-wise-PCC/DP_celltypes.R")


source("~/Benchmark-realdata-main/fig2/results_analysis/fig/celltype-wise-PCC/all_celltypes.R")




c2_df = as.data.frame(vec_noDP_2noDP)
c3_df = as.data.frame(vec_DP_2DP)

c2_df[,2] = 'vec_noDP_2noDP'
c3_df[,2] = 'vec_DP_2DP'

colnames(c2_df) = c('MRC','category')
colnames(c3_df) = c('MRC','category')

plot_df = rbind(c2_df,c3_df)


plot_df$category <- factor(
  plot_df$category,
  levels = c("vec_noDP_2noDP", "vec_DP_2DP")
)


p = ggplot(plot_df,aes(x = category,y = MRC,fill = category)) +
  geom_boxplot(aes(color = category),
               width = .6, size = .7, alpha = .5 ) +
  theme_classic() +
  scale_color_manual(values= c('#6793ea','#ca7491')) +
  scale_fill_manual(values= c('#6793ea','#ca7491')) +
  theme(
    axis.text.x = element_text(size = 14, face = "plain", angle = -45),
    axis.text.y = element_text(size = 14, face = "plain"),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.title = element_text(size = 14, face = "plain"),
    plot.title = element_text(size = 14, face = "plain", hjust = 0.5),
    plot.subtitle = element_text(size = 14, face = "plain", hjust = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    legend.text = element_text(size = 14),
    legend.spacing.x = unit(0.4, "cm"),
    # axis.title = element_text(size = 8)
    legend.position = "bottom"
  )+
  stat_compare_means(method.args = list(alternative = "greater"))

p



ggsave(
  filename = "~/Benchmark-realdata-main/fig2/results_analysis/fig/pdf/F2D.pdf",  
  plot = p,  
  width = 3,  
  height = 5,  
  units = "in", 
  dpi = 300,  
  device = "pdf" 
)



