
source("~/Benchmark-realdata-main/fig2/results_analysis/fig/F2A/preprocess.R")

rank_sample_average = apply(pcc_real_BRCA_df,2,mean)


source("~/Benchmark-realdata-main/fig2/results_analysis/fig/celltype-wise-PCC/all_celltypes.R")
cor_vec_realbulk_df = cor_vec_realbulk_df[,names(rank_sample_average)]


vec_sample_2celltypes = as.numeric(cor(rank_sample_average,t(cor_vec_realbulk_df),method = 'spearman'))
mean(vec_sample_2celltypes)
median(vec_sample_2celltypes)
quantile(vec_sample_2celltypes)






c1_df = as.data.frame(vec_sample_2celltypes)
c2_df = as.data.frame(vec_celltypes_2celltypes)

c1_df[,2] = 'vec_sample_2celltypes'
c2_df[,2] = 'vec_celltypes_2celltypes'
# c3_df[,2] = 'vec_DP_2DP'

colnames(c1_df) = c('MRC','category')
colnames(c2_df) = c('MRC','category')
# colnames(c3_df) = c('MRC','category')

plot_df = rbind(c1_df,c2_df)


plot_df$category <- factor(
  plot_df$category,
  levels = c("vec_sample_2celltypes", "vec_celltypes_2celltypes")
)


p = ggplot(plot_df,aes(x = category,y = MRC,fill = category)) +
  geom_boxplot(aes(color = category),
               width = .6, size = .7, alpha = .5 ) +
  theme_classic() +
  scale_color_manual(values= c('#eb746a','green4')) +
  scale_fill_manual(values= c('#eb746a','green4')) +
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
  )

p






ggsave(
  filename = "~/Benchmark-realdata-main/fig2/results_analysis/fig/pdf/F2C.pdf",  # 保存的文件名
  plot = p,  # 要保存的图形对象
  width = 3,  # 宽度（默认单位为英寸）
  height = 5,  # 高度
  units = "in",  # 单位：英寸
  dpi = 300,  # 分辨率（每英寸点数）
  device = "pdf"  # 设备类型
)



