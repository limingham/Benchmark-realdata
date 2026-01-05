library(tidyverse)
library(reshape2)
library(ggpubr)

source("~/Benchmark-realdata-main/fig1/results_analysis/fig/F1A/preprocess.R")


pcc_pseudo_BRCA_df$sample = row.names(pcc_pseudo_BRCA_df)
pcc_real_BRCA_df$sample = row.names(pcc_real_BRCA_df)


pcc_psedobulk <- melt(pcc_pseudo_BRCA_df, 
                      id.vars = "sample",        
                      variable.name = "Method",  
                      value.name = "Value")    
pcc_psedobulk$category = 'psedobulk'

pcc_realbulk <- melt(pcc_real_BRCA_df, 
                     id.vars = "sample",       
                     variable.name = "Method",  
                     value.name = "Value")    
pcc_realbulk$category = 'realbulk'


plot_df = rbind(pcc_psedobulk,pcc_realbulk)

pcc_psedobulk %>% group_by(Method) %>% summarise_each(mean)
pcc_realbulk %>% group_by(Method) %>% summarise_each(mean)

plot_df$Method = as.character(plot_df$Method)
p = ggplot(plot_df,aes(x = Method,y = Value,fill = category)) +
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
  stat_compare_means(method = 'wilcox.test',paired = TRUE,label = "p.signif")

p
ggsave(
  filename = "~/Benchmark-realdata-main/fig1/results_analysis/fig/pdf/F1A.pdf",  
  plot = p,  
  width = 5,  
  height = 4,  
  units = "in",  
  dpi = 300,  
  device = "pdf"  
)