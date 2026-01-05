

library(tidyverse)
library(reshape2)
library(ggpubr)


source('~/Benchmark-realdata-main/fig1/results_analysis/readTNBC_realbulk.R')
prd_list_TNBC_realbulk = prd_list



source('~/Benchmark-realdata-main/fig1/results_analysis/readER_realbulk.R')
prd_list_ER_realbulk = prd_list



key_TNBC <- readRDS("~/Benchmark-realdata-main/fig1/data_preprocess/key_TNBC.rds")


key_TNBC = key_TNBC[sort(row.names(key_TNBC)),sort(colnames(key_TNBC))]



key_ER <- readRDS("~/Benchmark-realdata-main/fig1/data_preprocess/key_ER.rds")

key_ER = key_ER[sort(row.names(key_ER)),sort(colnames(key_ER))]


key_all = cbind(key_TNBC,key_ER)


cor_vec_realbulk_list = list()

for(ii in names(prd_list_TNBC_realbulk)){
  prd_in_TNBC = prd_list_TNBC_realbulk[[ii]]
  prd_in_ER = prd_list_ER_realbulk[[ii]]
  
  prd_in_TNBC[setdiff(row.names(key_TNBC),row.names(prd_in_TNBC)),] = 0
  prd_in_ER[setdiff(row.names(key_ER),row.names(prd_in_ER)),] = 0
  
  prd_in_TNBC = prd_in_TNBC[row.names(key_all),sort(colnames(prd_in_TNBC))]
  prd_in_ER = prd_in_ER[row.names(key_all),sort(colnames(prd_in_ER))]
  prd_in = cbind(prd_in_TNBC,prd_in_ER)
  
  
  prd_in = prd_in[row.names(key_all),colnames(key_all)]
  

  name11 =row.names(prd_in)
  


  cor_in =c()
  for(kk in name11){
    cor_in = c(cor_in,cor(as.numeric(prd_in[kk,]),
                          as.numeric(key_all[kk,]),method = 'pearson'))

  }
  names(cor_in) = name11
  
  
  
  
  cor_in[is.na(cor_in)] =0
  cor_in =as.data.frame(cor_in)
  cor_in[,2] = ii
  cor_in[,3] = 'non-DP'
  cor_in[row.names(cor_in) %in% c('MonocyteMacrophage','Cycling_Myeloid','Endothelial_CXCL12','DCs'),3] ='DP'

  
  colnames(cor_in) = c('cor_vec','method_vec','category')
  cor_vec_realbulk_list[[ii]] = cor_in
  
}

names(cor_vec_realbulk_list) = names(prd_list_TNBC_realbulk)


ranked_list <- lapply(cor_vec_realbulk_list, function(sub_df) {
  sub_df$cor_rank <- rank(-sub_df$cor_vec, ties.method = "min")
  return(sub_df)
})

# 
combined_df <- do.call(rbind, ranked_list)

# 

combined_df =subset(combined_df,method_vec != 'bisque')

p_rank = ggplot(combined_df,aes(x = method_vec,y = cor_rank,fill = category)) +
  geom_boxplot(aes(color = category),
               width = .6, size = .7, alpha = .5 ) +
  theme_classic() +
  scale_color_manual(values= c('#6793ea','#ca7491')) +
  scale_fill_manual(values= c('#6793ea','#ca7491')) +
  # geom_hline(yintercept = 0.6, color = "red", linetype = "dashed", size = 1) +  # 添加辅助线
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
  stat_compare_means(
    method = "wilcox.test",
      method.args = list(alternative = "less")   
  )

p_rank





ggsave(
  filename = "~/Benchmark-realdata-main/fig1/results_analysis/fig/pdf/F1E.pdf", 
  plot = p_rank,  
  width = 4.5,  
  height = 4, 
  units = "in",  
  dpi = 300,  
  device = "pdf"  
)
