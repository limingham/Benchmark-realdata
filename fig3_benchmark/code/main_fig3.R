library(tidyverse)
library(stringr)
library(reshape2)
library(readr)
library(pheatmap)


setwd('~/ReCIDE/benchmark_syq/返修/github/')
source('./fig3_benchmark/code/function_use.R')


#### conservative permissive
crit = 'conservative'


disease_vec = c('BRCA','COAD','ESCA','HNSC','NSCLC','KIRC','STAD')


F1_list = list()
for (disease in disease_vec) {
  for (database in c('TCGA','GEO')) {
    key_trend = readRDS(paste0('./key_trend_in_fig3fig4/',disease,'.rds'))
    key_trend$PValue =  key_trend$PValue/2
    
    decon_results =  readRDS(paste0('./benchmark_deconvolution_results/',database,'-',disease,'.rds')) 
    
    F1_list[[paste0(database,'-',disease)]] = function_cal_F1(decon_results,
                                                              key_trend,
                                                              paste0(database,'-',disease),
                                                              criteria = crit)
  }
}




final_df <- do.call(rbind,F1_list)
final_df[,'method'] <- as.character(sapply(row.names(final_df),function(x){x=str_split(x,'\\.')[[1]][2]}))

melt_final_df <- melt(final_df)
melt_final_df[is.na(melt_final_df)] <- 0
melt_final_df = melt_final_df[melt_final_df$value> (-5),]
# 





melt_final_df = melt_final_df[melt_final_df$variable == 'F1',]

print(melt_final_df %>% group_by(method) %>% summarise_each(mean))
print(melt_final_df %>% group_by(disease) %>% summarise_each(mean))



miss_data = melt_final_df


mean_values <- data.frame()
for (method in unique(miss_data$method)){mean_values[method,'mean_value'] <- mean(miss_data[miss_data$method==method,'value'])}
mean_values["music",'color'] <- '#aad8f6'
mean_values["DWLS",'color'] <- "#a08ac1"
mean_values["CIBERSORT",'color'] <- '#ffa500'
# mean_values["bisque",'color'] <- '#F58840'
mean_values["bayes",'color'] <- '#b5d98f'
mean_values["ReCIDE",'color'] <- '#e79bbd'

miss_data[,'method']=factor(miss_data[,'method'],levels=row.names(mean_values[order(mean_values$mean_value),]))
miss_data[,'rank'] <- ave(-miss_data$value, miss_data$disease, FUN = function(x) rank(x, ties.method = "min"))
miss_data[miss_data[,'rank']==1,'significant'] <- TRUE
miss_data[miss_data[,'rank']!=1,'significant'] <- FALSE

vec_color <-mean_values[order(mean_values$mean_value),'color']
p_box = ggplot(miss_data,aes(x = method,y = value,fill = method)) +
  geom_boxplot(aes(color = method),
               width = .6, size = .7, alpha = .5 ) +
  theme_classic() +
  scale_color_manual(values= vec_color) +
  scale_fill_manual(values= vec_color) +
  theme(
    axis.text.x = element_text(size = 10, face = "plain", angle = -45),
    axis.text.y = element_text(size = 10, face = "plain"),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.title = element_text(size = 8, face = "plain"),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "plain", hjust = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.4, "cm"),
    # axis.title = element_text(size = 8)
    legend.position = "bottom"
  )+
  scale_y_continuous(limits = c(0, 1)) +
  coord_flip()

p_box
# 转换因子顺序以确保热图的顺序
miss_data$source <- ifelse(grepl("TCGA", miss_data$disease), "TCGA", "GEO")
# miss_data[is.na(miss_data)] =0


miss_data[,'method']=factor(miss_data[,'method'],levels=row.names(mean_values[order(mean_values$mean_value),]))##level和boxplot保持一致
miss_data$source <- factor(miss_data$source, levels = c("TCGA","GEO"))
miss_data <- miss_data %>%
  arrange(source)
miss_data$disease <- factor(miss_data$disease, levels = sort(unique(miss_data$disease)))

p1 <- ggplot(miss_data, aes(x = disease, y = method, fill = value)) +
  geom_tile(color = NA) +  # 去除内部单元格边框
  scale_fill_gradientn(colors =c('white',"#fed486","#cc0000"), values = c(0, 0.5, 1), limits = c(0, 1)) +
  # labs(title = "Heatmap of Algorithm Values by Dataset", fill = "Value") +
  theme_minimal() +
  #geom_text(aes(label = ifelse(significant, "*", "")), color = "black", size = 3) +
  #coord_fixed()+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()  # 去除背景网格线
  ) +
  coord_fixed(ratio = 1.3)+
  # 添加外部大框
  annotate("rect", xmin = 0.5, xmax = length(unique(miss_data$disease)) + 0.5,
           ymin = 0.5, ymax = length(unique(miss_data$method)) + 0.5,
           color = "black", fill = NA, size = 0.7)

p1


if(crit == 'permissive'){
  ggsave(
    filename = paste0("./fig3_benchmark/figure/permissive_boxplot.pdf"),  # 保存的文件名
    plot = p_box,  # 要保存的图形对象
    width = 3,  # 宽度（默认单位为英寸）
    height = 4,  # 高度
    units = "in",  # 单位：英寸
    dpi = 300,  # 分辨率（每英寸点数）
    device = "pdf"  # 设备类型
  )
  
  
  ggsave(
    filename = "./fig3_benchmark/figure/permissive_heatmap.pdf",  # 保存的文件名
    plot = p1,  # 要保存的图形对象
    width = 8,  # 宽度（默认单位为英寸）
    height = 4,  # 高度
    units = "in",  # 单位：英寸
    dpi = 300,  # 分辨率（每英寸点数）
    device = "pdf"  # 设备类型
  )
  
  saveRDS(melt_final_df,file = "./fig3_benchmark/results/permissive.rds")
}

if(crit == 'conservative'){
  ggsave(
    filename = paste0("./fig3_benchmark/figure/conservative_boxplot.pdf"),  # 保存的文件名
    plot = p_box,  # 要保存的图形对象
    width = 3,  # 宽度（默认单位为英寸）
    height = 4,  # 高度
    units = "in",  # 单位：英寸
    dpi = 300,  # 分辨率（每英寸点数）
    device = "pdf"  # 设备类型
  )
  
  
  ggsave(
    filename = "./fig3_benchmark/figure/conservative_heatmap.pdf",  # 保存的文件名
    plot = p1,  # 要保存的图形对象
    width = 8,  # 宽度（默认单位为英寸）
    height = 4,  # 高度
    units = "in",  # 单位：英寸
    dpi = 300,  # 分辨率（每英寸点数）
    device = "pdf"  # 设备类型
  )
  
  saveRDS(melt_final_df,file = "./fig3_benchmark/results/conservative.rds")
}