library(tidyverse)
library(reshape2)
setwd('~/ReCIDE/benchmark_syq/返修/github/')
source("./fig4_benchmark/code/function_2_bulk_total.R")

crit ='conservative'
if(crit == 'conservative'){
  source('./fig4_benchmark/code/function_conservative.R')
  
}else if(crit == 'permissive'){
  source('./fig4_benchmark/code/function_permissive.R')
}

disease_vec = c('BRCA','COAD','ESCA','HNSC','NSCLC','KIRC','STAD','PRAD')


F1_list = list()
for (disease in disease_vec) {

    key_trend = readRDS(paste0('./key_trend_in_fig3fig4/',disease,'.rds'))
    key_trend$PValue =  key_trend$PValue/2
    
    decon_results_TCGA =  readRDS(paste0('./benchmark_deconvolution_results/','TCGA-',disease,'.rds')) 
    decon_results_GEO =  readRDS(paste0('./benchmark_deconvolution_results/','GEO-',disease,'.rds')) 
    
    TCGA_pvalue = get_pvalue(decon_results_TCGA,
                                    key_trend)
    
    GEO_pvalue = get_pvalue(decon_results_GEO,
                                    key_trend)
    
    F1_list[[disease]] = compare_bulk2bulk(TCGA_pvalue,
                                          GEO_pvalue,
                                          disease,
                                          key_trend)

}







final <- do.call(rbind,F1_list)
final[,'method'] <- as.character(sapply(row.names(final),function(x){x=str_split(x,'\\.')[[1]][2]}))
final[is.na(final)] <- 0

final_melt <- melt(final)

miss_data_Jaccard = final_melt[final_melt$variable == 'Jaccard',]
print(miss_data_Jaccard %>% group_by(method) %>% summarise_each(mean))


miss_data_F1 = final_melt[final_melt$variable == 'F1',]

miss_data_F1 = miss_data_F1[miss_data_F1$disease !='PRAD',]
print(miss_data_F1 %>% group_by(method) %>% summarise_each(mean))






for(metric in c('Jaccard','F1')){
  if(metric == 'Jaccard'){
miss_data <- final_melt[final_melt$variable=='Jaccard',]
  }else{
    miss_data <- final_melt[final_melt$variable=='F1',]
    miss_data = miss_data[miss_data$disease !='PRAD',]
    
    
  }
miss_data$value[is.na(miss_data$value)] <- 0
miss_data$value <- as.numeric(miss_data$value)
mean_values <- data.frame()
for (method in unique(miss_data$method)){mean_values[method,'mean_value'] <- mean(miss_data[miss_data$method==method,'value'])}
mean_values["music",'color'] <- '#aad8f6'
mean_values["DWLS",'color'] <- "#a08ac1"
mean_values["CIBERSORT",'color'] <- '#ffa500'
# mean_values["bisque",'color'] <- '#F58840'
mean_values["bayes",'color'] <- '#b5d98f'#e79bbd
mean_values["ReCIDE",'color'] <- '#e79bbd'

miss_data[,'method']=factor(miss_data[,'method'],levels=row.names(mean_values[order(mean_values$mean_value),]))
# miss_data[,'method']=factor(miss_data[,'method'],levels=c("music","CIBERSORT","DWLS","bisque","bayes","ReCIDE"))


miss_data[,'rank'] <- ave(-miss_data$value, miss_data$disease, FUN = function(x) rank(x, ties.method = "min"))
miss_data[miss_data[,'rank']==1,'significant'] <- TRUE
miss_data[miss_data[,'rank']!=1,'significant'] <- FALSE

vec_color <-mean_values[order(mean_values$mean_value),'color']
p1 = ggplot(miss_data, aes(x = method, y = value, fill = method)) +
  geom_boxplot(aes(color = method),
               width = .6, size = .7, alpha = .5) +
  theme_classic() +
  scale_color_manual(values = vec_color) +
  scale_fill_manual(values = vec_color) +
  theme(
    axis.text.x = element_text(size = 10, face = "plain", angle = -45),
    axis.text.y = element_text(size = 10, face = "plain"),
    axis.title = element_text(size = 8, face = "plain"),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "plain", hjust = 0.5),
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.4, "cm"),
    legend.position = "bottom"
  ) +coord_flip() +
  scale_y_continuous(limits = c(0, 1))+  stat_summary(fun = mean, geom = "text", 
                                                      aes(label = round(..y.., 2)), vjust = -0.5, size = 3, color = "black")

# scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25))
p1

miss_data %>% group_by(method) %>% summarise_each(mean)

## 7 5 4 5
# 转换因子顺序以确保热图的顺序
miss_data$method <- factor(miss_data$method, levels = unique(miss_data$method))
miss_data$source <- ifelse(grepl("TCGA", miss_data$disease), "TCGA", "GEO")
#miss_data$disease <- factor(miss_data$disease, levels = unique(miss_data$disease))

miss_data[,'method']=factor(miss_data[,'method'],levels=row.names(mean_values[order(mean_values$mean_value),]))##level和boxplot保持一致
miss_data$source <- factor(miss_data$source, levels = c("TCGA","GEO"))


# 绘制热图 low='#54a1cd', high="#cc0000", midpoint=0.5, mid="#ffff53"
p2 = ggplot(miss_data, aes(x = disease, y = method, fill = value)) +
  geom_tile(color = NA) +  # 去除内部单元格边框
  # scale_fill_gradientn(colors =c('#68c8ff',"#ffff53","#cc0000"),limits = c(0, 1)) +
  scale_fill_gradientn(colors =c('white',"#fed486","#cc0000"), values = c(0, 0.5, 1), limits = c(0, 1)) +
  # scale_fill_gradientn(colors =c('white',"#f4c8b8","#c72120")) +
  # labs(title = "Heatmap of Algorithm Values by Dataset", fill = "Value") +
  theme_minimal() + 
  # geom_text(aes(label = ifelse(significant, "*", "")), color = "black", size = 3) + 
  #coord_fixed()+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()  # 去除背景网格线
  ) +
  coord_fixed(ratio = 0.8)+
  # 添加外部大框
  annotate("rect", xmin = 0.5, xmax = length(unique(miss_data$disease)) + 0.5,
           ymin = 0.5, ymax = length(unique(miss_data$method)) + 0.5,
           color = "black", fill = NA, size = 0.7)

p1+p2



ggsave(
  filename = paste0('./fig4_benchmark/results_and_fig/',crit,'_',metric,'_boxplot.pdf'),  # 保存的文件名
  plot = p1,  # 要保存的图形对象
  width = 3,  # 宽度（默认单位为英寸）
  height = 4,  # 高度
  units = "in",  # 单位：英寸
  dpi = 300,  # 分辨率（每英寸点数）
  device = "pdf"  # 设备类型
)


ggsave(
  filename = paste0('./fig4_benchmark/results_and_fig/',crit,'_',metric,'_heatmap.pdf'),  # 保存的文件名
  plot = p2,  # 要保存的图形对象
  width = 8,  # 宽度（默认单位为英寸）
  height = 4,  # 高度
  units = "in",  # 单位：英寸
  dpi = 300,  # 分辨率（每英寸点数）
  device = "pdf"  # 设备类型
)
}