library(survival)
library(survminer)
library(reshape2)
library(stringr)
library(tidyverse)
library(readr)
library(pheatmap)

setwd('~/Benchmark-realdata-main/')
source("./fig5_benchmark/code/function_used.R")
list_clinical = readRDS("./clinical_in_fig5/clinical_in_fig5.rds")

disease_vec = c('COAD','HNSC','NSCLC','PAAD')


prd_list = list()
for (disease in disease_vec) {
  for (database in c('TCGA-','GEO-')) {
    
 
  key_trend = readRDS(paste0('./key_trend_in_fig5/',disease,'.rds'))
  
  decon_results =  readRDS(paste0('./benchmark_deconvolution_results/',database,disease,'.rds')) 
  
  decon_results_use = read_deconvolurtion_results(decon_results,
                           key_trend,
                           paste0(database,disease))
  

  # 
  prd_list[[paste0(database,disease)]] <- decon_results_use  # 追加单个元素
  
  }
}

names(prd_list) = c('TCGA-COAD','GEO-COAD','TCGA-HNSC_neg','GEO-HNSC_neg','TCGA-','GEO-','TCGA-PAAD','GEO-PAAD')



list_prd = list()

for (dis in names(prd_list)) {
  
  if(dis != 'TCGA-' & dis != 'GEO-'){
    list_in = list()
    for (dis_in in 1:length(prd_list[[dis]])) {
      list_in[[dis_in]] = prd_list[[dis]][[dis_in]][[2]]
    }
    names(list_in) = names(prd_list[[dis]])
    
    list_prd[[dis]] = list_in
  }else{
    list_in = list()
    for (dis_in in 1:length(prd_list[[dis]])) {
      list_in[[dis_in]] = prd_list[[dis]][[dis_in]][[1]]
    }
    names(list_in) = names(prd_list[[dis]])
    
    list_prd[[paste0(dis,'LUSC')]] = list_in
    
    list_in = list()
    for (dis_in in 1:length(prd_list[[dis]])) {
      list_in[[dis_in]] = prd_list[[dis]][[dis_in]][[2]]
    }
    names(list_in) = names(prd_list[[dis]])
    
    list_prd[[paste0(dis,'LUAD')]] = list_in
    
  }}

for (mm in 1:length(list_prd)) {
  list_prd[[mm]][['bayes']] = list_prd[[mm]][['bayes']][sort(row.names(list_prd[[mm]][['bayes']])),
                                                        sort(colnames(list_prd[[mm]][['bayes']]))]
  
  list_prd[[mm]][['CIBERSORT']] = list_prd[[mm]][['CIBERSORT']][sort(row.names(list_prd[[mm]][['CIBERSORT']])),
                                                        sort(colnames(list_prd[[mm]][['CIBERSORT']]))]
  list_prd[[mm]][['DWLS']] = list_prd[[mm]][['DWLS']][sort(row.names(list_prd[[mm]][['DWLS']])),
                                                                sort(colnames(list_prd[[mm]][['DWLS']]))]
  list_prd[[mm]][['music']] = list_prd[[mm]][['music']][sort(row.names(list_prd[[mm]][['music']])),
                                                      sort(colnames(list_prd[[mm]][['music']]))]
  list_prd[[mm]][['ReCIDE']] = list_prd[[mm]][['ReCIDE']][sort(row.names(list_prd[[mm]][['ReCIDE']])),
                                                          sort(colnames(list_prd[[mm]][['ReCIDE']]))]
  
  colnames(list_prd[[mm]][['CIBERSORT']]) = colnames(list_prd[[mm]][['bayes']])
  colnames(list_prd[[mm]][['DWLS']]) = colnames(list_prd[[mm]][['bayes']])
  colnames(list_prd[[mm]][['music']]) = colnames(list_prd[[mm]][['bayes']])
  colnames(list_prd[[mm]][['ReCIDE']]) = colnames(list_prd[[mm]][['bayes']])
  
  for (jj in 1:length(list_prd[[mm]])) {
    colnames(list_prd[[mm]][[jj]]) = substr(colnames(list_prd[[mm]][[jj]]),1,12)
    
  }
}







list_p = list()
for (dise in names(list_prd)) {
  list_dise = list()
  for (ii in names(list_prd[[dise]])) {
    
    list_in = list()
    
    # cell_table = lasso[[dise]][[ii]] %>% table() %>% as.data.frame()
    # cell_table =cell_table[cell_table[,2]>10,]
    
    cell_table = list_prd[[dise]][[ii]]
    cell_table[,1] =row.names(cell_table)
    
    
    if(nrow(list_prd[[dise]][[ii]])>0){
      for(cc in cell_table[,1]){
        
        prd_in = list_prd[[dise]][[ii]]
        
        # prd_in = prd_in[cell_table[,1],]
        
        clinical_in = list_clinical[[dise]]
        
        
        if(nrow(prd_in)>0){
          rn =intersect(colnames(prd_in),row.names(clinical_in))
          prd_in = prd_in[,rn]
          clinical_in = clinical_in[rn,]
          
          
          gene_exp2<-as.data.frame(t(prd_in[cc,]))
          gene_exp2[is.na(gene_exp2)] = 0
          
          gene_exp2 = gene_exp2[gene_exp2[,1]>0.001,,drop =FALSE]
          clinical_in = clinical_in[row.names(gene_exp2),]
          # row.names(gene_exp2) = patient_inter
          if(substr(dise,1,4) == 'TCGA'){
            list_in[[cc]] = func_calculate_hr_TCGA(gene_exp2,clinical_in)
          }else if (substr(dise,1,4) == 'GEO-'){
            list_in[[cc]] = func_calculate_hr_GEO(gene_exp2,clinical_in,dise)
          }
          
        }
      }
      names(list_in) = cell_table[,1]
      # names(list_in) = row.names(list_prd[[dise]][[ii]])
      
    }
    list_dise[[ii]] = list_in
    
    # names(list_in) = row.names(list_prd[[dise]][[ii]])
    # list_dise[[ii]] = list_in
  }
  names(list_dise) =  names(list_prd[[dise]])
  
  list_p[[dise]] = list_dise
}
names(list_p) = names(list_prd)

dis2 =  c('COAD','HNSC_neg','LUSC','LUAD','PAAD')


cellname = c()


results_list = list()
# p_com_list = list()
# n_i = dis2[1]
# mm = names(list_p_in_TCGA)[1]
for (n_i in dis2) {
  results_list_in = list()
  
  list_p_in_TCGA = list_p[[paste0('TCGA-',n_i)]]
  list_p_in_GEO = list_p[[paste0('GEO-',n_i)]]
  
  for (mm in names(list_p_in_TCGA)) {
    
    
    df1 = as.data.frame(list_p_in_TCGA[[mm]])
    # colnames(df1) = names(list_p[[n_i]])
    df1 = df1[,sort(colnames(df1))]
    
    df2 = as.data.frame(list_p_in_GEO[[mm]])
    # colnames(df2) = names(GEO[[n_i]])
    df2 = df2[,sort(colnames(df2))]
    
    cn = intersect(colnames(df1),colnames(df2))
    
    
    if(length(cn)>0){
      df1 = df1[,cn,drop =FALSE]
      df2 = df2[,cn,drop =FALSE]  
      df = rbind(df1,df2)
      df['celltype',] = colnames(df)
      df['methods',] = mm
      
      row.names(df) = c('direction_TCGA','p_TCGA','prop_median_TCGA','direction_GEO','p_GEO','prop_median_GEO','celltype','methods')
      df['direction_TCGA',] <- ifelse(df['direction_TCGA',] > 1, "+", 
                                      ifelse(df['direction_TCGA',] < 1, "-", df['direction_TCGA',]))
      
      df['direction_GEO',] <- ifelse(df['direction_GEO',] > 1, "+", 
                                     ifelse(df['direction_GEO',] < 1, "-", df['direction_GEO',]))
      
      df['p_TCGA',] =p.adjust(df['p_TCGA',], method = "fdr")
      df['p_GEO',] =p.adjust(df['p_GEO',], method = "fdr")
      
      # df= df[c('p_TCGA','direction_TCGA','p_GEO','direction_GEO','celltype','methods'),]
      
      
      results_list_in[[mm]] = df
      
    }
    
    if(length(results_list_in) >0){
      result_df_in = do.call(cbind, results_list_in)
      results_list[[n_i]] =  as.data.frame(t(result_df_in))
      
    }
  }}








thr = 0.05



# names(results_list)
ratio_TT_list = list()
ratio_TT_list1 = list()
ratio_TT_list2 = list()
for (di in 1:length(results_list)) {
  TNBC_df <- results_list[[di]]
  # TNBC_df <- readRDS(file_dir1)
  TNBC_df[,'p_TCGA'] = as.numeric(TNBC_df[,'p_TCGA'])
  TNBC_df[,"p_GEO"] = as.numeric(TNBC_df[,"p_GEO"])
  TNBC_df[,'prop_median_TCGA'] = as.numeric(TNBC_df[,'prop_median_TCGA'])
  TNBC_df[,"prop_median_GEO"] = as.numeric(TNBC_df[,"prop_median_GEO"])
  

  
  TNBC_df[TNBC_df$prop_median_TCGA ==0,'direction_TCGA'] ='0'
  TNBC_df[TNBC_df$prop_median_GEO ==0,'direction_GEO'] ='0'
  TNBC_df = TNBC_df[TNBC_df$prop_median_TCGA>0 | TNBC_df$prop_median_GEO>0,]
  # 
  TNBC_df1 = TNBC_df
  
  
  
  ##筛选出单显著
  TNBC_df = TNBC_df[TNBC_df[,"p_TCGA"] < thr | TNBC_df[,"p_GEO"] < thr,]
  
  
  
  ###都显著且同向
  ratio_TT = c()
  ratio_TT1 = c()
  ratio_TT2 =c()
  for (method in unique(TNBC_df1[,'methods'])) {
    TNBC_df_in = TNBC_df[TNBC_df[,'methods'] %in% method,]
    if(nrow(TNBC_df_in)>0){
      
      df_inter = TNBC_df_in[((TNBC_df_in[,"p_TCGA"] < thr) | (TNBC_df_in[,"p_GEO"] < thr)) & 
                              (TNBC_df_in[,"direction_TCGA"] == TNBC_df_in[,"direction_GEO"])  ,]
      
      # & 
      
      ratio_TT1 = c(ratio_TT1,nrow(df_inter))
      ratio_TT2 = c(ratio_TT2,nrow(TNBC_df_in))
      
      ratio_TT = c(ratio_TT,nrow(df_inter)/nrow(TNBC_df_in))
    }else{ratio_TT = c(ratio_TT,NA)
    ratio_TT1 = c(ratio_TT1,NA)
    ratio_TT2 = c(ratio_TT2,NA)
    }
  }
  
  names(ratio_TT) = unique(TNBC_df1[,'methods'])
  ratio_TT_list[[di]] = as.data.frame(ratio_TT)
  
  names(ratio_TT1) = unique(TNBC_df1[,'methods'])
  ratio_TT_list1[[di]] = as.data.frame(ratio_TT1)
  names(ratio_TT2) = unique(TNBC_df1[,'methods'])
  ratio_TT_list2[[di]] = as.data.frame(ratio_TT2)
  
}







ratio_TT_df1 <- do.call(cbind, ratio_TT_list1)
colnames(ratio_TT_df1) = names(results_list)
ratio_TT_df1[is.na(ratio_TT_df1)] = 0



ratio_TT_df2 <- do.call(cbind, ratio_TT_list2)
colnames(ratio_TT_df2) = names(results_list)
ratio_TT_df2[is.na(ratio_TT_df2)] = 0


aa1 = apply(ratio_TT_df1, 1, sum)/apply(ratio_TT_df2, 1, sum)



ratio_TT_df <- do.call(cbind, ratio_TT_list)
colnames(ratio_TT_df) = names(results_list)
ratio_TT_df[is.na(ratio_TT_df)] = 0



ratio_TT_df$method = row.names(ratio_TT_df)
ratio_TT_df = ratio_TT_df[ratio_TT_df$method!='bisque',]

miss_data =melt(ratio_TT_df)
colnames(miss_data) = c('method','disease','value')

miss_data$value[is.na(miss_data$value)] <- 0
miss_data$value <- as.numeric(miss_data$value)



mean_values <- data.frame()
for (method in unique(miss_data$method)){mean_values[method,'mean_value'] <- mean(miss_data[miss_data$method==method,'value'])}
mean_values["music",'color'] <- '#aad8f6'
mean_values["DWLS",'color'] <- "#a08ac1"
mean_values["CIBERSORT",'color'] <- '#ffa500'
# mean_values["bisque",'color'] <- '#F58840'
mean_values["bayes",'color'] <- '#b5d98f'
mean_values["ReCIDE",'color'] <- '#e79bbd'

# miss_data[,'method']=factor(miss_data[,'method'],levels=row.names(mean_values[order(mean_values$mean_value),]))
miss_data[,'method']=factor(miss_data[,'method'],levels=rev(c('ReCIDE','bayes','CIBERSORT','DWLS','music')))##level和boxplot保持一致

miss_data[,'rank'] <- ave(-miss_data$value, miss_data$disease, FUN = function(x) rank(x, ties.method = "min"))


# saveRDS(miss_data,file = '~/ReCIDE/benchmark_syq/返修/汇总图/data/prognostic_loose.rds')

miss_data[miss_data[,'rank']==1,'significant'] <- TRUE
miss_data[miss_data[,'rank']!=1,'significant'] <- FALSE



##############################
#############################
#######################
miss_data_mean = miss_data %>% group_by(method) %>% summarise_each(mean) %>% as.data.frame() 


miss_data_mean[,'value'] = aa1[as.character(miss_data_mean$method)]
miss_data_mean[,'method']=factor(miss_data_mean[,'method'],levels=rev(c('ReCIDE','bayes','CIBERSORT','DWLS','music')))##level和boxplot保持一致

vec_color <-mean_values[rev(c('ReCIDE','bayes','CIBERSORT','DWLS','music')),'color']
p_box = ggplot(miss_data_mean,aes(x = method,y = value,fill = method)) +
  geom_col(aes(color = method), 
           width = 0.6, size = 0.7,  position = "dodge") +
  theme_classic() +
  scale_color_manual(values= vec_color) +
  scale_fill_manual(values= vec_color) +
  scale_y_continuous(limits = c(0,1))+
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
  coord_flip()

p_box

# 转换因子顺序以确保热图的顺序
miss_data$source <- ifelse(grepl("TCGA", miss_data$disease), "TCGA", "GEO")

# miss_data[,'method']=factor(miss_data[,'method'],levels=row.names(mean_values[order(mean_values$mean_value),]))##level和boxplot保持一致
miss_data[,'method']=factor(miss_data[,'method'],levels=rev(c('ReCIDE','bayes','CIBERSORT','DWLS','music')))##level和boxplot保持一致


miss_data$source <- factor(miss_data$source, levels = c("TCGA","GEO"))
miss_data <- miss_data %>%
  arrange(source)
miss_data$disease <- factor(miss_data$disease, levels = unique(miss_data$disease))

heatmap_plot <- ggplot(miss_data, aes(x = disease, y = method, fill = value)) +
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


heatmap_plot


# 
# 


ggsave(
  filename = paste0('./fig5_benchmark/results_and_fig/barplot.pdf'),  # 保存的文件名
  plot = p_box,  # 要保存的图形对象
  width = 3,  # 宽度（默认单位为英寸）
  height = 5,  # 高度
  units = "in",  # 单位：英寸
  dpi = 300,  # 分辨率（每英寸点数）
  device = "pdf"  # 设备类型
)


#
ggsave(
  filename = paste0('./fig5_benchmark/results_and_fig/heatmap.pdf'),  # 保存的文件名
  plot = heatmap_plot,  # 要保存的图形对象
  width = 8,  # 宽度（默认单位为英寸）
  height = 5,  # 高度
  units = "in",  # 单位：英寸
  dpi = 300,  # 分辨率（每英寸点数）
  device = "pdf"  # 设备类型
)


