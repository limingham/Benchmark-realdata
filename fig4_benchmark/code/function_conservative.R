# # # # #
#   disease = 'BRCA'
# celltype_use = intersect(row.names(prop),row.names(key_trend))

compare_bulk2bulk <- function(df,df2,disease,key_trend){
  # check_list <- list()
  # key_trend <- readRDS(key_trend_dir)
  row.names(key_trend) <- key_trend[,1]
  # key_trend = key_trend[celltype_use,]
  key_trend_signifcance <- key_trend[key_trend$PValue < 0.05,]
  
  
  # key_trend_signifcance = key_trend_signifcance[celltype_use,]
  df[is.na(df)]=1
  df2[is.na(df2)]=1
  
  # row.names(df) = paste0(df[,1],df[,5])
  # row.names(df2) = paste0(df2[,1],df2[,5])
  # df = df[sort(row.names(df)),]
  # df2 = df2[sort(row.names(df2)),]
  # 
  row.names(df) = paste0(df[,1],df[,5])
  row.names(df2) = paste0(df2[,1],df2[,5])
  df = df[sort(row.names(df)),]
  df2 = df2[sort(row.names(df2)),]
  
  
  df2_plus = setdiff(row.names(df),row.names(df2))
  df1_plus = setdiff(row.names(df2),row.names(df))
  
  df2[df2_plus,] = df[df2_plus,]
  df2[df2_plus,'PValue'] = 1
  df2[df2_plus,'Trend'] = '0'
  
  df[df1_plus,] = df2[df1_plus,]
  df[df1_plus,'PValue'] = 1
  df[df1_plus,'Trend'] = '0'
  
  df = df[sort(row.names(df)),]
  df2 = df2[sort(row.names(df2)),]
  
  
  combined_df = cbind(df,df2)
  colnames(combined_df) <- c('Celltype','p_GEO', 'GEO_trend','GEO_diff','orig','GEO_Trend_pvalue',
                             'Celltype_TCGA','p_TCGA', 'TCGA_trend','TCGA_diff','orig_TCGA','TCGA_Trend_pvalue')
  

  # combined_df = combined_df[combined_df$Celltype %in% celltype_use,]
  combined_df_ = combined_df[combined_df[,'p_GEO'] < 0.05 & combined_df[,'p_TCGA'] < 0.05,]
  combined_df_1 = combined_df[combined_df[,'p_GEO'] < 0.05 | combined_df[,'p_TCGA'] < 0.05,]
  
  # 
  # val_df = rbind(df,df2)
  # val_df = val_df[val_df$PValue<0.05,]
  ###质控

  df_result <- data.frame()
  for (j in unique(df$orig)){
    # j = 'ReCIDE'
    # j = 'bayes'
    TNBC_df_in = combined_df_[combined_df_[,'orig'] %in% j,]
    TNBC_df_in_1 = combined_df_1[combined_df_1[,'orig'] %in% j,]
    
    
    row.names(TNBC_df_in) = TNBC_df_in$Celltype
    
    # cell_inter = intersect(row.names(df_inter),celltype_use)
    # df_inter = TNBC_df_in[cell_inter,]
    
    df_inter = TNBC_df_in[(TNBC_df_in[,"p_TCGA"] <0.05 & TNBC_df_in[,"p_GEO"]<0.05) & 
                            (TNBC_df_in[,"TCGA_trend"] == TNBC_df_in[,"GEO_trend"])& 
                            (TNBC_df_in[,"TCGA_trend"] != "0"),]
    
    row.names(df_inter) = df_inter$Celltype
    
    
    # cell_inter = intersect(row.names(df_inter),celltype_use)
    # df_inter_sc = df_inter[cell_inter,]
    # 
    df_inter_sc = df_inter
    
    
    
    # key_trend_signifcance = key_trend[key_trend$PValue<0.05,]
    # key_trend_signifcance = key_trend_signifcance[key_trend_signifcance$Celltype %in% celltype_use,]
    ###筛选出单细胞中显著且拿来用的
    
    
    # key_trend_signifcance
    
    val_celltype = paste0(df_inter_sc$Celltype,'-',df_inter_sc$TCGA_trend)
    
    trend_in_key_sign = paste0(key_trend_signifcance$Celltype,'-',key_trend_signifcance$Trend)
    rr_num <- sum(val_celltype %in% trend_in_key_sign)
    
    trend_in_key = paste0(key_trend$Celltype,'-',key_trend$Trend)
    ar_num <- sum(val_celltype %in% trend_in_key)
    
    # 


    df_result[j,'Jaccard'] <-nrow(df_inter)/nrow(TNBC_df_in_1)
    df_result[j,'cell_num'] <-nrow(df_inter)
    df_result[j,'RR'] <- rr_num/nrow(key_trend_signifcance)
    df_result[j,'AR'] <- ar_num/nrow(df_inter_sc)
    df_result[j,'F1'] <-2*df_result[j,'RR']*df_result[j,'AR']/(df_result[j,'RR']+df_result[j,'AR'])
    
    df_result$disease <- disease
    

    
    
    
  }
  return(df_result)
}



