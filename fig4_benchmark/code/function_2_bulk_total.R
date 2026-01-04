# disease = 'TCGA-COAD'
# inter1 = intersect(colnames())

get_pvalue <- function(results_list,key_trend){


  
  for (i in 1:length(results_list)){
    df_category1 = as.data.frame(t(results_list[[i]][[1]]))
    df_category2 = as.data.frame(t(results_list[[i]][[2]]))
    
    col_inter = intersect(colnames(df_category1),colnames(df_category2))
    col_inter = intersect(col_inter,row.names(key_trend))

    
    df_category1 = df_category1[,col_inter]
    df_category2 = df_category2[,col_inter]
    
    
    
    df_category1 = df_category1[sort(row.names(df_category1)),sort(colnames(df_category1))]
    df_category2 = df_category2[sort(row.names(df_category2)),sort(colnames(df_category2))]
    

    
    
    p_values <- numeric(length = ncol(df_category1))
    trends <- character(length = ncol(df_category2))
    diff <- numeric(length = ncol(df_category2))
    
    
    
    for (j in seq_along(colnames(df_category1))) {
      col_name <- colnames(df_category1)[j]
      
      
      df_category1_in = as.numeric(df_category1[, col_name])
      df_category1_in = df_category1_in[df_category1_in>0.001]
      
      df_category2_in = as.numeric(df_category2[, col_name])
      df_category2_in = df_category2_in[df_category2_in>0.001]
      

      
      diff[j] <- abs(median(df_category1_in)-median(df_category2_in))
      
      if(length(df_category1_in) <3 | length(df_category2_in)<3){
        
        trends[j] <- '0' 
        p_values[j] = 1
        
        
      }else{
        pvalue <- wilcox.test(df_category1_in, df_category2_in, alternative = "two.sided")$p.value
        p_values[j] <- pvalue
        if (median(df_category1_in) < median(df_category2_in)) {
          trends[j] <- '+'  ##normal smaller than tumor,tumor high
        } else if (median(df_category1_in) > median(df_category2_in)) {
          trends[j] <- '-'
        }else {
          trends[j] <- '0'
        } 
      }}
    


    adjusted_p_values <- p.adjust(p_values, method = "fdr")
    
    
    
    # ?????????????????????
    re_df <- data.frame(
      Celltype = colnames(df_category1),
      PValue = adjusted_p_values,
      Trend = trends,
      Diff=diff
    )
    
    re_df$orig=names(results_list)[i]
    if(i == 1){
      RE_DF = re_df
    }else{
      RE_DF <- rbind(RE_DF,re_df)}
  }
 
  
  RE_DF[RE_DF$Trend=='+','Trend_pvalue'] <- RE_DF[RE_DF$Trend=='+','PValue']
  RE_DF[RE_DF$Trend=='-','Trend_pvalue'] <- 0 - RE_DF[RE_DF$Trend=='-','PValue']
  RE_DF[RE_DF$Trend==0,'Trend_pvalue'] <- 2

  return(RE_DF)
}
