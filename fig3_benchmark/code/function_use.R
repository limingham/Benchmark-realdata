# disease = 'GEO-PAAD'
# celltype_use = row.names(prop)
# criteria ='strict'

function_cal_F1 <- function(results_list,key_trend,disease,
                            criteria = crit){

  
  # i ='music'
  for (i in 1:length(results_list)){
    df_category1 = as.data.frame(t(results_list[[i]][[1]]))
    df_category2 = as.data.frame(t(results_list[[i]][[2]]))
    
    df_category1 <- df_category1 %>%
      # 保留NA占比≤50%的行
      filter(rowSums(is.na(.)) / ncol(.) <= 0.5) %>%
      # 保留NA占比≤50%的列
      select(where(~ sum(is.na(.x)) / length(.x) <= 0.5))
    
    
    df_category2 <- df_category2 %>%
      # 保留NA占比≤50%的行
      filter(rowSums(is.na(.)) / ncol(.) <= 0.5) %>%
      # 保留NA占比≤50%的列
      select(where(~ sum(is.na(.x)) / length(.x) <= 0.5))
    
    
    col_inter = intersect(colnames(df_category1),colnames(df_category2))
    col_inter = intersect(col_inter,row.names(key_trend))
    
    df_category1 = df_category1[,col_inter, drop = FALSE]
    df_category2 = df_category2[,col_inter, drop = FALSE]
    
    
    p_values <- numeric(length = ncol(df_category1))
    trends <- character(length = ncol(df_category2))
    diff <- numeric(length = ncol(df_category2))
    
    for (j in seq_along(colnames(df_category1))) {
      col_name <- colnames(df_category1)[j]
      
      df_category1_in = as.numeric(df_category1[, col_name])
      df_category1_in = df_category1_in[df_category1_in>0.001]
      
      df_category2_in = as.numeric(df_category2[, col_name])
      df_category2_in = df_category2_in[df_category2_in>0.001]
      
      
      if(length(df_category1_in) <3 | length(df_category2_in)<3){
        
        trends[j] <- '0' 
        p_values[j] = 1
        
        
      }else{
        

        pvalue <- wilcox.test(df_category1_in, df_category2_in, alternative = "two.sided")$p.value
        p_values[j] <- pvalue
        diff[j] <- abs(median(df_category1_in)-median(df_category2_in))
        if (diff[j]==0){
          diff[j] <- abs(mean(df_category1_in)-mean(df_category2_in))
        }

        
        
        if (median(df_category1_in)<median(df_category2_in)) {
          trends[j] <- '+'  ##normal smaller than tumor,tumor high
        } else if (median(df_category1_in)>median(df_category2_in)) {
          trends[j] <- '-'
        } else {
          trends[j] <- '0' 
          p_values[j] = 1
          
        }
      }
      
      
    }
    
    
    # adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
    if (names(results_list)[i]=='ReCIDE'){
      adjusted_p_values <- p.adjust(p_values, method = "fdr")
    }else{
      adjusted_p_values <- p.adjust(p_values, method = "fdr")
    }
    # 
    # # ??????????????????

    # }
    
    
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
  

  
  
  
  
  df_result = data.frame()
  for (methods in c("CIBERSORT","DWLS","ReCIDE","bayes","music")) {
    # methods = "music"
    RE_DF$PValue[is.nan(RE_DF$PValue)] <- 1
    df_test = RE_DF[RE_DF[,'orig'] == methods,]
    row.names(df_test) <- df_test[,1]
    df_test[is.na(df_test)]=1
    
    key_trend$PValue = key_trend$PValue
    df2_test = key_trend #scRNA
    
    ############################
    
    df2_test_sig1 <- df2_test[df2_test$PValue<0.05,]
    if(length(intersect(df2_test_sig1$Celltype,row.names(df2_test))) > 0){
      c_ = intersect(row.names(df2_test),row.names(df_test))

      
      
      
      df_test = df_test[c_,]
      df2_test = df2_test[c_,]
      

      
      df_test_sig <- df_test[df_test$PValue<0.05,]
      df2_test_sig <- df2_test[df2_test$PValue<0.05,]
      
      df_test_nsig <- df_test[df_test$PValue>0.05,]
      df2_test_nsig <- df2_test[df2_test$PValue>0.05,]
      
      trend_consist <- df_test[which(df_test$Trend==df2_test$Trend),]##趋势一致的细胞类型
      
      sig_2 <- intersect(df_test_sig$Celltype,df2_test_sig$Celltype)##2个都显著的细胞类型
      sig_2_consis <- intersect(trend_consist$Celltype,sig_2)##2个都显著且同向细胞类型
      

      
      
      ###########################| ((df2_test[row.names(df_test_sig),'PValue'] == 0.5))
      if(criteria == 'permissive'){
        AR_num=sum((df_test_sig[,'Trend']==df2_test[row.names(df_test_sig),'Trend' ]))##宽松版本
      }else if(criteria == 'conservative'){
        AR_num=length(sig_2_consis)##严格版本
      }
      

      
      df_result[methods,'RR'] <- length(sig_2_consis)/dim(df2_test_sig)[1]##RR 
      df_result[methods,'Jaccard'] <-length(sig_2_consis)#bulk_scRNA_trend2/dim(df_test_sig)[1] ##AR
      df_result[methods,'AR_num'] <-AR_num
      df_result[methods,'all_num']<- dim(df_test_sig)[1] 
      df_result[methods,'trend_match_ratio_correctCD'] <- -5
      df_result[methods,'AR'] <-AR_num/dim(df_test_sig)[1]##和单细胞趋势相同的
      df_result[methods,'F1'] <- 2*(df_result[methods,'AR'] *df_result[methods,'RR'])/(df_result[methods,'AR']+df_result[methods,'RR'])
      
      df_result$disease <- disease
      # list_tt[[disease]] <- df_result
    }else{
      df_result[methods,'RR'] <- -5
      df_result[methods,'Jaccard'] <- -5
      df_result[methods,'AR_num'] <--5
      df_result[methods,'all_num']<- -5
      df_result[methods,'trend_match_ratio_correctCD'] <- -5
      df_result[methods,'AR'] <--5
      df_result[methods,'F1'] <- -5
      
      df_result$disease <- disease
      # list_tt[[disease]] <- df_result
    }}
  return(df_result)
  
}