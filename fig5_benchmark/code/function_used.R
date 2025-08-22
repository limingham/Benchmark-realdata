
read_deconvolurtion_results <- function(results_list,key_trend,disease){

  
  
  
  
  if(!disease %in% c('TCGA-NSCLC','GEO-NSCLC')){
    for (r1 in 1:length(results_list)) {
      for (r2 in 1:length(results_list[[r1]])) {

        cellaa = intersect(row.names(results_list[[r1]][[r2]]),row.names(key_trend))
        results_list[[r1]][[r2]] = results_list[[r1]][[r2]][cellaa,]
      }

    }}else if(disease %in% c('TCGA-NSCLC','GEO-NSCLC')){
          for (r1 in 1:length(results_list)) {
            for (r2 in 1:length(results_list[[r1]])) {
              if(r2 == 1){
                key_trend1 =  key_trend[['LUSC']]
                cellbb = intersect(row.names(results_list[[r1]][[r2]]),row.names(key_trend1))
                results_list[[r1]][[r2]] = results_list[[r1]][[r2]][cellbb,]

              }else{
                key_trend2 =  key_trend[['LUAD']]

                cellbb = intersect(row.names(results_list[[r1]][[r2]]),row.names(key_trend2))
                results_list[[r1]][[r2]] = results_list[[r1]][[r2]][cellbb,]
              }

            }

          }

    }
  
  
  # 
  return(results_list) 
}






func_calculate_hr_TCGA<- function(gene_exp2,Clinical){
  
  gene_exp<-as.data.frame(gene_exp2)
  ###这里需要注意，gene_exp是一列数，而不是一行数
  
  ##<=的话小于的被cut掉，成为TRUE
  cutoff= (gene_exp2[,1]<=quantile(as.numeric(gene_exp2[,1]))[3])#以中位数为界分为高表达及低表达组
  if(all(cutoff) == 'TRUE' | length(cutoff)<15){
    return(c(1,1,0))
  }else{
    gene_exp$age <- as.numeric(Clinical[row.names(gene_exp), 'age_at_index'])
    gene_exp$gender <- as.factor(Clinical[row.names(gene_exp), 'gender'])
    gene_exp$race <- as.factor(Clinical[row.names(gene_exp), 'race'])
    gene_exp$primary_diagnosis <- as.factor(Clinical[row.names(gene_exp), 'primary_diagnosis'])
    gene_exp$tissue_or_organ_of_origin <- as.factor(Clinical[row.names(gene_exp), 'tissue_or_organ_of_origin'])
    gene_exp$prior_malignancy <- as.factor(Clinical[row.names(gene_exp), 'prior_malignancy'])
    gene_exp$prior_treatment <- as.factor(Clinical[row.names(gene_exp), 'prior_treatment'])
    gene_exp$treatments_pharmaceutical_treatment_or_therapy <- as.factor(Clinical[row.names(gene_exp), 'treatments_pharmaceutical_treatment_or_therapy'])
    # gene_exp$ajcc_pathologic_stage <- as.factor(Clinical[row.names(gene_exp), 'ajcc_pathologic_stage'])
    gene_exp$treatments_radiation_treatment_or_therapy <- as.factor(Clinical[row.names(gene_exp), 'treatments_radiation_treatment_or_therapy'])
    gene_exp$alcohol_history=as.factor(Clinical[row.names(gene_exp), 'alcohol_history'])
    
    
    # aaa=as.numeric(Clinical$years_smoked)
    # gene_exp$smoke <- as.numeric(Clinical[row.names(gene_exp), 'years_smoked'])
    # levels(aaa) <- c(levels(aaa), "aa")
    # 
    # # # 再替换 NA
    # aaa[is.na(aaa)] <- 0
    # gene_exp$years_smoked <- aaa

    
    
    # gene_exp$stage <- as.factor(Clinical[row.names(gene_exp), 'stage'])
    gene_exp$cutoff1 = cutoff
    
    surstat <- as.numeric(Clinical[row.names(gene_exp), 'vital_status'] == "Dead")
    surstat[is.na(surstat)] = 0
    
    surtime = c()
    for (kk in 1:length(surstat)) {
      if(surstat[kk] !=1){
        surtime[kk] = Clinical[row.names(gene_exp)[kk],'days_to_last_follow_up']
      }else{
        surtime[kk] = Clinical[row.names(gene_exp)[kk],'days_to_death']
      }
      
    }

    
    # 'cell_prop',
    colnames(gene_exp)<-c('cell_prop','age','gender','race','primary_diagnosis', 
                           'tissue_or_organ_of_origin','prior_malignancy', 'prior_treatment',
                          'treatments_pharmaceutical_treatment_or_therapy',
                          'treatments_radiation_treatment_or_therapy','alcohol_history','cutoff1')
    
    surv_object <- Surv(time = surtime, event = surstat)
    
    # 进行 Cox 回归分析
    # cox_model <- coxph(surv_object ~ cell_prop + age + gender, data = gene_exp)
    
    
    filter_vars <- function(data, vars) {
      keep <- sapply(vars, function(var) {
        if (is.factor(data[[var]]) || is.character(data[[var]])) {
          length(unique(na.omit(data[[var]]))) >= 2  # 至少2个非缺失水平
        } else {
          TRUE  # 连续变量保留
        }
      })
      vars[keep]
    }
    # 'cell_prop',,'cutoff1'
    
    
    all_vars <- c('age','gender','race','primary_diagnosis', 
                  'tissue_or_organ_of_origin','prior_malignancy',
                  'treatments_pharmaceutical_treatment_or_therapy',
                  'treatments_radiation_treatment_or_therapy','alcohol_history','cutoff1')
    
    valid_vars <- filter_vars(gene_exp, all_vars)
    
    
    
    
    formula_str <- paste("surv_object ~", paste(valid_vars, collapse = " + "))
    

###
      # formula_str = paste("surv_object ~", 'cutoff1')
###
    
    cox_model <- coxph(as.formula(formula_str), data = gene_exp)
    
    
    # 获取 Cox 回归模型的详细结果
    summary_cox <- summary(cox_model)
    aa = as.data.frame(summary_cox[["coefficients"]])
    # 提取危险比（HR）
    hr_values <- summary_cox$coefficients['cutoff1TRUE', "exp(coef)"]
    hr_p <- summary_cox$coefficients["cutoff1TRUE", "Pr(>|z|)"]
    # hr_values <- summary_cox$coefficients['cell_prop', "exp(coef)"]
    # hr_p <- summary_cox$coefficients["cell_prop", "Pr(>|z|)"]
    # 
    
    return(c(hr_values,hr_p,median(as.numeric(gene_exp2[,1]))))
    
  }
}








func_calculate_hr_GEO<- function(gene_exp2,Clinical,dise){
  
  
  if(dise == 'GEO-HNSC_neg'){
    # Clinical <- readRDS("/home/syq/DEC_benmark/TCGA/processed/HNSC/Bulk/GSE65858/meta.rds")
    

    Clinical[,"age"] = as.numeric(Clinical[,"age:ch1"])
    Clinical[,"treatment"] = as.factor(Clinical$`treatment:ch1`)

    Clinical[,"characteristics_ch1.1"] = as.factor(Clinical[,"characteristics_ch1.1"])
    Clinical[,"characteristics_ch1.7"] = as.factor(Clinical[,"characteristics_ch1.7"])
    
    Clinical[,"characteristics_ch1.3"] = as.factor(Clinical[,"characteristics_ch1.3"])
    Clinical[,"characteristics_ch1.5"] = as.factor(Clinical[,"characteristics_ch1.5"])
    # Clinical[,"characteristics_ch1.21"] = as.factor(Clinical[,"characteristics_ch1.21"])
    Clinical[,"characteristics_ch1.48"] = as.factor(Clinical[,"characteristics_ch1.48"])
    Clinical[,"characteristics_ch1.19"] = as.factor(Clinical[,"characteristics_ch1.19"])
    Clinical[,"dna_rna"] = as.factor(Clinical[,"hpv16_dna_rna:ch1"])
    
    # confounders = c("characteristics_ch1.1","age",
    #                 "characteristics_ch1.7","characteristics_ch1.3",
    #                 "characteristics_ch1.5","characteristics_ch1.48",
    # "characteristics_ch1.19","dna_rna","treatment")
    confounders = c("characteristics_ch1.1","age",
                    "characteristics_ch1.7",
                    "characteristics_ch1.3",
                    "characteristics_ch1.5",
                    "characteristics_ch1.48",
                    # "characteristics_ch1.19",
                    # "dna_rna",
                    "treatment"
                    )
    # #
    
    
  }else if(dise == 'GEO-LUAD'){
    Clinical[,"age"] = as.numeric(Clinical[,"age at surgery:ch1"])
    Clinical[,"had_adjuvant_chemo"] = Clinical[,"had_adjuvant_chemo:ch1"]
    Clinical[,"had_adjuvant_chemo"] = as.factor(Clinical[,"had_adjuvant_chemo"])
    Clinical[,"characteristics_ch1.4"]  = as.factor(Clinical[,"characteristics_ch1.4"])
    Clinical[,"pat_stage"] = Clinical[,"final.pat.stage:ch1"]
    
    confounders = c("characteristics_ch1.4","age","had_adjuvant_chemo")  
    # confounders = c()
  }else if(dise == 'GEO-COAD'){
    
    
    Clinical[,'age'] = as.numeric(Clinical[,"age.at.diagnosis (year):ch1"])
    
    # 分别转换每一列
    Clinical[,"tp53.mutation_ch1"] <- as.factor(Clinical[,"tp53.mutation:ch1"])
    Clinical[,"kras.mutation_ch1"] <- as.factor(Clinical[,"kras.mutation:ch1"])
    Clinical[,"braf.mutation_ch1"] <- as.factor(Clinical[,"braf.mutation:ch1"])
    
    Clinical[,"rna_ext"] <- as.factor(Clinical[,"characteristics_ch1.32"])
    # Clinical[,"cit"] <- as.factor(Clinical[,"characteristics_ch1.30"])
    
    Clinical[,"mmr.status_ch1"] <- as.factor(Clinical[,"mmr.status:ch1"])
    Clinical[,'tum_loc'] = as.factor(Clinical[,"tumor.location:ch1"])
    Clinical[,'characteristics_ch1.2'] = as.factor(Clinical[,"characteristics_ch1.2"])
    # Clinical[,'d_v'] = as.factor(Clinical[,"dataset:ch1"])
    # Clinical[,'braf'] = as.factor(Clinical[,"characteristics_ch1.26"])
    Clinical[,'chem'] = as.factor(Clinical[,"characteristics_ch1.9"])
    Clinical[,'characteristics_ch1.16'] = as.factor(Clinical[,"characteristics_ch1.16"])
    Clinical[,'characteristics_ch1.17'] = as.factor(Clinical[,"characteristics_ch1.17"])
    Clinical[,'cit'] = as.factor(Clinical[,"cit.molecularsubtype:ch1"])
    Clinical[,'rfs'] = as.numeric(Clinical[,"rfs.event:ch1"])
    
    # "characteristics_ch1.8", ,"kras.mutation_ch1","chem" ,'characteristics_ch1.17'
    confounders = c("characteristics_ch1.2","tum_loc","tp53.mutation_ch1","age","kras.mutation_ch1",'braf.mutation_ch1',
                    "mmr.status_ch1",'characteristics_ch1.17','characteristics_ch1.16')
    # #
    # confounders = c("characteristics_ch1.2","tp53.mutation_ch1","age","kras.mutation_ch1",
    #                 "mmr.status_ch1")
 
    # 
  }else if(dise == 'GEO-LUSC'){
    
    Clinical[,"age"] = as.numeric(Clinical[,"age at surgery:ch1"])
    Clinical[,"had_adjuvant_chemo"] = Clinical[,"had_adjuvant_chemo:ch1"]
    
    Clinical[,"had_adjuvant_chemo"] = as.factor(Clinical[,"had_adjuvant_chemo"])
    Clinical[,"characteristics_ch1.4"]  = as.factor(Clinical[,"characteristics_ch1.4"])
    
    Clinical[,"pat_stage"] = Clinical[,"final.pat.stage:ch1"]
    confounders = c("characteristics_ch1.4","age","had_adjuvant_chemo")  
    
  }else if(dise == 'GEO-PAAD'){
    Clinical[,"stroma"] = as.numeric(Clinical[,"stroma_subtype_0na_1low_2normal_3activated.ch2"])
    Clinical[,"tumor_subtype_0na_1classical_2basal2"] = as.numeric(Clinical[,"tumor_subtype_0na_1classical_2basal.ch2"])
    #
    # confounders=c('stroma','tumor_subtype_0na_1classical_2basal2')
    confounders=c()
  }
  
  gene_exp<-as.data.frame(gene_exp2)
  ###这里需要注意，gene_exp是一列数，而不是一行数
  
  ##<=的话小于的被cut掉，成为TRUE
  cutoff= (gene_exp[,1]<=quantile(as.numeric(gene_exp[,1]))[3])#以中位数为界分为高表达及低表达组
  if(all(cutoff) == 'TRUE' | length(cutoff)<15){
    return(c(1,1,0))
  }else{
    
    gene_exp[,confounders] = Clinical[row.names(gene_exp), confounders]
    
    # gene_exp$cutoff1 = cutoff
    
    surstat <- Clinical[row.names(gene_exp), 'OS']
    surtime =Clinical[row.names(gene_exp), 'OS.time']
    
    gene_exp$cutoff1 = cutoff
    colnames(gene_exp)<-c('cell_prop',confounders,'cutoff1')
    
    surv_object <- Surv(time = surtime, event = surstat)
    
    # 进行 Cox 回归分析
    # cox_model <- coxph(surv_object ~ cell_prop + age + gender, data = gene_exp)
    
    
    filter_vars <- function(data, vars) {
      keep <- sapply(vars, function(var) {
        if (is.factor(data[[var]]) || is.character(data[[var]])) {
          length(unique(na.omit(data[[var]]))) >= 2  # 至少2个非缺失水平
        } else {
          TRUE  # 连续变量保留
        }
      })
      vars[keep]
    }
    
    # colnames(colnames)[1] = 'cell_prop'
    all_vars <- colnames(gene_exp)[2:ncol(gene_exp)]
    # all_vars <- colnames(gene_exp)[1:(ncol(gene_exp)-1)]
    
    valid_vars <- filter_vars(gene_exp, all_vars)
    
    
    formula_str <- paste("surv_object ~", paste(valid_vars, collapse = " + "))
    cox_model <- coxph(as.formula(formula_str), data = gene_exp)
    
    
    # 获取 Cox 回归模型的详细结果
    summary_cox <- summary(cox_model)
    
    # 提取危险比（HR）
    hr_values <- summary_cox$coefficients['cutoff1TRUE', "exp(coef)"]
    hr_p <- summary_cox$coefficients["cutoff1TRUE", "Pr(>|z|)"]
    
    # hr_values <- summary_cox$coefficients['cell_prop', "exp(coef)"]
    # hr_p <- summary_cox$coefficients["cell_prop", "Pr(>|z|)"]
    
    return(c(hr_values,hr_p,median(as.numeric(gene_exp2[,1]))))
  }
  
}
