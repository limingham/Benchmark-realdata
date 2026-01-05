

library(tidyverse)
library(reshape2)
library(ggpubr)
library(pheatmap)


source('~/Benchmark-realdata-main/fig1/results_analysis/readTNBC_realbulk.R')
prd_list_TNBC_realbulk = prd_list



source('~/Benchmark-realdata-main/fig1/results_analysis/readER_realbulk.R')
prd_list_ER_realbulk = prd_list



key_TNBC <- readRDS("~/Benchmark-realdata-main/fig1/data_preprocess/key_TNBC.rds")


key_TNBC = key_TNBC[sort(row.names(key_TNBC)),sort(colnames(key_TNBC))]



key_ER <- readRDS("~/Benchmark-realdata-main/fig1/data_preprocess/key_ER.rds")

key_ER = key_ER[sort(row.names(key_ER)),sort(colnames(key_ER))]

# 


key_all = cbind(key_TNBC,key_ER)

key_all = key_all[ ,colnames(key_all) != 'CID44991']


cor_vec_realbulk_list = list()

for(ii in names(prd_list_TNBC_realbulk)){
  prd_in_TNBC = prd_list_TNBC_realbulk[[ii]]
  prd_in_ER = prd_list_ER_realbulk[[ii]]
  
  
  prd_in_TNBC = prd_in_TNBC[,sort(colnames(prd_in_TNBC))]
  prd_in_ER = prd_in_ER[,sort(colnames(prd_in_ER))]
  prd_in = cbind(prd_in_TNBC,prd_in_ER)
  
  
  prd_in = prd_in[row.names(key_all),colnames(key_all)]
  # prd_in = prd_in_ER[row.names(key_ER),colnames(key_ER)]
  
  cor_in =c()
  for(kk in 1:nrow(prd_in)){
    

      
      
      
      cor_in = c(cor_in,cor(as.numeric(prd_in[kk,colnames(key_all)]),
                            as.numeric(key_all[kk,colnames(key_all)]),method = 'pearson'))
      
    }
  names(cor_in) = row.names(prd_in)
  
  
  # cor_vec_realbulk = c(cor_vec_realbulk,mean(cor_in))
  cor_vec_realbulk_list[[ii]] = cor_in
}

names(cor_vec_realbulk_list) = names(prd_list_TNBC_realbulk)


# cor(cor_vec_psedobulk[c(1,3,4,5,6)],cor_vec_realbulk[c(1,3,4,5,6)])



cor_vec_realbulk_df = as.data.frame(cor_vec_realbulk_list)
cor_vec_realbulk_df[is.na(cor_vec_realbulk_df)] =0

cor1 = cor(t(cor_vec_realbulk_df[,c(1,3,4,5,6)]),method = 'spearman')

diag(cor1) =NA
median(cor1,na.rm = TRUE)

vec_celltypes_2celltypes =cor1[upper.tri(cor1)]


cor2 = cor1[!row.names(cor1) %in% c('MonocyteMacrophage','Endothelial_CXCL12','Cycling_Myeloid','DCs'),
            !colnames(cor1) %in% c('MonocyteMacrophage','Endothelial_CXCL12','Cycling_Myeloid','DCs')]
median(cor2,na.rm = TRUE)

vec_noDP_2noDP =cor2[upper.tri(cor2)]

