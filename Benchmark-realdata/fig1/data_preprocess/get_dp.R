

###Some samples were excluded from DP cell type identification due to the presence of a large number of cell proportions with zero values.

key_TNBC <- readRDS("~/Benchmark-realdata-main/fig2/data_preprocess/key_TNBC_findDP.rds")
key_ER <- readRDS("~/Benchmark-realdata-main/fig2/data_preprocess/key_ER_findDP.rds")


key_TNBC = key_TNBC[sort(row.names(key_TNBC)),]



key_ER = key_ER[sort(row.names(key_ER)),]
# 

p_value = c()
fc_scRNA = c()
for(mm in 1:nrow(key_ER)){
  m1 = as.numeric(key_ER[mm,])
  m2 = as.numeric(key_TNBC[mm,])
  
  # 
  aa = wilcox.test(m1,m2)
  
  p_value =c(p_value, aa[["p.value"]])
}
p_adjust = p.adjust(p_value,method = 'fdr')
names(p_adjust) = row.names(key_ER)

sort(p_adjust)


