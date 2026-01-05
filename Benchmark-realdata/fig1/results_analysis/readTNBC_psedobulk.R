
for (ll in 1:10) {
  

file_list = list.files('~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk/')

###########recide
###########
###########
###########
recide_= readRDS(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[1],"/results_ReCIDE.rds"))
recide_df <- data.frame(recide_[,ll])
rownames(recide_df) <- rownames(recide_)  
 

for(mm in 2:length(file_list)){
  recide_in_= readRDS(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[mm],"/results_ReCIDE.rds"))
  recide_in_df <- data.frame(recide_in_[,ll])
  rownames(recide_in_df) <- rownames(recide_in_) 
  recide_df <- merge(
    recide_df, recide_in_df, 
    by = "row.names",    
    all = TRUE           
  )
  rownames(recide_df) <- recide_df$Row.names
  recide_df$Row.names <- NULL
  recide_df[is.na(recide_df)] <- 0
}

colnames(recide_df) = file_list





###########DWLS
###########
###########
###########
DWLS1= readRDS(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[1],"/results_DWLS.rds"))

 

DWLS_df <- data.frame(DWLS1[,ll])
rownames(DWLS_df) <- rownames(DWLS1)    

for(mm in 2:length(file_list)){
  
  DWLS_in= readRDS(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[mm],"/results_DWLS.rds"))
  DWLS_in_df <- data.frame(DWLS_in[,ll])
  rownames(DWLS_in_df) <- rownames(DWLS_in)    
  DWLS_df <- merge(
    DWLS_df, DWLS_in_df, 
    by = "row.names",    
    all = TRUE           
  )
  rownames(DWLS_df) <- DWLS_df$Row.names
  DWLS_df$Row.names <- NULL
  DWLS_df[is.na(DWLS_df)] <- 0
}

colnames(DWLS_df) = file_list

DWLS_df = DWLS_df[row.names(DWLS_df)!="1",]



###########music
###########
###########
###########
music1= readRDS(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[1],"/results_music.rds"))
music1 = as.data.frame(t(music1))
 

music_df <- data.frame(music1[,ll])
rownames(music_df) <- rownames(music1)    

for(mm in 2:length(file_list)){
  music_in= readRDS(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[mm],"/results_music.rds"))
  music_in = as.data.frame(t(music_in))
  music_in_df <- data.frame(music_in[,ll])

   
  
  rownames(music_in_df) <- rownames(music_in)    
  music_df <- merge(
    music_df, music_in_df, 
    by = "row.names",    
    all = TRUE           
  )
  rownames(music_df) <- music_df$Row.names
  music_df$Row.names <- NULL
  music_df[is.na(music_df)] <- 0
}

colnames(music_df) = file_list





###########bisque
###########
###########
bisque1= readRDS(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[1],"/results_bisque.rds"))
bisque1 = as.data.frame(bisque1)
 

bisque_df <- data.frame(bisque1[,ll])
rownames(bisque_df) <- rownames(bisque1)    

for(mm in 2:length(file_list)){
  bisque_in= readRDS(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[mm],"/results_bisque.rds"))
  bisque_in = as.data.frame(bisque_in)
  bisque_in_df <- data.frame(bisque_in[,ll])

   
  
  rownames(bisque_in_df) <- rownames(bisque_in)    
  bisque_df <- merge(
    bisque_df, bisque_in_df, 
    by = "row.names",    
    all = TRUE           
  )
  rownames(bisque_df) <- bisque_df$Row.names
  bisque_df$Row.names <- NULL
  bisque_df[is.na(bisque_df)] <- 0
}

colnames(bisque_df) = file_list






###########bayes
###########
###########
bayes1= readRDS(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[1],"/results_bayes.rds"))
bayes1 = as.data.frame(t(bayes1))
 

bayes_df <- data.frame(bayes1[,ll])
rownames(bayes_df) <- rownames(bayes1)    

for(mm in 2:length(file_list)){
  bayes_in= readRDS(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[mm],"/results_bayes.rds"))
  bayes_in = as.data.frame(t(bayes_in))
  bayes_in_df <- as.data.frame(bayes_in[,ll])

  rownames(bayes_in_df) <- rownames(bayes_in)    
  bayes_df <- merge(
    bayes_df, bayes_in_df, 
    by = "row.names",    
    all = TRUE           
  )
  rownames(bayes_df) <- bayes_df$Row.names
  bayes_df$Row.names <- NULL
  bayes_df[is.na(bayes_df)] <- 0
}

colnames(bayes_df) = file_list




###########CIBERSORT
###########
###########
CIBERSORT1= read.table(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[1],"/results_CIBERSORT.txt"),sep = '\t',row.names = 1,header = TRUE)
CIBERSORT1 = as.data.frame(t(CIBERSORT1))
 
CIBERSORT1 = CIBERSORT1[1:(nrow(CIBERSORT1)-3),]
CIBERSORT_df <- data.frame(CIBERSORT1[,ll])
rownames(CIBERSORT_df) <- rownames(CIBERSORT1)    

for(mm in 2:length(file_list)){
  CIBERSORT_in= read.table(paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/TNBC_pseudobulk//",file_list[mm],"/results_CIBERSORT.txt"),sep = '\t',row.names = 1,header = TRUE)
  CIBERSORT_in = as.data.frame(t(CIBERSORT_in))
  CIBERSORT_in = CIBERSORT_in[1:(nrow(CIBERSORT_in)-3),]
  
   
  
 CIBERSORT_in_df <- data.frame(CIBERSORT_in[,ll])

  
  rownames(CIBERSORT_in_df) <- rownames(CIBERSORT_in)    
  CIBERSORT_df <- merge(
    CIBERSORT_df, CIBERSORT_in_df, 
    by = "row.names",    
    all = TRUE           
  )
  rownames(CIBERSORT_df) <- CIBERSORT_df$Row.names
  CIBERSORT_df$Row.names <- NULL
  CIBERSORT_df[is.na(CIBERSORT_df)] <- 0
}

colnames(CIBERSORT_df) = file_list

prd_list = list()
prd_list[['ReCIDE']] = recide_df
prd_list[['bisque']] = bisque_df
prd_list[['music']] = music_df
prd_list[['bayes']] = bayes_df
prd_list[['CIBERSORT']] = CIBERSORT_df
prd_list[['DWLS']] = DWLS_df

saveRDS(prd_list,file = paste0('~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pseudo_TNBC',ll,'.rds'))






key_TNBC <- readRDS("~/Benchmark-realdata-main/fig1/data_preprocess/key_TNBC.rds")

mae_list = list()

rmse_vec_all =c()
rmse_all_list =list()

# mm = names(prd_list)[1]
for(mm in names(prd_list)){
  rmse_vec =c()
  prd_in = prd_list[[mm]]
  # key_TNBC = key_TNBC[,colnames(prd_in)]
  names_add = setdiff(row.names(key_TNBC),row.names(prd_in))
  prd_in[names_add,] = 0
  prd_in = prd_in[sort(row.names(key_TNBC)),sort(colnames(prd_in))]
  key_TNBC = key_TNBC[sort(row.names(key_TNBC)),sort(colnames(prd_in))]
  
  rmse_vec = c()
  for(cc in colnames(prd_in)){

    
    results0 = as.numeric(prd_in[,cc])

      
      
      rmse_vec = c(rmse_vec,cor(as.numeric(prd_in[,cc]),
                                as.numeric(key_TNBC[,cc]),method = 'pearson'))
    
    
    
  }
  
  rmse_vec_all =c(rmse_vec_all,mean(rmse_vec))
  names(rmse_vec) = colnames(prd_in)
  rmse_all_list[[mm]] = rmse_vec
  
}


saveRDS(rmse_all_list,file = paste0('~/Benchmark-realdata-main/fig1/results_analysis/TNBC_pseudobulk/pcc_pseudo_TNBC',ll,'.rds'))






}
