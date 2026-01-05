library(Seurat)
library(fastSave)


bulk_data <- readRDS("~/Benchmark-realdata-main/fig1/data_preprocess/psedobulk.rds")
sc_ER = readRDS.lbzip2('~/Benchmark-realdata-main/fig1/data_preprocess/seurat_use.rdsFS')
sc_ER = subset(sc_ER,subtype == 'ER+')


inter_name = intersect(colnames(bulk_data),unique(sc_ER@meta.data$orig.ident))


for(sample in inter_name){
  
  
  ER_bulk = as.data.frame(bulk_data[[sample]])
  # ER_bulk[,2] = ER_bulk[,1]
  # row.names(ER_bulk) = row.names(bulk_data)
  seurat_ER = subset(sc_ER,orig.ident != sample)
  # ER_bulk[] <- lapply(ER_bulk, function(x) as.integer(as.character(x)))
  
  
  results_path <- paste0("~/Benchmark-realdata-main/fig1/deconvolution_results/ER_pseudobulk/",sample)
  if (dir.exists(results_path)) {
  } else {
    dir.create(results_path)
  }
  source("~/Benchmark-realdata-main/deconvolution_code/ReCIDE_all_function_beta.R")
  source("~/Benchmark-realdata-main/fig1/deconvolution_code/run_all_methods.R")
  
  # 
  # 
  seurat_ER@meta.data[["celltype_use"]] <- str_replace_all(seurat_ER@meta.data[["celltype_use"]], " ", "_")
  seurat_ER@meta.data[["celltype_use"]] <- str_replace_all(seurat_ER@meta.data[["celltype_use"]], "\\/", "_")
  seurat_ER@meta.data[["celltype_use"]] <- str_replace_all(seurat_ER@meta.data[["celltype_use"]], "\\.", "_")
  seurat_ER@meta.data[["celltype_use"]] <- str_replace_all(seurat_ER@meta.data[["celltype_use"]], "\\(", "")
  seurat_ER@meta.data[["celltype_use"]] <- str_replace_all(seurat_ER@meta.data[["celltype_use"]], "\\)", "")
  seurat_ER@meta.data[["celltype_use"]] <- str_replace_all(seurat_ER@meta.data[["celltype_use"]], ",", "_")
  seurat_ER@meta.data[["celltype_use"]] <- str_replace(seurat_ER@meta.data[["celltype_use"]], "-", "")
  seurat_ER@meta.data[["celltype_use"]] <- str_replace(seurat_ER@meta.data[["celltype_use"]], "\\+", "")
  seurat_ER@meta.data[["celltype_use"]] <- str_replace(seurat_ER@meta.data[["celltype_use"]], "_III", "3")
  seurat_ER@meta.data[["celltype_use"]] <- str_replace(seurat_ER@meta.data[["celltype_use"]], "_II", "2")
  seurat_ER@meta.data[["celltype_use"]] <- str_replace(seurat_ER@meta.data[["celltype_use"]], "_I", "1")
  seurat_ER@meta.data[["celltype_use"]] <- str_replace_all(seurat_ER@meta.data[["celltype_use"]], "\\-", "_")
  
  
  
  celltype1 = "celltype_use"
  subject1 = "orig.ident"
  
  dir_ref_recide = paste0(results_path,'/ref_ReCIDE.rds')
  dir_results_recide = paste0(results_path,'/results_ReCIDE.rds')
  dir_results_bisque = paste0(results_path,'/results_bisque.rds')
  dir_results_music = paste0(results_path,'/results_music.rds')
  dir_results_bayes = paste0(results_path,'/results_bayes.rds')
  dir_ref_CIBERSORT = paste0(results_path,'/ref_CIBERSORT.rds')
  dir_results_CIBERSORT = paste0(results_path,'/results_CIBERSORT.txt')
  dir_ref_DWLS = paste0(results_path,'/ref_DWLS.rds')
  dir_results_DWLS = paste0(results_path,'/results_DWLS.rds')

  ######run_ReCIDE
  #
  source('~/Benchmark-realdata-main/fig1/deconvolution_code/run_all_methods.R')
  func_run_ReCIDE(SC_ref = seurat_ER,
                  EXP_df = ER_bulk,
                  celltype = celltype1,
                  subject = subject1,
                  dir_ref = dir_ref_recide,
                  dir_results = dir_results_recide)


  ######run_bisque
  source('~/Benchmark-realdata-main/fig1/deconvolution_code/run_all_methods.R')

  func_run_bisque(ref_seurat = seurat_ER,
                  bulkdata = ER_bulk,
                  celltype = celltype1,
                  subject = subject1,
                  dir_results = dir_results_bisque)
  gc()

  ######run_music
  source('~/Benchmark-realdata-main/fig1/deconvolution_code/run_all_methods.R')

  func_run_music(ref_seurat = seurat_ER,
                 EXP_df = ER_bulk,
                 celltype = celltype1,
                 subject = subject1,
                 dir_results = dir_results_music)
  gc()




  ######run_CIBERSORT
  source('~/Benchmark-realdata-main/fig1/deconvolution_code/run_all_methods.R')

  func_run_CIBERSORT(combined.data = seurat_ER,
                     bulk.mtx = ER_bulk,
                     celltype = celltype1,
                     subject = subject1,
                     dir_ref = dir_ref_CIBERSORT,
                     dir_results = dir_results_CIBERSORT)
  gc()


  ######run_DWLS
  source('~/Benchmark-realdata-main/fig1/deconvolution_code/run_all_methods.R')

  func_run_DWLS(scdata_test = seurat_ER,
                bulkdata = ER_bulk,
                celltype = celltype1,
                dir_DWLS_ref = dir_ref_DWLS,
                dir_DWLS_results = dir_results_DWLS)
  gc()

  # 
  ######run_bayes
  source('~/Benchmark-realdata-main/fig1/deconvolution_code/run_all_methods.R')
  
  func_run_bayes(ref_seurat = seurat_ER,
                 bulk.mtx = ER_bulk,
                 celltype = celltype1,
                 subject = subject1,
                 dir_results = dir_results_bayes,
                 key1 = 'Epithelial')
  # rm(list = ls())
  gc()
}