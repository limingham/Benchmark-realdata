library(Seurat)
library(fastSave)


bulk_data <- readRDS("~/Benchmark-realdata-main/fig2/data_preprocess/bulk_data.rds")
sc_TNBC = readRDS.lbzip2('~/Benchmark-realdata-main/fig2/data_preprocess/seurat_use.rdsFS')
sc_TNBC = subset(sc_TNBC,subtype == 'ER+')

inter_name = intersect(colnames(bulk_data),unique(sc_TNBC@meta.data$orig.ident))

for(sample in inter_name){
  
  
  TNBC_bulk = as.data.frame(bulk_data[,sample])
  TNBC_bulk[,2] = TNBC_bulk[,1]
  row.names(TNBC_bulk) = row.names(bulk_data)
  seurat_TNBC = subset(sc_TNBC,orig.ident != sample)

  results_path <- paste0("~/Benchmark-realdata-main/fig2/deconvolution_results/ER_realbulk/",sample)
  if (dir.exists(results_path)) {
  } else {
    dir.create(results_path)
  }
  source("/home/lmh/ReCIDE_git/ReCIDE/R/ReCIDE_all_function_beta.R")
  # source("~/revision/base_code/dwls_function.R")
  
  # 
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace_all(seurat_TNBC@meta.data[["celltype_use"]], " ", "_")
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace_all(seurat_TNBC@meta.data[["celltype_use"]], "\\/", "_")
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace_all(seurat_TNBC@meta.data[["celltype_use"]], "\\.", "_")
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace_all(seurat_TNBC@meta.data[["celltype_use"]], "\\(", "")
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace_all(seurat_TNBC@meta.data[["celltype_use"]], "\\)", "")
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace_all(seurat_TNBC@meta.data[["celltype_use"]], ",", "_")
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace(seurat_TNBC@meta.data[["celltype_use"]], "-", "")
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace(seurat_TNBC@meta.data[["celltype_use"]], "\\+", "")
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace(seurat_TNBC@meta.data[["celltype_use"]], "_III", "3")
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace(seurat_TNBC@meta.data[["celltype_use"]], "_II", "2")
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace(seurat_TNBC@meta.data[["celltype_use"]], "_I", "1")
  seurat_TNBC@meta.data[["celltype_use"]] <- str_replace_all(seurat_TNBC@meta.data[["celltype_use"]], "\\-", "_")
  
  
  
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
  
  source('~/ReCIDE/benchmark_syq/main_run_all_revise.R')
  func_run_ReCIDE(SC_ref = seurat_TNBC,
                  EXP_df = TNBC_bulk,
                  celltype = celltype1,
                  subject = subject1,
                  dir_ref = dir_ref_recide,
                  dir_results = dir_results_recide)
  
  
  ######run_bisque
  source('~/ReCIDE/benchmark_syq/main_run_all_revise.R')
  
  func_run_bisque(ref_seurat = seurat_TNBC,
                  bulkdata = bulk_data,
                  celltype = celltype1,
                  subject = subject1,
                  dir_results = dir_results_bisque)
  gc()
  
  ######run_music
  source('~/ReCIDE/benchmark_syq/main_run_all_revise.R')
  
  func_run_music(ref_seurat = seurat_TNBC,
                 EXP_df = TNBC_bulk,
                 celltype = celltype1,
                 subject = subject1,
                 dir_results = dir_results_music)
  gc()
  
  
  
  
  ######run_CIBERSORT
  source('~/ReCIDE/benchmark_syq/main_run_all_revise.R')
  
  func_run_CIBERSORT(combined.data = seurat_TNBC,
                     bulk.mtx = TNBC_bulk,
                     celltype = celltype1,
                     subject = subject1,
                     dir_ref = dir_ref_CIBERSORT,
                     dir_results = dir_results_CIBERSORT)
  gc()
  
  
  ######run_DWLS
  source('~/ReCIDE/benchmark_syq/main_run_all_revise.R')
  
  func_run_DWLS(scdata_test = seurat_TNBC,
                bulkdata = TNBC_bulk,
                celltype = celltype1,
                dir_DWLS_ref = dir_ref_DWLS,
                dir_DWLS_results = dir_results_DWLS)
  gc()
  
  
  ######run_bayes
  source('~/ReCIDE/benchmark_syq/main_run_all_revise.R')
  
  func_run_bayes(ref_seurat = seurat_TNBC,
                 bulk.mtx = TNBC_bulk,
                 celltype = celltype1,
                 subject = subject1,
                 dir_results = dir_results_bayes,
                 key1 = 'Epithelial')
  
  # rm(list = ls())
  gc()
}