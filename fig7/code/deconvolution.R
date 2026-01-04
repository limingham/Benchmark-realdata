

library(Seurat)
library(fastSave)

results_path <- "~/Benchmark-realdata-main/fig7/results/ESCA/"
if (dir.exists(results_path)) {
} else {
  dir.create(results_path)
}
source("~/Benchmark-realdata-main/fig7/code/ReCIDE_all_function_beta.R")



seurat_all = readRDS.lbzip2('~/Benchmark-realdata-main/fig7/reference/ESCA.rdsFS',n.cores = 100)
seurat_all@meta.data[["anno_pancancer"]] = str_replace_all(seurat_all@meta.data[["anno_pancancer"]], " ", "_")
seurat_all@meta.data[["anno_pancancer"]] = str_replace_all(seurat_all@meta.data[["anno_pancancer"]], "\\+", "") 
seurat_all@meta.data[["anno_pancancer"]] = str_replace_all(seurat_all@meta.data[["anno_pancancer"]], "-", "_") 
cell_table = as.data.frame(table(seurat_all@meta.data[["anno_pancancer"]]))
cell_table = cell_table[cell_table[,2]>3,]
seurat_ER = subset(seurat_all,anno_pancancer %in% cell_table[,1])

seurat_ER = subset(seurat_ER,Tissue %in% "Tumor")


ER_bulk <- readRDS("~/Benchmark-realdata-main/fig7/bulk_data/ESCA.rds")
ER_bulk[] <- lapply(ER_bulk, function(x) as.numeric(as.character(x)))


celltype1 = "anno_pancancer"
subject1 = "Patient"

dir_ref_recide = paste0(results_path,'/ref_ReCIDE.rds')
dir_results_recide = paste0(results_path,'/results_ReCIDE.rds')
dir_results_bayes = paste0(results_path,'/results_bayes.rds')



source('~/Benchmark-realdata-main/fig7/code/run_all_methods.R')
func_run_ReCIDE(SC_ref = seurat_ER,
                EXP_df = ER_bulk,
                celltype = celltype1,
                subject = subject1,
                dir_ref = dir_ref_recide,
                dir_results = dir_results_recide)



