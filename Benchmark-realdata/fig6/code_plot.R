library(Seurat)
library(fastSave)
library(harmony)



B_cell = readRDS.lbzip2('~/Benchmark-realdata-main/fig6/Seurat_B_anno.rdsFS',n.cores = 50)


B_cell_plot = DimPlot(B_cell,group.by = 'anno_pancancer',raster=FALSE,reduction = 'umap')



Endo = readRDS.lbzip2('~/Benchmark-realdata-main/fig6/Seurat_Endo_anno.rdsFS',n.cores = 40)

Endo_plot = DimPlot(Endo,reduction = "umap",group.by = 'anno_pancancer',raster=FALSE,label = TRUE)


Fibroblast = readRDS.lbzip2('~/Benchmark-realdata-main/fig6/Seurat_F_anno.rdsFS',n.cores = 40)

Fibroblast_plot = DimPlot(Fibroblast,reduction = "umap",group.by = 'anno_pancancer',raster=FALSE,label = TRUE)


myeloid = readRDS.lbzip2('~/Benchmark-realdata-main/fig6/Seurat_M_anno.rdsFS',n.cores = 50)

myeloid_plot = DimPlot(myeloid,group.by = 'anno_pancancer',raster=FALSE,reduction = 'umap')


T_cell = readRDS.lbzip2('~/Benchmark-realdata-main/fig6/Seurat_T_anno.rdsFS',n.cores = 50)

T_cell_plot = DimPlot(T_cell,group.by = 'anno_pancancer',raster=FALSE,reduction = 'umap')
