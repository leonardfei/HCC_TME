# Use HCC_TME_DATA_DIR for external inputs; no local working directory is embedded.
data_dir <- Sys.getenv("HCC_TME_DATA_DIR", unset = "data")

library(qs)
library(Seurat)
all=qread(file.path(data_dir, "allcell_anno.qs"))

C1_xue=all
C1_xue=subset(C1_xue,Cancer_type=="PB",invert=T)
Idents(C1_xue)=C1_xue$celltype
C1_xue$cohort=rep("c0",969485)

C2_hbv=qread(file.path(data_dir, "orign_tissue.qs"))
C2_hbv$cohort=rep("c5",106592)
C2_hbv$Cancer_type=C2_hbv$hbv_state1
C2_hbv$Sample=C2_hbv$orig.ident

C3_healthy=qread(file.path(data_dir, "C3_seurat.qs"))
C3_healthy$cohort=rep("c3",8419)
C3_healthy$Cancer_type=rep("Healthy",8419)

C4_inte=qread(file.path(data_dir, "allcell_merge.qs"))
C4_inte$Cancer_type=C4_inte$sampletype
C4_inte@meta.data <- C4_inte@meta.data[, !(names(C4_inte@meta.data) %in% c("RNA_snn_res.0.1", "RNA_snn_res.0.2","RNA_snn_res.0.3","RNA_snn_res.0.4","RNA_snn_res.0.5","RNA_snn_res.0.6","RNA_snn_res.0.7","RNA_snn_res.0.8","RNA_snn_res.0.9","RNA_snn_res.1","RNA_snn_res.1.1","RNA_snn_res.1.2","RNA_snn_res.1.3","RNA_snn_res.1.4","RNA_snn_res.1.5","RNA_snn_res.1.6","RNA_snn_res.1.7","RNA_snn_res.1.8","RNA_snn_res.1.9","RNA_snn_res.2"))]


all=merge(C1_xue,C2_hbv)
all=merge(all,C3_healthy)
all=merge(all,C4_inte)

library(sceasy)
time.R2py = system.time({
  sceasy::convertFormat(all, from="seurat", to="anndata",
                        outFile='./all.h5ad')
})
print(time.R2py)

time.R2py = system.time({
  sceasy::convertFormat(C4_inte, from="seurat", to="anndata",
                        outFile='./C4_inte.h5ad')
})
print(time.R2py)

# Run `python 1_Data_merge_annotation/merge_anndata.py` before importing
# `all_merge.h5ad` back into R.


library(Seurat)
library(dior)
# Analysis note.
time.py2R = system.time({
  seurat.data <- read_h5ad(file = './all_merge.h5ad',
                           assay_name = 'RNA',
                           target.object = 'seurat')
})
print(time.py2R)

# Merge the data.
# Use HCC_TME_DATA_DIR for external inputs; no local working directory is embedded.
Idents(all)=all$Cancer_type
all_hcc=subset(all,idents = c("Healthy","HBV","AL","HCC"))
C4_inte=qread("C4_inte.qs")
C4_inte_hcc=subset(C4_inte,Cancer_type=="ICC",invert=T)
all_merge_hcc=merge(all_hcc,C4_inte_hcc)
qsave(all_merge_hcc,"all_merge_hcc_orig.qs")

all_merge_hcc[["percent.mt"]] <- PercentageFeatureSet(all_merge_hcc, pattern = "^MT-")
all_filt <- subset(all_merge_hcc, subset = nFeature_RNA > 100 & percent.mt < 50 & nFeature_RNA < 6000 & nCount_RNA < 30000)
all_filt$Cancer_type=factor(all_filt$Cancer_type,levels = c("Healthy","HBV","AL","HCC"))


all_filt <- NormalizeData(all_filt, normalization.method = "LogNormalize", scale.factor = 10000)
all_filt <- FindVariableFeatures(all_filt, selection.method = "vst", nfeatures = 3000)
all_filt <- ScaleData(all_filt)

all_filt <- RunPCA (all_filt, features = VariableFeatures(object = all_filt), ndims.print = 1:2)

library(harmony)

all_filt2 = all_filt %>% RunHarmony("cohort", plot_convergence = TRUE)
set.resolutions <- seq(0.1, 1, by = 0.2)
pdf(file = "PCA-harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

all_filt2 <- FindNeighbors(all_filt2, dims = 1:50,reduction = "harmony")
all_filt2  <- FindClusters(object = all_filt2 , resolution = set.resolutions, verbose = FALSE)
all_filt2  <- RunUMAP(all_filt2 , dims = 1:50,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = all_filt2, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
#plan("multisession", workers = 1)
qsave(all_filt2,"all_filt2_preanno.qs")


library(paletteer)
color_name=palettes_d_names
pal=paletteer_d("khroma::smoothrainbow")

FeaturePlot(all_filt2,c("ACTB","PTPRC","ALB","AFP","EPCAM",
                        "APOA2","COL1A1","CDH5","CD3D","KLRF1",
                        "GNLY","MS4A1","MZB1","CD68","CD1C",
                        "FCGR3B","CLEC4C","TPSAB1"),cols = pal,raster = T,ncol = 5)


Idents(all_filt2)=all_filt2$RNA_snn_res.0.1
DotPlot(all_filt2,features = "CD3D")
marker=data.frame(cluster=0:25,cell=0:25)
marker[marker$cluster %in% c(9,5,18,0,3,6,11,25),2]<-"T_NK"
marker[marker$cluster %in% c(8,12),2]<-"B_PCs"
marker[marker$cluster %in% c(1,14,16,20,22),2]<-"Myeloid"
marker[marker$cluster %in% c(7,4),2]<-"Stromal"
marker[marker$cluster %in% c(24,17,21,13,19,2,15,23,10),2]<-"Malignant_Epi"

all_filt2@active.ident=factor(all_filt2@active.ident,levels = 0:25)
all_filt2@meta.data$cl1=sapply(all_filt2@active.ident,function(x){marker[x,2]})
DimPlot(all_filt2,reduction = "umap",group.by = "cl1",label = F,cols = SC_color[10:30],raster = T)
qsave(all_filt2,"all_filt2_cl1.qs")
rm(all_filt)
