data_dir <- Sys.getenv("HCC_TME_DATA_DIR", unset = "data")

library(qs)
library(Seurat)
# Use HCC_TME_DATA_DIR for external inputs; no local working directory is embedded.
mph=qread(file.path(data_dir, "Mo_Mph.qs"))
Idents(mph)=mph$Cancer_type
mph_filt=subset(mph,idents = c("AL","HCC","PB"))

mph_filt <- NormalizeData(mph_filt, normalization.method = "LogNormalize", scale.factor = 10000)
mph_filt <- FindVariableFeatures(mph_filt, selection.method = "vst", nfeatures = 1500)
mph_filt <- ScaleData(mph_filt)

mph_filt <- RunPCA (mph_filt, features = VariableFeatures(object = mph_filt), ndims.print = 1:2)
ElbowPlot(object = mph_filt, ndims = 100)

library(harmony)
mph_filt2 = mph_filt %>% RunHarmony("Sample", plot_convergence = TRUE)

set.resolutions <- seq(0.2, 1, by = 0.1)
pdf(file = "PCA-mph_filt2.pdf")
mph_filt2 <- FindNeighbors(mph_filt2, dims = 1:20,reduction = "harmony")
mph_filt2  <- FindClusters(object = mph_filt2 , resolution = set.resolutions, verbose = FALSE)
mph_filt2  <- RunUMAP(mph_filt2 , dims = 1:20,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = mph_filt2, reduction = 'umap',raster =T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()

DimPlot(mph_filt2,group.by = "Cancer_type",raster = T,cols = c("#d0ddef","#6d75a5","#bb605f"))
library(Nebulosa)
plot_density(mph_filt2, "IL1B",size = 0.1,reduction = "umap")
A=mph_filt2$new_clusters
A=as.data.frame(A)
A$type=A$A
A[A$A=="Mph_SLC40A1",]$type="Mph_IGF1"
mph_filt2$new_clusters=A$type
qsave(mph_filt2,"mph_filt2.qs")
datafilt <- mph_filt2[,sample(colnames(mph_filt2), 20000)]
#########VIA##############
seurat_object <- DietSeurat(
  object = mph_filt2,
  counts = TRUE, # Analysis note.
  data = TRUE, # Analysis note.
  scale.data = FALSE, # Analysis note.
  features = rownames(mph_filt2), # Analysis note.
  assays = "RNA", # Analysis note.
  dimreducs = c("pca","umap","harmony"), # Perform dimensionality reduction.
  graphs = c("RNA_nn", "RNA_snn"), # Analysis note.
  misc = TRUE # Analysis note.
)

seurat_object <- DietSeurat(
  object = datafilt,
  counts = TRUE, # Analysis note.
  data = TRUE, # Analysis note.
  scale.data = FALSE, # Analysis note.
  features = rownames(datafilt), # Analysis note.
  assays = "RNA", # Analysis note.
  dimreducs = c("pca","umap","harmony"), # Perform dimensionality reduction.
  graphs = c("RNA_nn", "RNA_snn"), # Analysis note.
  misc = TRUE # Analysis note.
)

i <- sapply(seurat_object@meta.data, is.factor)
seurat_object@meta.data[i] <- lapply(seurat_object@meta.data[i], as.character)

library(SeuratDisk)
SaveH5Seurat(seurat_object, filename = "mph_pb_downsample.h5Seurat", overwrite = TRUE)
Convert("mph_pb_downsample.h5Seurat", dest = "h5ad", assay = "RNA", overwrite = TRUE)
# Run `python 3_mph_analysis/mph_pb_trajectory.py` to perform the pyVIA
# trajectory analysis and generate the trajectory figures.
