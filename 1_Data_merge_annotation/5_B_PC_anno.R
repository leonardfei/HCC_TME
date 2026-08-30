#######B_PC#############
all_filt2=qread("all_filt2_cl1.qs")
B_PC=subset(all_filt2,cl1=="B_PCs")

B_PC <- NormalizeData(B_PC, normalization.method = "LogNormalize", scale.factor = 10000)
B_PC <- FindVariableFeatures(B_PC, selection.method = "vst", nfeatures = 1000)
B_PC <- ScaleData(B_PC)

B_PC <- RunPCA (B_PC, features = VariableFeatures(object = B_PC), ndims.print = 1:2)

library(harmony)

B_PC2 = B_PC %>% RunHarmony("cohort", plot_convergence = TRUE)
set.resolutions <- seq(0.5, 1.5, by = 0.1)
pdf(file = "PCA-B_PC_harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

B_PC2 <- FindNeighbors(B_PC2, dims = 1:10,reduction = "harmony")
B_PC2  <- FindClusters(object = B_PC2 , resolution = set.resolutions, verbose = FALSE)
B_PC2  <- RunUMAP(B_PC2 , dims = 1:10,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = B_PC2, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
qsave(B_PC2,"B_PC2_preanno.qs")
FeaturePlot(B_PC2,c("CD79A","CD3D","MS4A1","MZB1","MKI67"),raster = T,ncol = 4)

Idents(B_PC2)=B_PC2$RNA_snn_res.0.6
library(COSG)
marker_anno=cosg(B_PC2,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 30)
gene_anno=marker_anno$name

marker=data.frame(cluster=0:16,cell=0:16)
marker[marker$cluster %in% c(1),2]<-"B_naive_TCL1A"
marker[marker$cluster %in% c(5,12),2]<-"B_memory_APOA1"
marker[marker$cluster %in% c(2),2]<-"B_memory_BAG3"
marker[marker$cluster %in% c(0,16),2]<-"B_memory_GPR183"
marker[marker$cluster %in% c(4),2]<-"B_ABC_ITGAX"
marker[marker$cluster %in% c(10),2]<-"B_PC_cycling"
marker[marker$cluster %in% c(3,9,11,15,7,6),2]<-"PC_MZB1"
marker[marker$cluster %in% c(13,8,14),2]<-"B_PC_filt"

B_PC2@active.ident=factor(B_PC2@active.ident,levels = 0:16)
B_PC2@meta.data$cl2=sapply(B_PC2@active.ident,function(x){marker[x,2]})
DimPlot(B_PC2,reduction = "umap",group.by = "cl2",label = F,cols = SC_color[10:30],raster = T)
B_PC2_filt=subset(B_PC2,cl2=="B_PC_filt",invert=T)

qsave(B_PC2_filt,"B_PC2_filt_anno.qs")

#############marker_exp####################
B_PC2_filt=SetIdent(B_PC2_filt,value = "cl2")
library(COSG)
marker_anno=cosg(B_PC2_filt,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 5,expressed_pct = 0.2)
gene_anno=marker_anno$name
library(tidyr)
long_gene_anno<-gather(gene_anno, cell, gene, 1:7)
markers= markers=long_gene_anno$gene
markers=markers[!duplicated(markers)]
markers[26]="ITGAX"
markers[30]="TBX21"

library(ggplot2)
DotPlot(B_PC2_filt, features = markers)+
  coord_flip()+
  theme_bw()+ # Use a clean theme and rotate the plot.
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+ # Rotate axis labels by 90 degrees.
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#f6fbef","#b2e1ba","#49add0","#094687"))+ # Set the color gradient.
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) +scale_size(range = c(0,4))

##########cl3#########
Idents(B_PC2)=B_PC2$RNA_snn_res.0.6
marker=data.frame(cluster=0:16,cell=0:16)
marker[marker$cluster %in% c(1,5,12,2,0,16,4),2]<-"B_cell"
marker[marker$cluster %in% c(10),2]<-"B_PC_cycling"
marker[marker$cluster %in% c(3,9,11,15,7,6),2]<-"PC"
marker[marker$cluster %in% c(13,8,14),2]<-"B_PC_filt"

B_PC2@active.ident=factor(B_PC2@active.ident,levels = 0:16)
B_PC2@meta.data$cl3=sapply(B_PC2@active.ident,function(x){marker[x,2]})
DimPlot(B_PC2,reduction = "umap",group.by = "cl3",label = F,cols = SC_color[10:30],raster = T)
B_PC2_filt=subset(B_PC2,cl2=="B_PC_filt",invert=T)

qsave(B_PC2_filt,"B_PC2_filt_anno.qs")
