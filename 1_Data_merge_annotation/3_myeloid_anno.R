#######Myeloid#############
Myeloid=subset(all_filt2,cl1=="Myeloid")

Myeloid <- NormalizeData(Myeloid, normalization.method = "LogNormalize", scale.factor = 10000)
Myeloid <- FindVariableFeatures(Myeloid, selection.method = "vst", nfeatures = 1500)
Myeloid <- ScaleData(Myeloid)

Myeloid <- RunPCA (Myeloid, features = VariableFeatures(object = Myeloid), ndims.print = 1:2)

library(harmony)

Myeloid2 = Myeloid %>% RunHarmony("cohort", plot_convergence = TRUE)
set.resolutions <- seq(0.5, 1.5, by = 0.1)
pdf(file = "PCA-Myeloid_harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

Myeloid2 <- FindNeighbors(Myeloid2, dims = 1:10,reduction = "harmony")
Myeloid2  <- FindClusters(object = Myeloid2 , resolution = set.resolutions, verbose = FALSE)
Myeloid2  <- RunUMAP(Myeloid2 , dims = 1:10,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = Myeloid2, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
qsave(Myeloid2,"myeloid2_preanno.qs")
FeaturePlot(Myeloid2,c("FCGR3B","FCN1","CD68","CD1C","LAMP3","CLEC9A","TPSAB1","CLEC4C","CLC","MKI67"),raster = T,ncol = 4)

Idents(Myeloid2)=Myeloid2$RNA_snn_res.0.7
library(COSG)
marker_anno=cosg(Myeloid2,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 30)
gene_anno=marker_anno$name

marker=data.frame(cluster=0:25,cell=0:25)
marker[marker$cluster %in% c(9),2]<-"Neutrophil"
marker[marker$cluster %in% c(4,15),2]<-"DC_CD1C"
marker[marker$cluster %in% c(14),2]<-"DC_CLEC9A"
marker[marker$cluster %in% c(22),2]<-"DC_LAMP3"
marker[marker$cluster %in% c(16),2]<-"pDC_CLEC4C"
marker[marker$cluster %in% c(20),2]<-"Mast_TPSAB1"
marker[marker$cluster %in% c(25),2]<-"Eosinophil"
marker[marker$cluster %in% c(3,6,11,17,24),2]<-"Mono_like_Mph"
marker[marker$cluster %in% c(5,21),2]<-"Mph_MARCO"
marker[marker$cluster %in% c(0,2,8,7,18),2]<-"Mph_IGF1"
marker[marker$cluster %in% c(1,19),2]<-"Mph_SPP1"
marker[marker$cluster %in% c(12),2]<-"Mph_STARD13"
marker[marker$cluster %in% c(13,23),2]<-"Mph_CD68_CD247"
marker[marker$cluster %in% c(10),2]<-"Myeloid_cycling"

Myeloid2@active.ident=factor(Myeloid2@active.ident,levels = 0:25)
Myeloid2@meta.data$cl2=sapply(Myeloid2@active.ident,function(x){marker[x,2]})
DimPlot(Myeloid2,reduction = "umap",group.by = "cl2",label = F,cols = SC_color[10:30],raster = T)
qsave(Myeloid2,"myeloid2_anno.qs")

#############marker_exp####################
Myeloid2=SetIdent(Myeloid2,value = "cl2")
library(COSG)
marker_anno=cosg(Myeloid2,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 5,expressed_pct = 0.2)
gene_anno=marker_anno$name
library(tidyr)
long_gene_anno<-gather(gene_anno, cell, gene, 1:14)
markers= markers=long_gene_anno$gene
markers=markers[!duplicated(markers)]

library(ggplot2)
DotPlot(Myeloid2, features = markers)+
  coord_flip()+
  theme_bw()+ # Use a clean theme and rotate the plot.
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+ # Rotate axis labels by 90 degrees.
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#f6fbef","#b2e1ba","#49add0","#094687"))+ # Set the color gradient.
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) +scale_size(range = c(0,4))


#######cl3#########
Idents(Myeloid2)=Myeloid2$RNA_snn_res.0.7

marker=data.frame(cluster=0:25,cell=0:25)
marker[marker$cluster %in% c(9),2]<-"Neutrophil"
marker[marker$cluster %in% c(4,15,14,22,16),2]<-"DC"
marker[marker$cluster %in% c(20),2]<-"Mast"
marker[marker$cluster %in% c(25),2]<-"Eosinophil"
marker[marker$cluster %in% c(3,6,11,17,24,5,21,0,2,8,7,18,1,19,12,13,23),2]<-"Mph"
marker[marker$cluster %in% c(10),2]<-"Myeloid_cycling"

Myeloid2@active.ident=factor(Myeloid2@active.ident,levels = 0:25)
Myeloid2@meta.data$cl3=sapply(Myeloid2@active.ident,function(x){marker[x,2]})
DimPlot(Myeloid2,reduction = "umap",group.by = "cl3",label = F,cols = SC_color[10:30],raster = T)
qsave(Myeloid2,"myeloid2_anno.qs")
