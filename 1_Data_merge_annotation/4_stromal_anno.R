#######Stromal#############
Stromal=subset(all_filt2,cl1=="Stromal")

Stromal <- NormalizeData(Stromal, normalization.method = "LogNormalize", scale.factor = 10000)
Stromal <- FindVariableFeatures(Stromal, selection.method = "vst", nfeatures = 1500)
Stromal <- ScaleData(Stromal)

Stromal <- RunPCA (Stromal, features = VariableFeatures(object = Stromal), ndims.print = 1:2)

library(harmony)

Stromal2 = Stromal %>% RunHarmony("cohort", plot_convergence = TRUE)
set.resolutions <- seq(0.5, 1.5, by = 0.1)
pdf(file = "PCA-Stromal_harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

Stromal2 <- FindNeighbors(Stromal2, dims = 1:15,reduction = "harmony")
Stromal2  <- FindClusters(object = Stromal2 , resolution = set.resolutions, verbose = FALSE)
Stromal2  <- RunUMAP(Stromal2 , dims = 1:15,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = Stromal2, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
qsave(Stromal2,"Stromal2_preanno.qs")
FeaturePlot(Stromal2,c("FAP","ADIRF","APOC1","HLA-DRB1","RGS5","TOP2A","CFD","CD36"),raster = T,ncol = 4)
library(Nebulosa)
plot_density(Stromal2, "IL34")

Idents(Stromal2)=Stromal2$RNA_snn_res.0.8
library(COSG)
marker_anno=cosg(Stromal2,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 30)
gene_anno=marker_anno$name

#######Endothelial marker
FeaturePlot(Stromal2,c("PROX1","PDPN","ALCAM","CD44","VEGFA"),raster = T) # Analysis note.
FeaturePlot(Stromal2,c("EPHB4","NR2F2","ACKR1","MMRN1","SELP","VCAM1","POSTN","IGFBP7","CCL14","PRCP"),raster = T) # Analysis note.
FeaturePlot(Stromal2,c("GJA4","GJA5","EFNB2","VEGFC","SOX17","DKK2","LTBP4","FBLN5","FN1","MGP","SERPINE2","ENPP2","TFPI","NPR3"),raster = T) # Analysis note.
FeaturePlot(Stromal2,c("CXCL12","CXCR4","ACKR3","LYVE1","DLL4","KCNE3","ESM1","ANGPT2","APLN"),raster = T) # Analysis note.
FeaturePlot(Stromal2,c("CA4","CD36"),raster = T) # Analysis note.
FeaturePlot(Stromal2,c("ISG20","IFIT1","IFIT3"),raster = T) # Analysis note.
FeaturePlot(Stromal2,c("MKI67","TOP2A"),raster = T) # Analysis note.

marker=data.frame(cluster=0:31,cell=0:31)
marker[marker$cluster %in% c(8,29),2]<-"Fb_FAP"
marker[marker$cluster %in% c(25),2]<-"Fb_HLA_DRB1"
marker[marker$cluster %in% c(2),2]<-"Fb_SMC_ADIRF"
marker[marker$cluster %in% c(23,27,31),2]<-"Fb_pericyte_RGS5"
marker[marker$cluster %in% c(6,30),2]<-"Fb_APOC1"
marker[marker$cluster %in% c(5),2]<-"Fb_IL34"
marker[marker$cluster %in% c(15),2]<-"Fb_Epi_KRT19"
marker[marker$cluster %in% c(21),2]<-"EC_lymphatic_PDPN"
marker[marker$cluster %in% c(11,28,13),2]<-"EC_sinusoid_CLEC4G"
marker[marker$cluster %in% c(7),2]<-"EC_vein_ACKR1"
marker[marker$cluster %in% c(4,9),2]<-"EC_artery_GJA4"
marker[marker$cluster %in% c(18),2]<-"EC_T_CD3D"
marker[marker$cluster %in% c(14),2]<-"Stromal_cycling"
marker[marker$cluster %in% c(12),2]<-"EC_vein_PTGS1"
marker[marker$cluster %in% c(1,16,24),2]<-"EC_capillary_CA4"
marker[marker$cluster %in% c(3,22,26,10,19),2]<-"EC_capillary_VWA1"
marker[marker$cluster %in% c(17,0,20),2]<-"EC_APOA1"
marker[marker$cluster %in% c(26),2]<-"EC_Epi_KRT19"

Stromal2@active.ident=factor(Stromal2@active.ident,levels = 0:31)
Stromal2@meta.data$cl2=sapply(Stromal2@active.ident,function(x){marker[x,2]})
DimPlot(Stromal2,reduction = "umap",group.by = "cl2",label = F,cols = SC_color[10:40],raster = T)
qsave(Stromal2,"Stromal2_anno.qs")

#############marker_exp####################
Stromal2=SetIdent(Stromal2,value = "cl2")
library(COSG)
marker_anno=cosg(Stromal2,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 5,expressed_pct = 0.2)
gene_anno=marker_anno$name
library(tidyr)
long_gene_anno<-gather(gene_anno, cell, gene, 1:18)
markers= markers=long_gene_anno$gene
markers=markers[!duplicated(markers)]

library(ggplot2)
DotPlot(Stromal2, features = markers)+
  coord_flip()+
  theme_bw()+ # Use a clean theme and rotate the plot.
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+ # Rotate axis labels by 90 degrees.
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#f6fbef","#b2e1ba","#49add0","#094687"))+ # Set the color gradient.
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) +scale_size(range = c(0,4))



######cl3########
Idents(Stromal2)=Stromal2$RNA_snn_res.0.8

marker=data.frame(cluster=0:31,cell=0:31)
marker[marker$cluster %in% c(8,29,25,2,23,27,31,6,30,5,15),2]<-"Fb"
marker[marker$cluster %in% c(21,11,28,13,7,4,9,18,12,1,16,24,3,22,26,10,19,17,0,20,26),2]<-"EC"
marker[marker$cluster %in% c(14),2]<-"Stromal_cycling"

Stromal2@active.ident=factor(Stromal2@active.ident,levels = 0:31)
Stromal2@meta.data$cl3=sapply(Stromal2@active.ident,function(x){marker[x,2]})
DimPlot(Stromal2,reduction = "umap",group.by = "cl3",label = F,cols = SC_color[10:40],raster = T)
qsave(Stromal2,"Stromal2_anno.qs")
