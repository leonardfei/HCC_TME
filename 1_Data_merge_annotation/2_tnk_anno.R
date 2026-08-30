#######T_NK#############
T_NK=subset(all_filt2,cl1=="T_NK")

T_NK <- NormalizeData(T_NK, normalization.method = "LogNormalize", scale.factor = 10000)
T_NK <- FindVariableFeatures(T_NK, selection.method = "vst", nfeatures = 2000)
T_NK <- ScaleData(T_NK)

T_NK <- RunPCA (T_NK, features = VariableFeatures(object = T_NK), ndims.print = 1:2)

library(harmony)

T_NK2 = T_NK %>% RunHarmony("Sample", plot_convergence = TRUE)
set.resolutions <- seq(0.5, 1.5, by = 0.1)
pdf(file = "PCA-T_NK_harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

T_NK2 <- FindNeighbors(T_NK2, dims = 1:30,reduction = "harmony")
T_NK2  <- FindClusters(object = T_NK2 , resolution = set.resolutions, verbose = FALSE)
T_NK2  <- RunUMAP(T_NK2 , dims = 1:30,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = T_NK2, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
qsave(T_NK2,"T_NK2_preanno.qs")

FeaturePlot(T_NK2,c("CD4","CD8A"),raster = T)
Idents(T_NK2)=T_NK2$RNA_snn_res.1.2
DotPlot(T_NK2,features = c("CD4","CD8A","MKI67","CD3D","FCGR3A","CD160","ITGA1","FOXP3","TRDC","TRGV9","KIT","IL4I1","SLC4A10"))
DotPlot(T_NK2,features = c("CCR7","SELL","BAG3","CXCL13","FOXP3","GPR183","CD69","PLCG2"))

Idents(T_NK2)=T_NK2$RNA_snn_res.1.2
library(COSG)
marker_anno=cosg(T_NK2,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 30)
gene_anno=marker_anno$name

library(Nebulosa)
plot_density(T_NK2, "PDCD1")

marker=data.frame(cluster=0:33,cell=0:33)
marker[marker$cluster %in% c(9,14),2]<-"NK_FCGR3A"
marker[marker$cluster %in% c(4,13),2]<-"NK_CD160"
marker[marker$cluster %in% c(8,24),2]<-"CD4_Treg_FOXP3"
marker[marker$cluster %in% c(12,33),2]<-"CD4_Tn_LEF1"
marker[marker$cluster %in% c(0),2]<-"CD4_Tcm_CCR7"
marker[marker$cluster %in% c(20),2]<-"CD4_Th1_CXCL13"
marker[marker$cluster %in% c(7),2]<-"CD4_Tcm_FOS"
marker[marker$cluster %in% c(28),2]<-"ILC_KIT"
marker[marker$cluster %in% c(11,29),2]<-"gdT_TRGV9"
marker[marker$cluster %in% c(5,15),2]<-"CD8_NKT_FCGR3A"
marker[marker$cluster %in% c(6,26),2]<-"CD8_MAIT_SLC4A10"
marker[marker$cluster %in% c(10,1),2]<-"CD8_Tex_PDCD1"
marker[marker$cluster %in% c(3),2]<-"CD8_Tstr_HSPA1B"
marker[marker$cluster %in% c(22),2]<-"CD8_Trm_XCL1"
marker[marker$cluster %in% c(31,17,2),2]<-"CD8_Tem_GZMK_APOE"
marker[marker$cluster %in% c(16),2]<-"CD8_Tem_IFIT1"
marker[marker$cluster %in% c(18,27,25),2]<-"CD8_ATP5E"
marker[marker$cluster %in% c(21,23,30,32,19),2]<-"T_NK_cycling"

T_NK2@active.ident=factor(T_NK2@active.ident,levels = 0:33)
T_NK2@meta.data$cl2=sapply(T_NK2@active.ident,function(x){marker[x,2]})
DimPlot(T_NK2,reduction = "umap",group.by = "cl2",label = F,cols = SC_color[10:40],raster = T)
qsave(T_NK2,"TNK2_anno.qs")

#############marker_exp####################
T_NK2=SetIdent(T_NK2,value = "cl2")
library(COSG)
marker_anno=cosg(T_NK2,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 5,expressed_pct = 0.2)
gene_anno=marker_anno$name
library(tidyr)
long_gene_anno<-gather(gene_anno, cell, gene, 1:18)
markers= markers=long_gene_anno$gene
markers=markers[!duplicated(markers)]

library(ggplot2)
DotPlot(T_NK2, features = markers)+
  coord_flip()+
  theme_bw()+ # Use a clean theme and rotate the plot.
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+ # Rotate axis labels by 90 degrees.
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#f6fbef","#b2e1ba","#49add0","#094687"))+ # Set the color gradient.
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) +scale_size(range = c(0,4))

#########cl3############
Idents(T_NK2)=T_NK2$RNA_snn_res.1.2
marker=data.frame(cluster=0:33,cell=0:33)
marker[marker$cluster %in% c(9,14,4,13),2]<-"NK"
marker[marker$cluster %in% c(8,24,12,33,0,20,7),2]<-"CD4_T"
marker[marker$cluster %in% c(28),2]<-"ILC"
marker[marker$cluster %in% c(11,29),2]<-"gdT"
marker[marker$cluster %in% c(5,15,6,26,10,1,3,22,31,17,2,16,18,27,25),2]<-"CD8_T"
marker[marker$cluster %in% c(21,23,30,32,19),2]<-"T_NK_cycling"

T_NK2@active.ident=factor(T_NK2@active.ident,levels = 0:33)
T_NK2@meta.data$cl3=sapply(T_NK2@active.ident,function(x){marker[x,2]})
DimPlot(T_NK2,reduction = "umap",group.by = "cl3",label = F,cols = SC_color[10:40],raster = T)
qsave(T_NK2,"TNK2_anno.qs")
