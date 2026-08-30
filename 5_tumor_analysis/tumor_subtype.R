data_dir <- Sys.getenv("HCC_TME_DATA_DIR", unset = "data")

rm(list=ls())
# Use HCC_TME_DATA_DIR for external inputs; no local working directory is embedded.
library(qs)
library(Seurat)
tumor=qread(file.path(data_dir, "tumor_epi.qs"))

tumor <- NormalizeData(tumor, normalization.method = "LogNormalize", scale.factor = 10000)
tumor <- FindVariableFeatures(tumor, selection.method = "vst", nfeatures = 2000)
tumor <- ScaleData(tumor)
tumor <- RunPCA (tumor, features = VariableFeatures(object = tumor), ndims.print = 1:2)
ElbowPlot(tumor,ndims = 50)

library(harmony)
tumor2 = tumor %>% RunHarmony("Sample", plot_convergence = TRUE)
set.resolutions <- seq(0.2, 1, by = 0.1)
pdf(file = "PCA-tumor_harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

tumor2 <- FindNeighbors(tumor2, dims = 1:15,reduction = "harmony")
tumor2  <- FindClusters(object = tumor2 , resolution = set.resolutions, verbose = FALSE)
tumor2  <- RunUMAP(tumor2 , dims = 1:15,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = tumor2, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()


Idents(tumor2)=tumor2$RNA_snn_res.0.2
tumor3=subset(tumor2,idents=c("2","8"),invert=T)


tumor3 <- NormalizeData(tumor3, normalization.method = "LogNormalize", scale.factor = 10000)
tumor3 <- FindVariableFeatures(tumor3, selection.method = "vst", nfeatures = 2000)
tumor3 <- ScaleData(tumor3)
tumor3 <- RunPCA (tumor3, features = VariableFeatures(object = tumor3), ndims.print = 1:2)
ElbowPlot(tumor3,ndims = 50)

library(harmony)
tumor3 = tumor3 %>% RunHarmony("Sample", plot_convergence = TRUE)
set.resolutions <- seq(0.2, 1, by = 0.1)
pdf(file = "PCA-tumor3_harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

tumor3 <- FindNeighbors(tumor3, dims = 1:15,reduction = "harmony")
tumor3  <- FindClusters(object = tumor3 , resolution = set.resolutions, verbose = FALSE)
tumor3  <- RunUMAP(tumor3 , dims = 1:15,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = tumor3, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
Idents(tumor3)=tumor3$RNA_snn_res.0.3
library(COSG)
marker_anno=cosg(tumor3,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 30)
gene_anno=marker_anno$name

tumor4=subset(tumor3,idents=c("9","7","11"),invert=T)


tumor4 <- NormalizeData(tumor4, normalization.method = "LogNormalize", scale.factor = 10000)
tumor4 <- FindVariableFeatures(tumor4, selection.method = "vst", nfeatures = 2000)
tumor4 <- ScaleData(tumor4)
tumor4 <- RunPCA (tumor4, features = VariableFeatures(object = tumor4), ndims.print = 1:2)
ElbowPlot(tumor4,ndims = 50)

library(harmony)
tumor4 = tumor4 %>% RunHarmony("Sample", plot_convergence = TRUE)
set.resolutions <- seq(0.2, 1, by = 0.1)
pdf(file = "PCA-tumor4_harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

tumor4 <- FindNeighbors(tumor4, dims = 1:15,reduction = "harmony")
tumor4  <- FindClusters(object = tumor4 , resolution = set.resolutions, verbose = FALSE)
tumor4  <- RunUMAP(tumor4 , dims = 1:15,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = tumor4, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()

Idents(tumor4)=tumor4$RNA_snn_res.0.2
library(COSG)
marker_anno=cosg(tumor4,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 30)
gene_anno=marker_anno$name
tumor5=subset(tumor4,idents="6",invert=T)


marker=data.frame(cluster=0:5,cell=0:5)
marker[marker$cluster %in% c(0),2]<-"Ca_CFHR1"
marker[marker$cluster %in% c(1),2]<-"Ca_BICC1"
marker[marker$cluster %in% c(2),2]<-"Ca_SEPP1"
marker[marker$cluster %in% c(3),2]<-"Ca_EGR1"
marker[marker$cluster %in% c(4),2]<-"Ca_TOP2A"
marker[marker$cluster %in% c(5),2]<-"Epi_KRT19"

tumor5@active.ident=factor(tumor5@active.ident,levels = 0:5)
tumor5@meta.data$cl2=sapply(tumor5@active.ident,function(x){marker[x,2]})
DimPlot(tumor5,reduction = "umap",group.by = "cl2",label = F,cols = SC_color[10:40],raster = T)
qsave(tumor5,"tumor5_anno.qs")

library(Nebulosa)
plot_density(tumor5, "")

tumor_hcc=subset(tumor5,Cancer_type=="HCC")
A=as.data.frame(table(tumor_hcc$Sample))
module=read.csv(file.path(data_dir, "cluster4.csv"),header = T)
module=module[match(A$Var1,module$X),]

tumor_hcc@active.ident=factor(tumor_hcc@active.ident,levels = c(module$X))
tumor_hcc@meta.data$module=sapply(tumor_hcc@active.ident,function(x){module[x,2]})
DimPlot(tumor_hcc,reduction = "umap",group.by = "module",label = F,cols = SC_color[5:40],raster = T)



#############cellratio_module_ca#########
library(ggplot2)
library(dplyr)
library(ggalluvial)

Ratio <- tumor_hcc@meta.data %>%
  group_by(module, cl2) %>% # Group the observations.
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))

ggplot(Ratio, aes(x =module, y= relative_freq, fill = cl2,
                  stratum=cl2, alluvium=cl2)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.3, knot.pos=0.2)+ # Set knot.pos to control the curvature of the alluvial connections.
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  scale_fill_manual(values = SC_color[20:30])









# Calculate mean gene-set expression.
geneset=read.csv(file.path(data_dir, "cancersea_features.csv"),header = T)
#geneset=as.list(geneset)
#tumor_hcc_hbv=subset(tumor_hcc_hbv,cl4=="Epi_KRT19",invert=T)
aver_dt<- AverageExpression(tumor_hcc,
                            #features = gene,
                            group.by = 'cl2',
                            slot= 'data')
aver_dt<- as.data.frame(aver_dt$RNA)

C=1:6
C=as.data.frame(C)
C$id=1:6

for (i in 1:length(colnames(geneset))) {
  g1=aver_dt[geneset[,i],]
  g1=na.omit(g1)
  out1=apply(g1,2,mean)
  C=cbind(C,out1)
}
C=C[-(1:2)]
rownames(C)=colnames(aver_dt)
colnames(C)=names(geneset)
n=as.data.frame(t(C))

library(pheatmap)
color=colorRampPalette(c('#035C8C','white','#BD2131'))(100)
n=t(scale(t(n)))
n[n>2]=2
n[n<(-2)]=-2

pheatmap(n,color = color,cluster_rows = F,cluster_cols = F,border_color = NA)


aver_dt<- AverageExpression(tumor_hcc,
                            #features = gene,
                            group.by = 'module',
                            slot= 'data')
aver_dt<- as.data.frame(aver_dt$RNA)

C=1:4
C=as.data.frame(C)
C$id=1:4

for (i in 1:length(colnames(geneset))) {
  g1=aver_dt[geneset[,i],]
  g1=na.omit(g1)
  out1=apply(g1,2,mean)
  C=cbind(C,out1)
}
C=C[-(1:2)]
rownames(C)=colnames(aver_dt)
colnames(C)=names(geneset)
n=as.data.frame(t(C))

library(pheatmap)
color=colorRampPalette(c('#035C8C','white','#BD2131'))(100)
n=t(scale(t(n)))
n[n>2]=2
n[n<(-2)]=-2

pheatmap(n,color = color,cluster_rows = F,cluster_cols = F,border_color = NA)

###################
tumor_hcc=qread(file.path(data_dir, "tumor_hcc_infercnv.qs"))

tumor_hcc_noepi=subset(tumor_hcc,cl2=="Epi_KRT19",invert=T)


library(paletteer)
color_name=palettes_d_names
pal=paletteer_d("khroma::smoothrainbow")
FeaturePlot(epi,"cnvscore",raster = F,cols = pal)
DimPlot(tumor_hcc_noepi,group.by = "module",cols = c("#18a6b1","#deb316","#c85a52","#238bbc"),raster = T)

library(ggplot2)
DotPlot(tumor_hcc_noepi, features = "cnvscore",group.by = "cl2")+
  #coord_flip()+
  theme_bw()+ # Use a clean theme and rotate the plot.
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+ # Rotate axis labels by 90 degrees.
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FEE08B", "#FDAE61","#F46D43",
                                                          "#D53E4F", "#9E0142"))+ # Set the color gradient.
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
# Calculate mean gene-set expression.
geneset=read.csv(file.path(data_dir, "cancersea_features.csv"),header = T)
#geneset=as.list(geneset)
#tumor_hcc_hbv=subset(tumor_hcc_hbv,cl4=="Epi_KRT19",invert=T)
aver_dt<- AverageExpression(tumor_hcc_noepi,
                            #features = gene,
                            group.by = 'cl2',
                            slot= 'data')
aver_dt<- as.data.frame(aver_dt$RNA)

C=1:5
C=as.data.frame(C)
C$id=1:5

for (i in 1:length(colnames(geneset))) {
  g1=aver_dt[geneset[,i],]
  g1=na.omit(g1)
  out1=apply(g1,2,mean)
  C=cbind(C,out1)
}
C=C[-(1:2)]
rownames(C)=colnames(aver_dt)
colnames(C)=names(geneset)
n=as.data.frame(t(C))

library(pheatmap)
color=colorRampPalette(c('#035C8C','white','#BD2131'))(100)
n=t(scale(t(n)))
n[n>2]=2
n[n<(-2)]=-2

pheatmap(n,color = color,cluster_rows = F,cluster_cols = F,border_color = NA)

#############cellratio_module_ca_noepi#########
library(ggplot2)
library(dplyr)
library(ggalluvial)

Ratio <- tumor_hcc_noepi@meta.data %>%
  group_by(module, cl2) %>% # Group the observations.
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))

ggplot(Ratio, aes(x =module, y= relative_freq, fill = cl2,
                  stratum=cl2, alluvium=cl2)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.3, knot.pos=0.2)+ # Set knot.pos to control the curvature of the alluvial connections.
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  scale_fill_manual(values = SC_color[20:30])









##############marker_exp#######################
tumor5=SetIdent(tumor5,value = "cl2")
library(COSG)
marker_anno=cosg(tumor5,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 10,expressed_pct = 0.2)
gene_anno=marker_anno$name
#write.csv(gene_anno,"topgene.csv")
#gene_anno=read.csv("topgene.csv",header = T,row.names = 1)
library(tidyr)
long_gene_anno<-gather(gene_anno, cell, gene, 1:6)
markers= markers=long_gene_anno$gene
markers=markers[!duplicated(markers)]

library(ggplot2)
DotPlot(tumor5, features = markers)+
  coord_flip()+
  theme_bw()+ # Use a clean theme and rotate the plot.
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+ # Rotate axis labels by 90 degrees.
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#f6fbef","#b2e1ba","#49add0","#094687"))+ # Set the color gradient.
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) +scale_size(range = c(0,3))

library(Nebulosa)
plot_density(tumor5,features=c("BICC1","CFHR1","EGR1","SEPP1","TOP2A","KRT19"),raster = T)

######cytotrace############
rm(list=ls())
# Install CytoTRACE2 during environment setup; see the repository README.
library(CytoTRACE2)
library(tidyverse)
library(Seurat)
library(qs)
# Use HCC_TME_DATA_DIR for external inputs; no local working directory is embedded.

tumor_hcc_noepi_downsample=subset(tumor_hcc_noepi,downsample=20000)

# Provide the input Seurat object.
cytotrace2_result_sce <- cytotrace2(tumor_hcc_noepi_downsample,
                                    is_seurat = TRUE,
                                    slot_type = "counts",
                                    species = 'human',
                                    seed = 1234)
qsave(cytotrace2_result_sce,"cytotrace2_result_sce.qs")

# making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame(phenotype = tumor_hcc_noepi_downsample@meta.data$cl2) %>%
  set_rownames(., colnames(tumor_hcc_noepi_downsample))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result_sce,
                  annotation = annotation,
                  is_seurat = TRUE)
# Analysis note.
p1 <- plots$CytoTRACE2_UMAP
# Analysis note.
p2 <- plots$CytoTRACE2_Potency_UMAP
# Analysis note.
p3 <- plots$CytoTRACE2_Relative_UMAP
# Use the cell-type annotation.
p4 <- plots$CytoTRACE2_Boxplot_byPheno

library(patchwork)
(p1+p2+p3+p4) + plot_layout(ncol = 2)

qsave(cytotrace2_result_sce,"cytotrace2_result_sce.qs")


plot_density(cytotrace2_result_sce,features="CytoTRACE2_Score",raster = T)

#######monocle###############
##########monocle2#######################

tumor_hcc_noepi_downsample=tumor_hcc_noepi[,sample(colnames(tumor_hcc_noepi), 20000)]

library(monocle)
sub_monocle=tumor_hcc_noepi_downsample
data <- as(as.matrix(sub_monocle@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sub_monocle@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,#as(as.matrix(data),"sparseMatrix")
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
# Choose the expression family that matches the matrix type.
# Use negbinomial.size() for sparse counts, tobit() for FPKM, and gaussianff() for log-FPKM.

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)


monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))
k = fData(monocle_cds)$num_cells_expressed>=5;table(k) # Filter genes expressed in fewer than five cells.
monocle_cds=monocle_cds[k,]
monocle_cds

# Select representative genes.
# Select highly dispersed genes.
# Use the Seurat method.
#FindVariableFeatures(sub_monocle, selection.method = "vst", nfeatures = 2000)
var.genes <- VariableFeatures(sub_monocle)
monocle_cds <- setOrderingFilter(monocle_cds, var.genes)
p2 <- plot_ordering_genes(monocle_cds)

# Use the Monocle method.
#disp_table <- dispersionTable(monocle_cds)
#disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
#monocle_cds <- setOrderingFilter(monocle_cds, disp.genes)
#p3 <- plot_ordering_genes(monocle_cds)

# Identify differentially expressed genes.
#diff.wilcox = FindAllMarkers(sub_monocle)
#all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
#diff.genes <- subset(all.markers,p_val_adj<0.01)$gene
#monocle_cds <- setOrderingFilter(monocle_cds, diff.genes)
#p1 <- plot_ordering_genes(monocle_cds)

#_____________________________________________________________

# Perform dimensionality reduction.
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, method = 'DDRTree')
# Order cells along the trajectory.
monocle_cds <- orderCells(monocle_cds,reverse = T)
qsave(monocle_cds,"monocle2_result.qs")

# Visualize the trajectory.
library(ggsci)
p1 = plot_cell_trajectory(monocle_cds, color_by = "cl2")+scale_color_npg()

plot_cell_trajectory(monocle_cds, color_by = "AUC")+scale_color_steps(low = "blue",high = "red", breaks = c(0,0.01,0.02,0.03,0.04, 0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14, 0.15))

# Display the plot using facets.
p2 = plot_cell_trajectory(monocle_cds, color_by = "cl5") + facet_wrap(~cl5, nrow = 2)+scale_color_npg()
p1|p2

p3=plot_complex_cell_trajectory(monocle_cds,x=1,y=2,color_by = "RNA_snn_res.0.4")+
  scale_color_manual(values = pal)+theme(legend.title = element_blank())
library(ggpubr)
df=pData(monocle_cds)
p4=ggplot(df,aes(Pseudotime,colour=cl2,fill=cl2))+
  geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic2()+scale_color_npg()

# Visualize the trajectory.
p5 <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")
(p1|p5)/p4

# Visualize the results.
s.genes <- c("VEGFA","CXCR4")
# Analysis note.
plot_genes_jitter(monocle_cds[s.genes,], grouping = "Pseudotime", color_by = "Pseudotime")+scale_color_gradient(low = "white",high = "#D7232A")
# Analysis note.
plot_genes_violin(monocle_cds[s.genes,], grouping = "RNA_snn_res.0.4", color_by = "RNA_snn_res.0.4")+scale_color_manual(values = pal)
# Visualize pseudotime.
plot_genes_in_pseudotime(monocle_cds[s.genes,], color_by = "RNA_snn_res.0.4")+scale_color_manual(values = pal)


# Analysis note.
disp_table <- dispersionTable(monocle_cds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
mycds_sub <- monocle_cds[disp.genes,]
plot_cell_trajectory(mycds_sub, color_by = "State")
beam_res <- BEAM(mycds_sub, branch_point = 3, cores = 2)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-4)),]
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 3, num_clusters = 3, show_rownames = T)

# Select the top 100 results by q value.
BEAM_genes=top_n(beam_res,n=100,desc(qval)) %>% pull(gene_short_name) %>% as.character()
p=plot_genes_branched_heatmap(mycds_sub[BEAM_genes,],branch_point = 3, num_clusters = 3, show_rownames = T)

# Plot an individual gene.
genes=row.names(subset(fData(mycds_sub),gene_short_name %in% c("HSP90AA1","DNAJB1")))
plot_genes_branched_pseudotime(mycds_sub[genes,],branch_point = 3,color_by = "State",ncol = 1)

# Analyze pathways.
neu=qread("../neu_hcc_pb_HBV.qs")
monocle_cds=qread(file.path(data_dir, "monocle2_result.qs"))

library(msigdbr)
A=msigdbr_collections()
msigdbr_show_species() # Analysis note.
h.human <- msigdbr(species="Homo sapiens",category="H")

h.names <- unique(h.human$gs_name)

h.sets <- vector("list",length=length(h.names))
names(h.sets) <- h.names

for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
}

sce = AddModuleScore(object = neu, features =
                       list(test = h.sets$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION),
                     name = "EMT")
input.data = data.frame(celltype = monocle_cds$cl5,
                        Pseudotime = monocle_cds$Pseudotime,
                        ID=colnames(monocle_cds))
A=sce$EMT1
A=as.data.frame(A)
A$ID=1:nrow(A)
A=A[match(input.data$ID,rownames(A)),]
input.data$test_pathway = A$A

# Define the plot theme.
if(T){
  text.size = 12
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"),
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   panel.grid=element_blank(), # Remove grid lines.
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
  )
}
ggplot(input.data, aes(x = Pseudotime, y = test_pathway)) +
  labs(x="Pseudotime",y = "EMT")+
  ggpubr::stat_cor(label.sep = "\n",
                   label.y = 0.09,
                   label.y.npc = "top",
                   size = 4,
                   method = "sp")  +
  geom_smooth(method='loess',size=0.8, color = "black") + theme_bw() +
  ggfun::facet_set(label = "EMT")+
  mytheme + theme(legend.position = "none")



sce = AddModuleScore(object = neu, features =
                       list(test = h.sets$HALLMARK_HEDGEHOG_SIGNALING),
                     name = "Hedgehog_signaling")
input.data = data.frame(celltype = monocle_cds$cl5,
                        Pseudotime = monocle_cds$Pseudotime,
                        ID=colnames(monocle_cds))
A=sce$Hedgehog_signaling1
A=as.data.frame(A)
A$ID=1:nrow(A)
A=A[match(input.data$ID,rownames(A)),]
input.data$test_pathway = A$A
ggplot(input.data, aes(x = Pseudotime, y = test_pathway)) +
  labs(x="Pseudotime",y = "Hedgehog signaling")+
  ggpubr::stat_cor(label.sep = "\n",
                   label.y = 0.04,
                   label.y.npc = "top",
                   size = 4,
                   method = "sp")  +
  geom_smooth(method='loess',size=0.8, color = "black") + theme_bw() +
  ggfun::facet_set(label = "EMT")+
  mytheme + theme(legend.position = "none")



sce = AddModuleScore(object = neu, features =
                       list(test = h.sets$HALLMARK_APICAL_SURFACE),
                     name = "Apical_surface")
input.data = data.frame(celltype = monocle_cds$cl5,
                        Pseudotime = monocle_cds$Pseudotime,
                        ID=colnames(monocle_cds))
A=sce$Apical_surface1
A=as.data.frame(A)
A$ID=1:nrow(A)
A=A[match(input.data$ID,rownames(A)),]
input.data$test_pathway = A$A
ggplot(input.data, aes(x = Pseudotime, y = test_pathway)) +
  labs(x="Pseudotime",y = "Apical surface")+
  ggpubr::stat_cor(label.sep = "\n",
                   label.y = 0.1,
                   label.y.npc = "top",
                   size = 4,
                   method = "sp")  +
  geom_smooth(method='loess',size=0.8, color = "black") + theme_bw() +
  ggfun::facet_set(label = "EMT")+
  mytheme + theme(legend.position = "none")



sce = AddModuleScore(object = neu, features =
                       list(test = h.sets$HALLMARK_ANGIOGENESIS),
                     name = "Angiogenesis")
input.data = data.frame(celltype = monocle_cds$cl5,
                        Pseudotime = monocle_cds$Pseudotime,
                        ID=colnames(monocle_cds))
A=sce$Angiogenesis1
A=as.data.frame(A)
A$ID=1:nrow(A)
A=A[match(input.data$ID,rownames(A)),]
input.data$test_pathway = A$A
ggplot(input.data, aes(x = Pseudotime, y = test_pathway)) +
  labs(x="Pseudotime",y = "Agiogenesis")+
  ggpubr::stat_cor(label.sep = "\n",
                   label.y = 0.05,
                   label.y.npc = "top",
                   size = 4,
                   method = "sp")  +
  geom_smooth(method='loess',size=0.8, color = "black") + theme_bw() +
  ggfun::facet_set(label = "EMT")+
  mytheme + theme(legend.position = "none")



sce = AddModuleScore(object = neu, features =
                       list(test = h.sets$HALLMARK_NOTCH_SIGNALING),
                     name = "Notch_signaling")
input.data = data.frame(celltype = monocle_cds$cl5,
                        Pseudotime = monocle_cds$Pseudotime,
                        ID=colnames(monocle_cds))
A=sce$Notch_signaling1
A=as.data.frame(A)
A$ID=1:nrow(A)
A=A[match(input.data$ID,rownames(A)),]
input.data$test_pathway = A$A
ggplot(input.data, aes(x = Pseudotime, y = test_pathway)) +
  labs(x="Pseudotime",y = "Notch signaling")+
  ggpubr::stat_cor(label.sep = "\n",
                   label.y = -0.05,
                   label.x = 5,
                   label.y.npc = "top",
                   size = 4,
                   method = "sp")  +
  geom_smooth(method='loess',size=0.8, color = "black") + theme_bw() +
  ggfun::facet_set(label = "EMT")+
  mytheme + theme(legend.position = "none")

###############survive_analysis######################
# Use HCC_TME_DATA_DIR for external inputs; no local working directory is embedded.
library(qs)
library(Seurat)
tumor=qread("tumor5_anno.qs")

DimPlot(tumor,group.by = "Cancer_type",raster = T)

tumor=SetIdent(tumor,value = "cl2")

library(COSG)
marker_anno=cosg(tumor,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 15,expressed_pct = 0.2)
gene_anno=marker_anno$name


LIHC=read.table(file.path(data_dir, "TCGA_LIHC_matrix.txt"),header = T,row.names = 1,sep = "\t")

geneset=as.list(gene_anno)
exp=as.matrix(LIHC)
# Use log2-transformed microarray data; for RNA-seq counts use kcdf = 'Poisson'.
library(GSVA)
re <- ssgseaParam(exp,geneSets = geneset)
re=gsva(re)
ssgsea_neu=as.data.frame(t(re))
cli=read.csv(file.path(data_dir, "LIHC_cli.csv"),header = T,row.names = 1)

ssgsea_neu=ssgsea_neu[substr(rownames(ssgsea_neu),16,16)=="A",]
rownames(ssgsea_neu)=substr(rownames(ssgsea_neu),1,15)
rownames(cli)=chartr("-",".",rownames(cli))

cli=cli[rownames(cli) %in% rownames(ssgsea_neu),]
ssgsea_neu=ssgsea_neu[match(rownames(cli),rownames(ssgsea_neu)),]

df=cbind(ssgsea_neu,cli)
#write.csv(df,"ssgsea_cli_tumor.csv")


library(survival)
library(survminer)
df$exp <- ''
res.cut=surv_cutpoint(
  df,
  time = "OS.time",
  event = "OS",
  variables="Ca_SEPP1",
  minprop = 0.1,
  progressbar = TRUE
)

df[df$Ca_SEPP1 >= res.cut$cutpoint$cutpoint,]$exp <- "H"
df[df$Ca_SEPP1 <  res.cut$cutpoint$cutpoint,]$exp <- "L"
# Fit the survival model.
fit <- survfit(Surv(OS.time, OS)~exp, data=df) # Fit the survival model.
# Display the P value.
surv_pvalue(fit)$pval.txt
# Plot the results.
#ggsurvplot(fit,pval=TRUE)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#e2145b", "#0a6b9d")# palette = c("#FF4500", "#4682B4")
)


library(survival)
library(plyr)
aa=df
a1=aa[1:6]
a1=scale(a1)
aa=aa[-(1:6)]
aa=cbind(a1,aa)
covariates <- colnames(aa)[c(1:6)]
# Overall survival analysis.
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = aa)})

# Extract the hazard ratio, 95% confidence interval, and P value.
univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
# Extract the P value.
                         p.value<-signif(x$wald["pvalue"], digits=2)
# Extract the hazard ratio.
                         HR <-signif(x$coef[2], digits=2);
# Extract the 95% confidence interval.
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })
# Convert to a data frame and transpose it.
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
write.csv(res, "cox_tumor.csv")
# Disease-free interval analysis.
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(DFI.time, DFI)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = aa)})

# Extract the hazard ratio, 95% confidence interval, and P value.
univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
# Extract the P value.
                         p.value<-signif(x$wald["pvalue"], digits=2)
# Extract the hazard ratio.
                         HR <-signif(x$coef[2], digits=2);
# Extract the 95% confidence interval.
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })
# Convert to a data frame and transpose it.
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

############GSE14520##############
load(file.path(data_dir, "GSE14520_data_surv.Rdata"))
exp=as.matrix(GSE14520_data)
# Use log2-transformed microarray data; for RNA-seq counts use kcdf = 'Poisson'.
library(GSVA)
re <- ssgseaParam(exp,geneSets = geneset)
re=gsva(re)
ssgsea_neu=as.data.frame(t(re))
cli=GSE14520_surv

df=cbind(ssgsea_neu,cli)

library(survival)
library(survminer)
df$exp <- ''
res.cut=surv_cutpoint(
  df,
  time = "OS.time",
  event = "OS",
  variables="Ca_TOP2A",
  minprop = 0.1,
  progressbar = TRUE
)

df[df$Ca_TOP2A >= res.cut$cutpoint$cutpoint,]$exp <- "H"
df[df$Ca_TOP2A <  res.cut$cutpoint$cutpoint,]$exp <- "L"
# Fit the survival model.
fit <- survfit(Surv(OS.time, OS)~exp, data=df) # Fit the survival model.
# Display the P value.
surv_pvalue(fit)$pval.txt
# Plot the results.
#ggsurvplot(fit,pval=TRUE)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#e2145b", "#0a6b9d")# palette = c("#FF4500", "#4682B4")
)



############GSE116174##############
load(file.path(data_dir, "GSE116174_data_surv.Rdata"))
exp=as.matrix(GSE116174_data)
# Use log2-transformed microarray data; for RNA-seq counts use kcdf = 'Poisson'.
library(GSVA)
re <- ssgseaParam(exp,geneSets = geneset)
re=gsva(re)
ssgsea_neu=as.data.frame(t(re))
cli=GSE116174_surv

df=cbind(ssgsea_neu,cli)
df$OS.time=as.numeric(df$OS.time)
df$OS=as.numeric(df$OS)
library(survival)
library(survminer)
df$exp <- ''
res.cut=surv_cutpoint(
  df,
  time = "OS.time",
  event = "OS",
  variables="Ca_TOP2A",
  minprop = 0.1,
  progressbar = TRUE
)

df[df$Ca_TOP2A >= res.cut$cutpoint$cutpoint,]$exp <- "H"
df[df$Ca_TOP2A <  res.cut$cutpoint$cutpoint,]$exp <- "L"
# Fit the survival model.
fit <- survfit(Surv(OS.time, OS)~exp, data=df) # Fit the survival model.
# Display the P value.
surv_pvalue(fit)$pval.txt
# Plot the results.
#ggsurvplot(fit,pval=TRUE)

ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#e2145b", "#0a6b9d")# palette = c("#FF4500", "#4682B4")
)




############GSE40873##############
load(file.path(data_dir, "GSE40873_data_surv.Rdata"))
exp=as.matrix(GSE40873_data)
# Use log2-transformed microarray data; for RNA-seq counts use kcdf = 'Poisson'.
library(GSVA)
re <- ssgseaParam(exp,geneSets = geneset)
re=gsva(re)
ssgsea_neu=as.data.frame(t(re))
cli=GSE40873_surv

df=cbind(ssgsea_neu,cli)

library(survival)
library(survminer)
df$exp <- ''
res.cut=surv_cutpoint(
  df,
  time = "OS.time",
  event = "OS",
  variables="Ca_TOP2A",
  minprop = 0.1,
  progressbar = TRUE
)

df[df$Ca_TOP2A >= res.cut$cutpoint$cutpoint,]$exp <- "H"
df[df$Ca_TOP2A <  res.cut$cutpoint$cutpoint,]$exp <- "L"
# Fit the survival model.
fit <- survfit(Surv(OS.time, OS)~exp, data=df) # Fit the survival model.
# Display the P value.
surv_pvalue(fit)$pval.txt
# Plot the results.
#ggsurvplot(fit,pval=TRUE)

ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#e2145b", "#0a6b9d")# palette = c("#FF4500", "#4682B4")
)
