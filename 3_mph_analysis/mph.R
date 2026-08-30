# Use HCC_TME_DATA_DIR for external inputs; no local working directory is embedded.
data_dir <- Sys.getenv("HCC_TME_DATA_DIR", unset = "data")

myeloid=qread(file.path(data_dir, "myeloid2_anno.qs"))
DimPlot(myeloid,raster = T,group.by = "cl2",cols = SC_color)
mph=subset(myeloid,cl3=="Mph")

mph <- NormalizeData(mph, normalization.method = "LogNormalize", scale.factor = 10000)
mph <- FindVariableFeatures(mph, selection.method = "vst", nfeatures = 1500)
mph <- ScaleData(mph)

mph <- RunPCA (mph, features = VariableFeatures(object = mph), ndims.print = 1:2)
ElbowPlot(object = mph, ndims = 100)
library(harmony)

mph2 = mph %>% RunHarmony("Sample", plot_convergence = TRUE)
set.resolutions <- seq(0.5, 0.7, by = 0.1)
pdf(file = "PCA-mph_harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

mph2 <- FindNeighbors(mph2, dims = 1:10,reduction = "harmony")
mph2  <- FindClusters(object = mph2 , resolution = set.resolutions, verbose = FALSE)
mph2  <- RunUMAP(mph2 , dims = 1:10,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = mph2, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()

DimPlot(mph2,group.by = "cl2",raster = T,cols = SC_color)

##DoubletFinder############
library(DoubletFinder)
library(BiocParallel)
library(qs)
library(Seurat)

A=as.data.frame(table(mph2$Sample))
A=A[A$Freq<100,]
Idents(mph2)=mph2$Sample
mph2_filt=subset(mph2,idents=c(A$Var1),invert=T)
A=as.data.frame(table(mph2_filt$Sample))

# Process samples separately for DoubletFinder.
sce_list <- SplitObject(mph2_filt, split.by = "Sample")

pc.num <- 1:20
for (z in 1:nrow(A)) {
  DoubletRate = ifelse(ncol(sce_list[[as.character(A$Var1[z])]])<800,0.004,
                       ifelse(ncol(sce_list[[as.character(A$Var1[z])]])<1600,0.008,
                              ifelse(ncol(sce_list[[as.character(A$Var1[z])]])<3200,0.016,
                                     ifelse(ncol(sce_list[[as.character(A$Var1[z])]])<4800,0.023,
                                            ifelse(ncol(sce_list[[as.character(A$Var1[z])]])<6400,0.031,
                                                   ifelse(ncol(sce_list[[as.character(A$Var1[z])]])<8000,0.039,
                                                          ifelse(ncol(sce_list[[as.character(A$Var1[z])]])<9600,0.046,
                                                                 ifelse(ncol(sce_list[[as.character(A$Var1[z])]])<11200,0.054,
                                                                        ifelse(ncol(sce_list[[as.character(A$Var1[z])]])<12800,0.061,
                                                                               ifelse(ncol(sce_list[[as.character(A$Var1[z])]])<14400,0.069,
                                                                                      ifelse(ncol(sce_list[[as.character(A$Var1[z])]])<16000,0.076,0.076)))))))))))
# Identify the optimal pK value.
  sweep.res <- paramSweep_v3(sce_list[[as.character(A$Var1[z])]], PCs = pc.num,sct = F) # Analysis note.
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

# Estimate the homotypic doublet proportion and the expected number of doublets.
  homotypic.prop <- modelHomotypic(sce_list[[as.character(A$Var1[z])]]$seurat_clusters) # Analysis note.
  nExp_poi <- round(DoubletRate * ncol(sce_list[[as.character(A$Var1[z])]]))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Detect doublets using the selected parameters.
  sce_list[[as.character(A$Var1[z])]] <- doubletFinder_v3(sce_list[[as.character(A$Var1[z])]],
                                                          PCs = pc.num,
                                                          pN = 0.25,
                                                          pK = pK_bcmvn,
                                                          nExp = nExp_poi.adj,
                                                          reuse.pANN = F,
                                                          sct = F) # Analysis note.

# Visualize the results.
  #DimPlot(sce_list[[as.character(A$Var1[z])]],
  #        reduction = "umap",
  #        group.by = "DF.classifications_0.25_0.005_6")
  sce_list[[as.character(A$Var1[z])]]$doublet=sce_list[[as.character(A$Var1[z])]]@meta.data[ncol(mph@meta.data)+2]
# Assume the Seurat object list was created with SplitObject.

}
seurat_combined <- Reduce(function(x, y) merge(x, y, add.cell.ids = c("x", "y"), project = "combined_project"), sce_list)



mph=seurat_combined
#############cellratio_doublet#########
library(ggplot2)
library(dplyr)
library(ggalluvial)

Ratio <- mph2@meta.data %>%
  group_by(cl2, doublet) %>% # Group the observations.
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))

ggplot(Ratio, aes(x =cl2, y= relative_freq, fill = doublet,
                  stratum=doublet, alluvium=doublet)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.3, knot.pos=0.2)+ # Set knot.pos to control the curvature of the alluvial connections.
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  scale_fill_manual(values = SC_color[5:6])+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) # Rotate axis labels by 90 degrees.

qsave(mph2,"mph_doublet.qs")


mph2=SetIdent(mph2,value = "cl2")
library(COSG)
marker_anno=cosg(mph2,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 30)
gene_anno=marker_anno$name

# Correct batch effects by cohort.
mph2 = mph2 %>% RunHarmony("cohort", plot_convergence = TRUE)
set.resolutions <- seq(0.5, 0.7, by = 0.1)
pdf(file = "PCA-mph_harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

mph2 <- FindNeighbors(mph2, dims = 1:10,reduction = "harmony")
mph2  <- FindClusters(object = mph2 , resolution = set.resolutions, verbose = FALSE)
mph2  <- RunUMAP(mph2 , dims = 1:10,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = mph2, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
qsave(mph2,"mph2.qs")

#############marker_exp####################
markers=c("FCN1","CD247","IGF1","MARCO","SPP1","CXCL9")

library(ggplot2)
DotPlot(mph2, features = markers,group.by = "cl2")+coord_flip()+theme_bw()+ # Use a clean theme and rotate the plot.
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+ # Rotate axis labels by 90 degrees.
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#f6fbef","#b2e1ba","#49add0","#094687"))+ # Set the color gradient.
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
meta.tb_hbv$patient=meta.tb_hbv$Sample

########mph2_filt############
mph2_filt=subset(mph2,cl2=="Mph_CD68_CD247",invert=T)

mph2_filt <- NormalizeData(mph2_filt, normalization.method = "LogNormalize", scale.factor = 10000)
mph2_filt <- FindVariableFeatures(mph2_filt, selection.method = "vst", nfeatures = 2000)
mph2_filt <- ScaleData(mph2_filt)

mph2_filt <- RunPCA (mph2_filt, features = VariableFeatures(object = mph2_filt), ndims.print = 1:2)
ElbowPlot(object = mph2_filt, ndims = 100)
library(harmony)

mph2_filt2 = mph2_filt %>% RunHarmony("cohort", plot_convergence = TRUE)
set.resolutions <- seq(0.2, 0.7, by = 0.1)
pdf(file = "PCA-mph2_filt_harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

mph2_filt2 <- FindNeighbors(mph2_filt2, dims = 1:20,reduction = "harmony")
mph2_filt2  <- FindClusters(object = mph2_filt2 , resolution = set.resolutions, verbose = FALSE)
mph2_filt2  <- RunUMAP(mph2_filt2 , dims = 1:20,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = mph2_filt2, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()

Idents(mph2_filt2)=mph2_filt2$RNA_snn_res.0.3
library(COSG)
marker_anno=cosg(mph2_filt2,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 30)
gene_anno=marker_anno$name

Idents(mph2_filt2)=mph2_filt2$RNA_snn_res.0.3
marker=data.frame(cluster=0:19,cell=0:19)
marker[marker$cluster %in% c(16,4,14),2]<-"Mph_MARCO"
marker[marker$cluster %in% c(3,11,0,8,11,13,15,19),2]<-"Mph_IGF1"
marker[marker$cluster %in% c(5,12,18,2),2]<-"Mph_SPP1"
marker[marker$cluster %in% c(7),2]<-"Mph_CXCL9"
marker[marker$cluster %in% c(6,9,1,10,17),2]<-"Mono_like_Mph"


mph2_filt2@active.ident=factor(mph2_filt2@active.ident,levels = 0:19)
mph2_filt2@meta.data$cl3=sapply(mph2_filt2@active.ident,function(x){marker[x,2]})
DimPlot(mph2_filt2,reduction = "umap",group.by = "cl3",label = F,cols = SC_color[10:40],raster = T)

# Calculate mean gene-set expression.
geneset=read.csv(file.path(data_dir, "Myeloid.csv"),header = F)
geneset=geneset[-2]
geneset=as.data.frame(t(geneset))
colnames(geneset)=geneset[1,]
geneset=geneset[-1,]


aver_dt<- AverageExpression(mph2_filt2,
                            #features = gene,
                            group.by = 'cl3',
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


aver_dt<- AverageExpression(mph2_filt2,
                            #features = gene,
                            group.by = 'Cancer_type',
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

##################UCELL_analysis######################################
geneset=as.list(geneset)

# Check packages required for UCell analysis.
cran.packages <- c("msigdbr", "dplyr", "purrr", "stringr","magrittr",
                   "RobustRankAggreg", "tibble", "reshape2", "ggsci",
                   "tidyr", "aplot", "ggfun", "ggplotify", "ggridges",
                   "gghalves", "Seurat", "SeuratObject", "methods",
                   "devtools", "BiocManager","data.table","doParallel",
                   "doRNG")

# install packages from Bioconductor
bioconductor.packages <- c("GSEABase", "AUCell", "SummarizedExperiment",
                           "singscore", "GSVA", "ComplexHeatmap", "ggtree",
                           "Nebulosa")
ucell.packages <- c(cran.packages, bioconductor.packages, "UCell", "irGSEA")
missing.packages <- ucell.packages[
  !vapply(ucell.packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing.packages)) {
  stop("Install the missing packages listed in README.md: ",
       paste(missing.packages, collapse = ", "))
}


library(UCell)
library(irGSEA)
geneset_select=geneset
mph2_filt2$Cancer_type=factor(mph2_filt2$Cancer_type,levels = c("Healthy","HBV","AL","HCC"))

mph2_filt.final <- irGSEA.score(object = mph2_filt2, assay = "RNA",
                                slot = "data", seeds = 123, ncores = 10,
                                min.cells = 3, min.feature = 0,
                                custom = T, geneset = geneset_select, msigdb = F,
                                species = "Homo sapiens", category = "H",
                                subcategory = NULL, geneid = "symbol",
                                method = c("UCell"),
                                #method = c("AUCell", "UCell", "singscore","ssgsea"),
                                aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                                kcdf = 'Gaussian')
A=mph2_filt.final@assays$UCell@scale.data
A[1:5,1:5]
A=as.data.frame(t(A))
mph2_filt.final$M2=A$`M2-like`
mph2_filt.final$M1=A$`M1-like`

Idents(mph2_filt.final)=mph2_filt.final$cl3
Idents(mph2_filt.final)=mph2_filt.final$Cancer_type

DotPlot(mph2_filt.final,features = "M1",group.by = "cl2")
VlnPlot(mph2_filt.final,features = "M1",group.by = "cl2",pt.size = 0)
irGSEA.densityheatmap(object = mph2_filt.final,
                      method = "UCell",
                      show.geneset = "M2-like",
                      heatmap_width = 6,
                      heatmap_height = 6)
irGSEA.density.scatterplot(object =mph2_filt.final,
                           method = "UCell",
                           show.geneset = "M2-like",
                           reduction = "umap")

mph


Idents(mph2_filt.final_filt)=mph2_filt.final_filt$module
irGSEA.densityheatmap(object = mph2_filt.final_filt,
                      method = "UCell",
                      show.geneset = "M2-like",
                      heatmap_width = 6,
                      heatmap_height = 6)
Idents(mph2_filt.final)=mph2_filt.final$Sample
mph2_filt.final_filt=subset(mph2_filt.final,Cancer_type=="HCC")
A=as.data.frame(table(mph2_filt.final_filt$Sample))
module=read.csv(file.path(data_dir, "cluster4.csv"),header = T)
module=module[match(A$Var1,module$X),]
mph2_filt.final_filt@active.ident=factor(mph2_filt.final_filt@active.ident,levels = c(module$X))
mph2_filt.final_filt@meta.data$module=sapply(mph2_filt.final_filt@active.ident,function(x){module[x,2]})


##############OR_AL_HBVa_HBVi###########
library("sscVis")
library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")
# Set the figure output directory.
out.prefix <- "./"
meta.tb=mph2_filt.final_filt
meta.tb$celltype=meta.tb$cl3
meta.tb$patient=meta.tb$Sample
do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          celltype = cellInfo.tb$celltype,
                          colname.patient = "patient",
                          module = cellInfo.tb$module,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data
  library(data.table)
  dir.create(dirname(out.prefix),F,T)

  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$celltype = as.character(celltype)

  if(is.factor(module)){
    cellInfo.tb$module = module
  }else{cellInfo.tb$module = as.factor(module)}

  module.avai.vec <- levels(cellInfo.tb[["module"]])
  count.dist <- unclass(cellInfo.tb[,table(celltype,module)])[,module.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)

  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }

  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}

OR.M.list <- do.tissueDist(cellInfo.tb=meta.tb@meta.data,
                           out.prefix=sprintf("%s.STARTRAC.dist.M.baseline",out.prefix),
                           pdf.width=4,pdf.height=6,verbose=1)
OR=OR.M.list$OR.dist.mtx # OR > 1.5 indicates enrichment and OR < 0.5 indicates depletion; P values are BH-adjusted.
Pvalue=OR.M.list$p.dist.tb
library(pheatmap)
n=OR
n[n>2]=2
color=colorRampPalette(c("white",'#FEE7CE','#FEC18C','#E75519'))(100)
#color=colorRampPalette(c('#729bb8','white','#bb605f'))(100)
pheatmap(n,color = color,cluster_cols = F,cluster_rows = F,display_numbers = T)

par(las = 2) # Rotate plot labels.
par(mar = c(5, 8, 4, 2))
num_hbvstate=table(meta.tb$module)
num_cl3=table(meta.tb$cl3)
barplot(num_hbvstate,col = SC_color[c(9,8,7)])
barplot(num_cl3,col = SC_color[5:30],horiz = T)

library(ggplot2)
library(dplyr)
library(ggalluvial)

Ratio <- meta.tb@meta.data %>%
  dplyr::group_by(module, cl3) %>% # Group the observations.
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(relative_freq = n/sum(n))

ggplot(Ratio, aes(x =module, y= relative_freq, fill = cl3,
                  stratum=cl3, alluvium=cl3)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.3, knot.pos=0.2)+ # Set knot.pos to control the curvature of the alluvial connections.
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  scale_fill_manual(values = SC_color[10:30])
qsave(mph2_filt.final_filt,"mph_subtype.qs")

#######IL1B_analysis#########

mph2_filt.final <- NormalizeData(mph2_filt.final, normalization.method = "LogNormalize", scale.factor = 10000)
mph2_filt.final <- FindVariableFeatures(mph2_filt.final, selection.method = "vst", nfeatures = 1500)
mph2_filt.final <- ScaleData(mph2_filt.final)

mph2_filt.final <- RunPCA (mph2_filt.final, features = VariableFeatures(object = mph2_filt.final), ndims.print = 1:2)

ElbowPlot(object = mph2_filt.final, ndims = 100)
library(harmony)
mph2_filt.final2 = mph2_filt.final %>% RunHarmony("Sample", plot_convergence = TRUE)

mph2_filt.final2 <- FindNeighbors(mph2_filt.final2, dims = 1:20,reduction = "harmony")
mph2_filt.final2  <- FindClusters(object = mph2_filt.final2 , resolution = 0.5, verbose = FALSE)
mph2_filt.final2  <- RunUMAP(mph2_filt.final2 , dims = 1:20,reduction = "harmony")

DimPlot(mph2_filt.final2,group.by = "cl3",raster = T,cols = SC_color)
DimPlot(mph2_filt.final2,group.by = "cohort",raster = T,cols = SC_color)

library(COSG)
Idents(mph2_filt.final2)=mph2_filt.final2$RNA_snn_res.0.5
marker_anno=cosg(mph2_filt.final2,groups="all",assay="RNA",slot="data",mu=1,n_genes_user = 30)
gene_anno=marker_anno$name

#########IL1B
IL1B=mph2_filt.final2@assays$RNA$data["IL1B",]
IL1B=as.data.frame(IL1B)
meta=mph2_filt.final2$orig.ident
meta=as.data.frame(meta)
meta$IL1B=IL1B$IL1B
meta$IL1B_type=""
meta[meta$IL1B > median(meta$IL1B),]$IL1B_type="IL1B_high"
meta[meta$IL1B <= median(meta$IL1B),]$IL1B_type="IL1B_low"
mph2_filt.final2$IL1B_type=meta$IL1B_type
DimPlot(mph2_filt.final2,group.by = "IL1B_type",cols = c("#0a6b9d","grey"),raster = T)
VlnPlot(mph2_filt.final2,features = "IL1B",group.by = "IL1B_type",cols = c("#0a6b9d","grey"),pt.size = 0)
qsave(mph2_filt.final2,"mph2_filt_IL1Btype.qs")

meta$cl3=mph2_filt.final2$cl3
meta=meta[order(meta$IL1B,decreasing = T),]
write.csv(meta[1:100,],"mph_top100_IL1B.csv")
table(meta[1:100,]$cl3)
#############cellratio_IL1B#########
library(ggplot2)
library(dplyr)
library(ggalluvial)

Ratio <- mph2_filt.final2@meta.data %>%
  group_by(IL1B_type, cl3) %>% # Group the observations.
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))

ggplot(Ratio, aes(x =IL1B_type, y= relative_freq, fill = cl3,
                  stratum=cl3, alluvium=cl3)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.3, knot.pos=0.2)+ # Set knot.pos to control the curvature of the alluvial connections.
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  scale_fill_manual(values = SC_color[5:20])+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) # Rotate axis labels by 90 degrees.

########ucell##########
Idents(mph2_filt.final2)=mph2_filt.final2$IL1B_type
irGSEA.densityheatmap(object = mph2_filt.final2,
                                        method = "UCell",
                                        show.geneset = "M2-like",heatmap_width = 4,heatmap_height = 7)
irGSEA.densityheatmap(object = mph2_filt.final2,
                      method = "UCell",
                      show.geneset = "M1-like",heatmap_width = 4,heatmap_height = 7)

irGSEA.halfvlnplot(object = mph2_filt.final2,
                              method = "UCell",
                              show.geneset = "Hypoxia",color.cluster = c("#0a6b9d","grey"))

########DEG_HCC_IL1B###############
Idents(mph2_filt.final2)=mph2_filt.final2$IL1B_type
deg_mph2_filt.final2=FindMarkers(mph2_filt.final2,ident.1 = "IL1B_high",ident.2 = "IL1B_low")
write.csv(deg_mph2_filt.final2,"deg_mph2_IL1B.csv")
library(EnhancedVolcano)
EnhancedVolcano(deg_mph2_filt.final2,lab = rownames(deg_mph2_filt.final2),x="avg_log2FC",y="p_val_adj",labSize = 4,max.overlaps = 400,FCcutoff = 1,
                  selectLab = c("IL1B","IL1A","CCL20","MMP10","TNF","NLRP3","PTGS2","CXCL1","CXCL2","CCL3","VEGFA","THBS1","PDGFB"),
                 pCutoff = 0.05)
library(ggplot2)
deg_mph2_filt.final2$def=deg_mph2_filt.final2$pct.1-deg_mph2_filt.final2$pct.2

p=ggplot(deg_mph2_filt.final2,aes(x=avg_log2FC,y=def))+
  geom_point(aes(size=abs(avg_log2FC), color= abs(avg_log2FC)))+
  scale_color_gradientn(values = seq(0,1,0.2),colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  theme_bw()

# Add labels.
library(ggrepel)
deg_mph2_filt.final2$label=ifelse(deg_mph2_filt.final2$p_val_adj<=0.05 & abs(deg_mph2_filt.final2$avg_log2FC)>=0.5 & abs(deg_mph2_filt.final2$def)>=0,as.character(rownames(deg_mph2_filt.final2)),'')
p+geom_text_repel(data=deg_mph2_filt.final2,aes(x=avg_log2FC, y=def,label=label),size=2.5,max.overlaps = 50)

#####gsea##########
library(clusterProfiler)
gene=deg_mph2_filt.final2
geneList=gene$avg_log2FC
names(geneList)=rownames(gene)
geneList = sort(geneList, decreasing = TRUE)
kegmt<-read.gmt("c5.go.v2024.1.Hs.symbols.gmt")
KEGG<-GSEA(geneList,TERM2GENE = kegmt,eps = 0)

library(enrichplot)
gseaplot2(KEGG,1,color="red",pvalue_table = T)
out=KEGG@result
write.csv(out,"GSEA_IL1B.csv")

######enrich#########
#############enrich_prevo_up50gene_gprofiler#################
library(ggplot2)
select_pathway=read.csv("./select_GO.csv")
mytheme<-theme(axis.title = element_text(size = 13),

               axis.text = element_text(size = 11),

               plot.title = element_text(size = 14, hjust= 0.5,face= "bold"),

               legend.title = element_text(size = 13),

               legend.text = element_text(size = 11))
mytheme2<-mytheme + theme(axis.text.y = element_blank()) # Hide y-axis text labels in the custom theme.

select_pathway$text_x <- rep(0.001,12) # Add a helper column to align text labels.
select_pathway$ID=factor(select_pathway$ID,levels = c(select_pathway$ID))
p2<- ggplot(data = select_pathway,aes(x = NES, y = ID)) +

  geom_bar(aes(fill = -log10(p.adjust)), stat = "identity", width = 0.8, alpha = 0.7) +

  scale_fill_distiller(palette = "Blues", direction = 1) +

  labs(x = "NES", y = "Terms", title = "GSEA") +

  geom_text(aes(x = text_x, # Use the helper column to control the text-label position.
                label= ID),hjust= 0)+ # Left-align text.
  theme_classic()+ mytheme2

p2




###############survive_analysis######################

LIHC=read.table(file.path(data_dir, "TCGA_LIHC_matrix.txt"),header = T,row.names = 1,sep = "\t")
deg=read.csv("deg_mph2_IL1B.csv",header = T)
geneset=list(up=deg$label[1:50],down=tail(deg$label))
exp=as.matrix(LIHC)
# Use log2-transformed microarray data; for RNA-seq counts use kcdf = 'Poisson'.
library(GSVA)
re=ssgseaParam(exp,geneset)
re2=gsva(re)

ssgsea_neu=as.data.frame(t(re2))
cli=read.csv(file.path(data_dir, "LIHC_cli.csv"),header = T,row.names = 1)

ssgsea_neu=ssgsea_neu[substr(rownames(ssgsea_neu),16,16)=="A",]
rownames(ssgsea_neu)=substr(rownames(ssgsea_neu),1,15)
rownames(cli)=chartr("-",".",rownames(cli))

cli=cli[rownames(cli) %in% rownames(ssgsea_neu),]
ssgsea_neu=ssgsea_neu[match(rownames(cli),rownames(ssgsea_neu)),]

df=cbind(ssgsea_neu,cli)
write.csv(df,"ssgsea_TCGALIHC_IL1B_top50gene.csv")


library(survival)
library(survminer)
df$exp <- ''
#res.cut=surv_cutpoint(
#  df,
#  time = "OS.time",
#  event = "OS",
#  variables="up",
#  minprop = 0.5,
#  progressbar = TRUE
#)

#df[df$up >= res.cut$cutpoint$cutpoint,]$exp <- "H"
#df[df$up <  res.cut$cutpoint$cutpoint,]$exp <- "L"

df[df$up >= median(df$up),]$exp <- "H"
df[df$up <  median(df$up),]$exp <- "L"

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
           palette = c("#0a6b9d", "gray")# palette = c("#FF4500", "#4682B4")
)

########DFI

# Fit the survival model.
fit <- survfit(Surv(DFI.time, DFI)~exp, data=df) # Fit the survival model.
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
           palette = c("#0a6b9d", "gray")# palette = c("#FF4500", "#4682B4")
)


########PFI
# Fit the survival model.
fit <- survfit(Surv(PFI.time, PFI)~exp, data=df) # Fit the survival model.
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
           palette = c("#0a6b9d", "gray")# palette = c("#FF4500", "#4682B4")
)

########DFI

# Fit the survival model.
fit <- survfit(Surv(DSS.time, DSS)~exp, data=df) # Fit the survival model.
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
           palette = c("#0a6b9d", "gray")# palette = c("#FF4500", "#4682B4")
)

############GSE14520##############
load(file.path(data_dir, "GSE14520_data_surv.Rdata"))
exp=as.matrix(GSE14520_data)
# Use log2-transformed microarray data; for RNA-seq counts use kcdf = 'Poisson'.
library(GSVA)
re=ssgseaParam(exp,geneset)
re2=gsva(re)
ssgsea_neu=as.data.frame(t(re2))
cli=GSE14520_surv

df=cbind(ssgsea_neu,cli)

library(survival)
library(survminer)
df$exp <- ''
res.cut=surv_cutpoint(
  df,
  time = "OS.time",
  event = "OS",
  variables="up",
  minprop = 0.1,
  progressbar = TRUE
)

df[df$up >= res.cut$cutpoint$cutpoint,]$exp <- "H"
df[df$up <  res.cut$cutpoint$cutpoint,]$exp <- "L"
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
           palette = c("#0a6b9d", "gray")# palette = c("#FF4500", "#4682B4")
)






############GSE40873##############
load(file.path(data_dir, "GSE40873_data_surv.Rdata"))
exp=as.matrix(GSE40873_data)
# Use log2-transformed microarray data; for RNA-seq counts use kcdf = 'Poisson'.
library(GSVA)
re=ssgseaParam(exp,geneset)
re2=gsva(re)
ssgsea_neu=as.data.frame(t(re2))
cli=GSE40873_surv

df=cbind(ssgsea_neu,cli)

library(survival)
library(survminer)
df$exp <- ''
res.cut=surv_cutpoint(
  df,
  time = "OS.time",
  event = "OS",
  variables="up",
  minprop = 0.1,
  progressbar = TRUE
)

df[df$up >= res.cut$cutpoint$cutpoint,]$exp <- "H"
df[df$up <  res.cut$cutpoint$cutpoint,]$exp <- "L"
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
           palette = c("#0a6b9d", "gray")# palette = c("#FF4500", "#4682B4")
)
############GSE116174##############
load(file.path(data_dir, "GSE116174_data_surv.Rdata"))
exp=as.matrix(GSE116174_data)
# Use log2-transformed microarray data; for RNA-seq counts use kcdf = 'Poisson'.
library(GSVA)
re=ssgseaParam(exp,geneset)
re2=gsva(re)
ssgsea_neu=as.data.frame(t(re2))
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
  variables="up",
  minprop = 0.1,
  progressbar = TRUE
)

df[df$up >= res.cut$cutpoint$cutpoint,]$exp <- "H"
df[df$up <  res.cut$cutpoint$cutpoint,]$exp <- "L"
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
           palette = c("#0a6b9d", "gray")# palette = c("#FF4500", "#4682B4")
)





####CIBERSORT########
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor","timeROC","pracma",
          "EPIC", "MCPcounter", "estimate", "IOBR")
missing.packages <- depens[
  !vapply(depens, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing.packages)) {
  stop("Install the missing packages listed in README.md: ",
       paste(missing.packages, collapse = ", "))
}

library(IOBR)
LIHC=read.table(file.path(data_dir, "TCGA_LIHC_matrix.txt"),header = T,row.names = 1,sep = "\t")
samplesTP<-LIHC[,substr(colnames(LIHC),14,15)<10]
samplesTP2=samplesTP[,substr(colnames(samplesTP),16,16)=="A"]
eset_lihc<-samplesTP2

# In this process, *biomaRt* R package is utilized to acquire the gene length of each Ensembl ID and calculate the TPM of each sample. If identical gene symbols exists, these genes would be ordered by the mean expression levels. The gene symbol with highest mean expression level is selected and remove others.
# NOTE: This process may take a few minutes which depends on the internet connection speed. Please wait for its completion.
eset_lihc<-count2tpm(countMat = eset_lihc, idType = "Symbol", org="hsa")
eset_lihc[1:5,1:5]
eset_lihc=log2(eset_lihc+1)

cibersort<-deconvo_tme(eset = eset_lihc, method = "cibersort", arrays = FALSE, perm = 1000 )



# Group the observations.
Ratio_sample <- mph2_filt.final2@meta.data %>%
  dplyr::group_by(IL1B_type, Sample) %>% # Group the observations.
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(relative_freq = n/sum(n))

IL1B_high=Ratio_sample[Ratio_sample$IL1B_type=="IL1B_high",]
c4=read.csv(file.path(data_dir, "cluster4.csv"),header = T)
c4=c4[c4$X %in% IL1B_high$Sample,]
IL1B_high=IL1B_high[match(c4$X,IL1B_high$Sample),]
IL1B_high=cbind(IL1B_high,c4)
write.csv(IL1B_high,"IL1B_high_rank_sample.csv")
