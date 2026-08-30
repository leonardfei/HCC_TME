# Use HCC_TME_DATA_DIR for external inputs; no local working directory is embedded.
data_dir <- Sys.getenv("HCC_TME_DATA_DIR", unset = "data")

library(qs)
library(Seurat)
TNK=qread("TNK2_anno.qs")
mye=qread("myeloid2_anno.qs")
stromal=qread("Stromal2_anno.qs")
BPC=qread("B_PC2_filt_anno.qs")
tumor=qread("tumor_epi.qs")

all=merge(TNK,y=c(mye,stromal,BPC,tumor))


all_subtype_anno <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all_subtype_anno <- FindVariableFeatures(all_subtype_anno, selection.method = "vst", nfeatures = 3000)
all_subtype_anno <- ScaleData(all_subtype_anno)

all_subtype_anno <- RunPCA (all_subtype_anno, features = VariableFeatures(object = all_subtype_anno), ndims.print = 1:2)

library(harmony)

all_subtype_anno2 = all_subtype_anno %>% RunHarmony("cohort", plot_convergence = TRUE)
set.resolutions <- seq(0.1, 0.3, by = 0.1)
pdf(file = "PCA_allsubtype-harmony.pdf")

#set.seed(101)
#library(future)
# Inspect the total number of available cores.
# Inspect the number of active workers.
#plan("multisession", workers = 2)
#options(future.globals.maxSize = 100 * 1024^3) # set 200G RAM

all_subtype_anno2 <- FindNeighbors(all_subtype_anno2, dims = 1:50,reduction = "harmony")
all_subtype_anno2  <- FindClusters(object = all_subtype_anno2 , resolution = set.resolutions, verbose = FALSE)
all_subtype_anno2  <- RunUMAP(all_subtype_anno2 , dims = 1:50,reduction = "harmony")

seurat_data.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = all_subtype_anno2, reduction = 'umap',raster = T,label = TRUE, group.by = paste0("RNA_snn_res.", x))
  print(p)
})
dev.off()
#plan("multisession", workers = 1)


A=all_subtype_anno2$cl3
A$cl2=all_subtype_anno2$cl2
A$cl1=all_subtype_anno2$cl1
A$A[A$cl1=="Malignant_Epi"]="Tumor_Epi"
all_subtype_anno2$cl3=A$A
A$cl2[A$cl1=="Malignant_Epi"]="Tumor_Epi"
all_subtype_anno2$cl2=A$cl2

DimPlot(all_subtype_anno2,group.by = "cl2",raster = T,cols = SC_color[3:65])
qsave(all_subtype_anno2,"all_subtype_anno2.qs")
library(paletteer)
color_name=palettes_d_names
pal=paletteer_d("khroma::smoothrainbow")

FeaturePlot(all_subtype_anno2,c("ACTB","PTPRC","ALB","AFP","EPCAM",
                        "APOA2","COL1A1","CDH5","CD3D","KLRF1",
                        "GNLY","MS4A1","MZB1","CD68","CD1C",
                        "FCGR3B","CLEC4C","TPSAB1"),cols = pal,raster = T,ncol = 6)

#######all_tumor_TME##############
all_tumor=subset(all_subtype_anno2,Cancer_type=="HCC")
all_tumor_TME=subset(all_tumor,cl3=="Tumor_Epi",invert=T)
library(ggplot2)
library(dplyr)
library(ggalluvial)

Ratio <- all_tumor_TME@meta.data %>%
  group_by(Sample, cl2) %>% # Group the observations.
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))

library(reshape2)
data_wide_d<-dcast(Ratio, Sample~Ratio$cl2,
                   value.var = 'relative_freq')
write.csv(data_wide_d,"cellratio_TME.csv")

library(pheatmap)
n=data_wide_d
rownames(n)=n$Sample
n=n[-1]
n=as.data.frame(t(n))

n=t(scale(t(n)))
n[n>1]=1
n[n<(-1)]=(-1)
color=colorRampPalette(c('#485494','white','#d6252a'))(100)
n[is.na(n)]=0
pheatmap(n,color = color,clustering_method = "ward.D")
write.csv(n,"heatmap_scale_allcelltype.csv")

annotation_row<-read.csv("heatmap_cell_anno.csv",header=T,row.names = 1)
# Import the annotated dataset.
annotation_row <- as.data.frame(annotation_row)
annotation_col<-read.csv("heatmap_sample_anno.csv",header = T,row.names = 1)
annotation_col <- as.data.frame(annotation_col)
ann_colors = list(
  CD45_type= c(CD45_neg="white", CD45_pos="tomato4"),
  cluster4= c(s1="#00AFBB",s2="#E7B800",s3="#d15c54",s4="#1e93c9"),
  cohort = c(c0=SC_color[21],c1=SC_color[22],c2=SC_color[23],c4=SC_color[24]),
  major_celltype = c(B=SC_color[5],CD4T=SC_color[6],CD8T=SC_color[7],
                     Cycling=SC_color[8],DC=SC_color[9],EC=SC_color[10],
                     Eos=SC_color[11],Fb=SC_color[12],gdT=SC_color[13],
                     ILC=SC_color[14],Mast=SC_color[15],Mph=SC_color[16],
                     Neu=SC_color[17],NK=SC_color[18],PC=SC_color[19])
)

pheatmap(n,color = color,clustering_method = "ward.D",annotation_row = annotation_row,annotation_col =annotation_col,annotation_colors = ann_colors )


# The data frame contains samples as columns and genes as rows.
A=pheatmap(n,color = color,clustering_method = "ward.D",annotation_row = annotation_row,annotation_col =annotation_col,annotation_colors = ann_colors )
df=n
df=as.data.frame(t(df))
# Perform hierarchical clustering.
hc <- A$tree_col

# Plot the clustering dendrogram.
dendrogram <- as.dendrogram(hc)
plot(dendrogram)

# Cut the dendrogram at the requested number of clusters.
# Assume the Seurat object list was created with SplitObject.
clusters <- cutree(hc, k = 3)

# Add and inspect cluster assignments.
df$patient_type <- clusters

# Add and inspect cluster assignments.
print(df$patient_type)
write.csv(df,"cluster_3_allcohort.csv")



#######xue_cohort##############
xue_cohort=subset(all_tumor_TME,cohort=="c0")
library(ggplot2)
library(dplyr)
library(ggalluvial)

Ratio <- xue_cohort@meta.data %>%
  group_by(Sample, cl2) %>% # Group the observations.
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))

library(reshape2)
data_wide_d<-dcast(Ratio, Sample~Ratio$cl2,
                   value.var = 'relative_freq')
#write.csv(data_wide_d,"cellratio_TME.csv")

library(pheatmap)
n=data_wide_d
rownames(n)=n$Sample
n=n[-1]
n=as.data.frame(t(n))

n=t(scale(t(n)))
n[n>1]=1
n[n<(-1)]=(-1)
color=colorRampPalette(c('#485494','white','#d6252a'))(100)
n[is.na(n)]=0
pheatmap(n,color = color,clustering_method = "ward.D")
#write.csv(n,"heatmap_scale_allcelltype.csv")

annotation_row<-read.csv("heatmap_cell_anno.csv",header=T,row.names = 1)
# Import the annotated dataset.
annotation_row <- as.data.frame(annotation_row)
annotation_col<-read.csv("heatmap_sample_anno.csv",header = T,row.names = 1)
annotation_col <- as.data.frame(annotation_col)
ann_colors = list(
  CD45_type= c(CD45_neg="white", CD45_pos="tomato4"),
  cohort = c(c0=SC_color[21],c1=SC_color[22],c2=SC_color[23],c4=SC_color[24]),
  major_celltype = c(B=SC_color[5],CD4T=SC_color[6],CD8T=SC_color[7],
                     Cycling=SC_color[8],DC=SC_color[9],EC=SC_color[10],
                     Eos=SC_color[11],Fb=SC_color[12],gdT=SC_color[13],
                     ILC=SC_color[14],Mast=SC_color[15],Mph=SC_color[16],
                     Neu=SC_color[17],NK=SC_color[18],PC=SC_color[19])
)

pheatmap(n,color = color,clustering_method = "ward.D",annotation_row = annotation_row,annotation_col =annotation_col,annotation_colors = ann_colors )

# The data frame contains samples as columns and genes as rows.
A=pheatmap(n,color = color,clustering_method = "ward.D",annotation_row = annotation_row,annotation_col =annotation_col,annotation_colors = ann_colors )
df=n
df=as.data.frame(t(df))
# Perform hierarchical clustering.
hc <- A$tree_col

# Plot the clustering dendrogram.
dendrogram <- as.dendrogram(hc)
plot(dendrogram)

# Cut the dendrogram at the requested number of clusters.
# Assume the Seurat object list was created with SplitObject.
clusters <- cutree(hc, k = 3)

# Add and inspect cluster assignments.
df$patient_type <- clusters

# Add and inspect cluster assignments.
print(df$patient_type)
write.csv(df,"cluster_3_xue.csv")

# No significant survival difference was observed.
cli=read.csv(file.path(data_dir, "patient_cli.csv"),header = T,row.names = 1)
cluster=read.csv("cluster_3_cli.csv",header = T,row.names = 1)
cli=cli[match(cluster$patient.1,rownames(cli)),]
out=cbind(cluster,cli)
write.csv(out,"out_cluster_cli.csv")

out=read.csv("out_cluster_cli.csv",header = T,row.names = 1)
library(survival)
library(survminer)
fit <- survfit(Surv(OS_time, OS_state)~cluster, data=out) # Fit the survival model.
fit <- survfit(Surv(FPS_time, Relapse_state)~cluster, data=out) # Fit the survival model.

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
           palette = c("#e2145b", "#0a6b9d","grey")# palette = c("#FF4500", "#4682B4")
)
