
data_dir <- Sys.getenv("HCC_TME_DATA_DIR", unset = "data")

library(qs)
library(Seurat)
# Use HCC_TME_DATA_DIR for external inputs; no local working directory is embedded.
all=qread(file.path(data_dir, "all_subtype_anno2.qs"))
DimPlot(all,group.by = "cl3",raster = T,cols = SC_color)
DimPlot(all,group.by = "Cancer_type",raster = T,cols = SC_color[5:10])
DimPlot(all,group.by = "cohort",raster = T,cols = SC_color[15:30])
DimPlot(all,group.by = "Sample",raster = T,cols = P_color)

########PCA#############
Idents(all)<- all$Sample
av.exp<- AverageExpression(all)$RNA
id=read.csv(file.path(data_dir, "sample.csv"),header = T)
av.exp=av.exp[,match(id$Var1,colnames(av.exp))]
group_list=id$sample_type
# Correct batch effects by cohort.
#library(sva)
#n=ComBat(av.exp,batch = id$cohort)
n=av.exp

## PCA
# BiocManager::install('ggfortify')
library(ggfortify)
exprSet=n
df=as.data.frame(t(exprSet))
df$group=group_list
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')

library("FactoMineR") # Load the required annotation package.
library("factoextra")
df=as.data.frame(t(exprSet))
dat.pca <- PCA(df, graph = FALSE) # Group the observations.
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = group_list, # color by groups
             palette = c("#00AFBB", "#E7B800","#d15c54","#1e93c9"),
             addEllipses = F, # Concentration ellipses
             legend.title = "Groups"
)

########percentage_boxplot_cancertype_majorcell################

meta.tb=all
meta.tb$celltype=meta.tb$cl1
meta.tb$Cancer_type=factor(meta.tb$Cancer_type,levels=c("Healthy","HBV","AL","HCC"))
meta.tb$sampleID = as.character(meta.tb$Sample)
meta.tb_hbv=meta.tb



library(ggplot2)
library(dplyr)
library(ggalluvial)
Ratio <- meta.tb_hbv@meta.data %>%
  group_by(Sample, cl1) %>% # Group cells by sample identifier and cell type.
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


cli=id
cli_neu=cli[match(Ratio$Sample,cli$Var1),]
Ratio_cli=cbind(Ratio,cli_neu)

# load packages
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)

Ratio_cli[1:4,1:4]
Exp=Ratio_cli[c(1,2,4,9)]

library(reshape2)
# Convert data from long to wide format.
data_wide_d<-dcast(Exp, Sample~Exp$cl1,
                   value.var = 'relative_freq')
group=Exp[match(data_wide_d$Sample,Exp$Sample),]
data_wide_d$group=group$sample_type
write.csv(data_wide_d,"cellratio_cancertype_majorcell.csv")









########percentage_boxplot_cancertype_subcell################

meta.tb=all
meta.tb$celltype=meta.tb$cl3
meta.tb$Cancer_type=factor(meta.tb$Cancer_type,levels=c("Healthy","HBV","AL","HCC"))
meta.tb$sampleID = as.character(meta.tb$Sample)
meta.tb_hbv=meta.tb



library(ggplot2)
library(dplyr)
library(ggalluvial)
Ratio <- meta.tb_hbv@meta.data %>%
  group_by(Sample, cl3) %>% # Group cells by sample identifier and cell type.
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


cli=id
cli_neu=cli[match(Ratio$Sample,cli$Var1),]
Ratio_cli=cbind(Ratio,cli_neu)

# load packages
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)

Ratio_cli[1:4,1:4]
Exp=Ratio_cli[c(1,2,4,9)]

library(reshape2)
# Convert data from long to wide format.
data_wide_d<-dcast(Exp, Sample~Exp$cl3,
                   value.var = 'relative_freq')
group=Exp[match(data_wide_d$Sample,Exp$Sample),]
data_wide_d$group=group$sample_type
write.csv(data_wide_d,"cellratio_cancertype_sub1cell.csv")

########percentage_boxplot_cancertype_sub1cell################
Idents(all)=all$cl1
meta.tb=subset(all,idents = c("B_PCs","Myeloid","T_NK"))
meta.tb$celltype=meta.tb$cl3
meta.tb$Cancer_type=factor(meta.tb$Cancer_type,levels=c("Healthy","HBV","AL","HCC"))
meta.tb$sampleID = as.character(meta.tb$Sample)
meta.tb_hbv=meta.tb



library(ggplot2)
library(dplyr)
library(ggalluvial)
Ratio <- meta.tb_hbv@meta.data %>%
  group_by(Sample, cl3) %>% # Group cells by sample identifier and cell type.
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))


cli=id
cli_neu=cli[match(Ratio$Sample,cli$Var1),]
Ratio_cli=cbind(Ratio,cli_neu)

# load packages
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)

Ratio_cli[1:4,1:4]
Exp=Ratio_cli[c(1,2,4,9)]

library(reshape2)
# Convert data from long to wide format.
data_wide_d<-dcast(Exp, Sample~Exp$cl3,
                   value.var = 'relative_freq')
group=Exp[match(data_wide_d$Sample,Exp$Sample),]
data_wide_d$group=group$sample_type
write.csv(data_wide_d,"cellratio_cancertype_sub1cell_remove_stromal.csv")


##############OR_cancertype_immune###########
all=qread(file.path(data_dir, "all_subtype_anno2.qs"))
Idents(all)=all$cl1
immune=subset(all,idents = c("B_PCs","Myeloid","T_NK"))
stromal=subset(all,idents="Stromal")
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
meta.tb=immune
meta.tb$celltype=meta.tb$cl2
meta.tb$Cancer_type=factor(meta.tb$Cancer_type,levels=c("Healthy","HBV","AL","HCC"))
do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          celltype = cellInfo.tb$celltype,
                          colname.patient = "patient",
                          Cancer_type = cellInfo.tb$Cancer_type,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data
  library(data.table)
  dir.create(dirname(out.prefix),F,T)

  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$celltype = as.character(celltype)

  if(is.factor(Cancer_type)){
    cellInfo.tb$Cancer_type = Cancer_type
  }else{cellInfo.tb$Cancer_type = as.factor(Cancer_type)}

  Cancer_type.avai.vec <- levels(cellInfo.tb[["Cancer_type"]])
  count.dist <- unclass(cellInfo.tb[,table(celltype,Cancer_type)])[,Cancer_type.avai.vec]
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
num_hbvstate=table(meta.tb$Cancer_type)
num_cl3=table(meta.tb$cl3)
barplot(num_hbvstate,col = SC_color[c(9,8,7)])
barplot(num_cl3,col = SC_color[5:30],horiz = T)

##############OR_cancertype_stromal###########
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
meta.tb=stromal
meta.tb$celltype=meta.tb$cl2
meta.tb$Cancer_type=factor(meta.tb$Cancer_type,levels=c("Healthy","HBV","AL","HCC"))
do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          celltype = cellInfo.tb$celltype,
                          colname.patient = "patient",
                          Cancer_type = cellInfo.tb$Cancer_type,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data
  library(data.table)
  dir.create(dirname(out.prefix),F,T)

  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$celltype = as.character(celltype)

  if(is.factor(Cancer_type)){
    cellInfo.tb$Cancer_type = Cancer_type
  }else{cellInfo.tb$Cancer_type = as.factor(Cancer_type)}

  Cancer_type.avai.vec <- levels(cellInfo.tb[["Cancer_type"]])
  count.dist <- unclass(cellInfo.tb[,table(celltype,Cancer_type)])[,Cancer_type.avai.vec]
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
num_hbvstate=table(meta.tb$Cancer_type)
num_cl3=table(meta.tb$cl3)
barplot(num_hbvstate,col = SC_color[c(9,8,7)])
barplot(num_cl3,col = SC_color[5:30],horiz = T)


######Roe_stromal############
if(T){
  text.size = 8
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   # axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"),
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   #axis.line = element_line(color = "black"),
                   #axis.ticks = element_line(color = "black"),
                   #panel.grid.minor.y = element_blank(),
                   #panel.grid.minor.x = element_blank(),
                   panel.grid=element_blank(), # Remove grid lines.
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
                   # panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")
                   # strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07",
  )
}

### 1.input data
library("Startrac")
library(data.table)
if (!is.null(out.prefix)) {
  dir.create(dirname(out.prefix),F,T)
}

cellInfo.tb=meta.tb@meta.data
cellInfo.tb = data.table(cellInfo.tb)
cellInfo.tb$cl2 = as.character(cellInfo.tb$cl2)
cellInfo.tb$Cancer_type=cellInfo.tb$Cancer_type %>% as.factor()
Cancer_type.avai.vec<-levels(cellInfo.tb[["Cancer_type"]])


# Calculate observed-to-expected ratios.
startrac.dist <- unclass(Startrac::calTissueDist(cellInfo.tb,byPatient = F,
                                                 colname.cluster="cl2",
                                                 colname.patient = "Sample",
                                                 colname.tissue = "Cancer_type"))
startrac.dist <- startrac.dist[,Cancer_type.avai.vec]
cuts <- c(0,0.0000001, 0.8, 1.5,2,Inf)
startrac.dist.bin.values <- factor(c("-","+/-", "+", "++","+++"),levels=c("-","+/-", "+", "++","+++"))
startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
                            ncol=ncol(startrac.dist))
colnames(startrac.dist.bin) <- colnames(startrac.dist)
rownames(startrac.dist.bin) <- rownames(startrac.dist)

# Visualize the results.
# Visualize the results.
# sscVis::plotMatrix.simple(startrac.dist,
#                           out.prefix=sprintf("%s.startrac.dist",out.prefix),
#                           show.number=F,
#                           clust.row=T,
#                           #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
#                           exp.name=expression(italic(R)[o/e]),
#                           z.hi=2,
#                           #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
#                           palatte=viridis::viridis(7),
#                           pdf.width = 4.5, pdf.height = 4.5)

# Visualize the results.
startrac.dist1 = ifelse(startrac.dist>2,2,startrac.dist)

invisible({
  plot_1 = pheatmap::pheatmap(
    startrac.dist1,
    color = viridis::viridis(7),
    cluster_rows = T,
    cluster_cols = F,
    fontsize = 10,
    border_color = "white",
    treeheight_col = 10,
    treeheight_row = 10
  )
})

# Visualize the results.
plot_data = reshape2::melt(startrac.dist1)
head(plot_data)
colnames(plot_data) = c("cl2", "Cancer_type", "Roe")
startrac.dist.bin2 = reshape2::melt(startrac.dist.bin)
plot_data$bin = startrac.dist.bin2$value %>% factor(levels = c("-","+/-", "+", "++","+++"))
plot_data$cl2 = factor(plot_data$cl2,
                       levels = rev(plot_1$tree_row$labels[plot_1$tree_row$order]))
head(plot_data)

p3 <- ggplot(data = plot_data, aes(x = Cancer_type, y = cl2, fill = Roe)) +
  geom_tile() +
  theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  scale_fill_gradientn(name = "Ro/e",
                       colours = viridis::viridis(7)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  mytheme

p3

colors= c(
  "-" = "white",
  "+/-" = "#fde6ce",
  "+" = "#fcc08b",
  "++" = "#f5904a",
  "+++" = "#e6540d"
)

labels = c("-" = 0, "+/-" = 0.8, "+" = 1.5, "++" = 2,"+++" = ">2")

# Visualize the results.
p4 <- ggplot(data = plot_data, aes(x = Cancer_type, y = cl2, fill = bin)) +
  geom_tile() +
  scale_fill_manual(
    name = "Ro/e",
    values =  colors[unique(plot_data$bin)],
    labels = labels[unique(plot_data$bin)]
  ) +
  geom_text(
    aes(label = bin),
    color = "black",
    size = 3,
    show.legend = T
  ) +
  guides(fill = guide_legend(title = "Ro/e",
                             override.aes = list(label = c(
                               "-","+/-", "+", "++","+++"
                             )))) +
  theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "cm")) + mytheme
p4

######Roe_immune############
if(T){
  text.size = 8
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   # axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"),
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   #axis.line = element_line(color = "black"),
                   #axis.ticks = element_line(color = "black"),
                   #panel.grid.minor.y = element_blank(),
                   #panel.grid.minor.x = element_blank(),
                   panel.grid=element_blank(), # Remove grid lines.
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
                   # panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")
                   # strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07",
  )
}

### 1.input data
library("Startrac")
library(data.table)
if (!is.null(out.prefix)) {
  dir.create(dirname(out.prefix),F,T)
}

meta.tb=immune
cellInfo.tb=meta.tb@meta.data
cellInfo.tb = data.table(cellInfo.tb)
cellInfo.tb$cl2 = as.character(cellInfo.tb$cl2)
cellInfo.tb$Cancer_type=cellInfo.tb$Cancer_type %>% as.factor()
Cancer_type.avai.vec<-levels(cellInfo.tb[["Cancer_type"]])


# Calculate observed-to-expected ratios.
startrac.dist <- unclass(Startrac::calTissueDist(cellInfo.tb,byPatient = F,
                                                 colname.cluster="cl2",
                                                 colname.patient = "Sample",
                                                 colname.tissue = "Cancer_type"))
startrac.dist <- startrac.dist[,Cancer_type.avai.vec]
cuts <- c(0,0.0000001, 0.8, 1.5,2,Inf)
startrac.dist.bin.values <- factor(c("-","+/-", "+", "++","+++"),levels=c("-","+/-", "+", "++","+++"))
startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
                            ncol=ncol(startrac.dist))
colnames(startrac.dist.bin) <- colnames(startrac.dist)
rownames(startrac.dist.bin) <- rownames(startrac.dist)

# Visualize the results.
# Visualize the results.
# sscVis::plotMatrix.simple(startrac.dist,
#                           out.prefix=sprintf("%s.startrac.dist",out.prefix),
#                           show.number=F,
#                           clust.row=T,
#                           #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
#                           exp.name=expression(italic(R)[o/e]),
#                           z.hi=2,
#                           #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
#                           palatte=viridis::viridis(7),
#                           pdf.width = 4.5, pdf.height = 4.5)

# Visualize the results.
startrac.dist1 = ifelse(startrac.dist>2,2,startrac.dist)

invisible({
  plot_1 = pheatmap::pheatmap(
    startrac.dist1,
    color = viridis::viridis(7),
    cluster_rows = T,
    cluster_cols = F,
    fontsize = 10,
    border_color = "white",
    treeheight_col = 10,
    treeheight_row = 10
  )
})

# Visualize the results.
plot_data = reshape2::melt(startrac.dist1)
head(plot_data)
colnames(plot_data) = c("cl2", "Cancer_type", "Roe")
startrac.dist.bin2 = reshape2::melt(startrac.dist.bin)
plot_data$bin = startrac.dist.bin2$value %>% factor(levels = c("-","+/-", "+", "++","+++"))
plot_data$cl2 = factor(plot_data$cl2,
                       levels = rev(plot_1$tree_row$labels[plot_1$tree_row$order]))
head(plot_data)

p3 <- ggplot(data = plot_data, aes(x = Cancer_type, y = cl2, fill = Roe)) +
  geom_tile() +
  theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  scale_fill_gradientn(name = "Ro/e",
                       colours = viridis::viridis(7)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  mytheme

p3

colors= c(
  "-" = "white",
  "+/-" = "#fde6ce",
  "+" = "#fcc08b",
  "++" = "#f5904a",
  "+++" = "#e6540d"
)

labels = c("-" = 0, "+/-" = 0.8, "+" = 1.5, "++" = 2,"+++" = ">2")

# Visualize the results.
p4 <- ggplot(data = plot_data, aes(x = Cancer_type, y = cl2, fill = bin)) +
  geom_tile() +
  scale_fill_manual(
    name = "Ro/e",
    values =  colors[unique(plot_data$bin)],
    labels = labels[unique(plot_data$bin)]
  ) +
  geom_text(
    aes(label = bin),
    color = "black",
    size = 3,
    show.legend = T
  ) +
  guides(fill = guide_legend(title = "Ro/e",
                             override.aes = list(label = c(
                               "-","+/-", "+", "++","+++"
                             )))) +
  theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "cm")) + mytheme
p4


######enrich_mph_hcc_nonhcc#########
mph=subset(all,cl3=="Mph")
table(mph$cl2)
table(mph$Cancer_type)

Idents(mph)=mph$Cancer_type
A=table(mph$Cancer_type)
A=as.data.frame(A)
marker=data.frame(cluster=c(A$Var1),cell=c(A$Var1))
marker$cell=c("non_HCC","non_HCC","HCC","non_HCC")

mph@active.ident=factor(mph@active.ident,levels = c(A$Var1))
mph@meta.data$cl4=sapply(mph@active.ident,function(x){marker[x,2]})
DimPlot(mph,reduction = "umap",group.by = "cl4",label = F,cols = SC_color[10:30],raster =T)

library(Seurat)
library(patchwork)
library(clusterProfiler)
# Load the required annotation package.
library(org.Hs.eg.db) # Load the required annotation package.
library(tidyverse)
# Annotate the data.
#mph <- readRDS("mph.rds")
# Identify marker genes for each cluster.
Idents(mph) <- 'cl4' # Calculate the statistic.
table(mph$cl4)

markers <- FindAllMarkers(mph)
sig_dge.all <- subset(markers, p_val_adj<0.05&abs(avg_log2FC)>0.2) # Retain all significant differentially expressed genes.
write.csv(sig_dge.all,"deg_hcc_nonhcc.csv")
#saveRDS(sig_dge.all, file = "deg.rds")
View(sig_dge.all)
# Select upregulated genes.
sig_dge.M1_up <- subset(sig_dge.all, avg_log2FC>0.5&cluster=='HCC') # Analysis note.
# Run enrichment analysis.
ego_M1 <- enrichGO(gene = sig_dge.M1_up$gene,
                   #universe = row.names(dge.celltype),
                   OrgDb = 'org.Hs.eg.db',
                   keyType = 'SYMBOL',
                   ont = "ALL", # Calculate the statistic.
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)
ego_M1 <- data.frame(ego_M1)
#write.csv(ego_CD4,'enrichGO_CD4.csv')
ego_M1$Group='HCC'
M1=ego_M1[1:5,]


# Select upregulated genes.
sig_dge.M2_up <- subset(sig_dge.all, avg_log2FC>0.5&cluster=='non_HCC') # Analysis note.
# Run enrichment analysis.
ego_M2 <- enrichGO(gene = sig_dge.M2_up$gene,
                   #universe = row.names(dge.celltype),
                   OrgDb = 'org.Hs.eg.db',
                   keyType = 'SYMBOL',
                   ont = "ALL", # Calculate the statistic.
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)
ego_M2 <- data.frame(ego_M2)
#write.csv(ego_CD4,'enrichGO_CD4.csv')
ego_M2$Group='non_HCC'
M2=ego_M2[1:5,]

all <- rbind(M1,M2)
all$Description=factor(all$Description,levels = c(all$Description[!duplicated(all$Description)]))
all$Group=factor(all$Group,levels = c(all$Group[!duplicated(all$Group)]))

ggplot(all, aes(Group, Description)) +
  geom_point(aes(color=-log10(p.adjust), size=GeneRatio))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  #scale_color_gradient(low='#1e93c9',high='#d15c54')+
  scale_color_gradientn(values = seq(0,1,0.2),colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
