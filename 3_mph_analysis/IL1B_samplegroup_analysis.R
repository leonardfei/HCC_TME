data_dir <- Sys.getenv("HCC_TME_DATA_DIR", unset = "data")

all=qread(file.path(data_dir, "all_subtype_anno2.qs"))
id=read.csv("IL1B_high_rank_sample.csv",header = T,row.names = 1)
all_filt=subset(all,subset=Sample %in% id$Sample)

all_filt=SetIdent(all_filt,value = "Sample")
marker=data.frame(cluster=id$Sample,cell=id$sample_type)
all_filt@active.ident=factor(all_filt@active.ident,levels = c(id$Sample))
all_filt@meta.data$IL1B_type=sapply(all_filt@active.ident,function(x){marker[x,2]})
DimPlot(all_filt,group.by = "cl1",cols = SC_color,raster = T)
DimPlot(all_filt,group.by = "IL1B_type",cols = c("#0a6b9d", "gray"),raster = T)

########percentage_boxplot################
Idents(all_filt)=all_filt$cl3
all_filt_immune=subset(all_filt,idents=c("Tumor_Epi","Stromal_cycling","Fb","EC"),invert=T)
library(ggplot2)
library(dplyr)
library(ggalluvial)
Ratio <- all_filt_immune@meta.data %>%
  dplyr::group_by(Sample, cl3) %>% # Group cells by sample identifier and cell type.
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(relative_freq = n/sum(n))


cli=id
cli_neu=cli[match(Ratio$Sample,cli$Sample),]
Ratio_cli=Ratio
Ratio_cli$IL1B_type=cli_neu$sample_type

# load packages
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)

Ratio_cli[1:4,1:5]
Exp=Ratio_cli[c(1,2,4,5)]

library(reshape2)
# Convert data from long to wide format.
data_wide_d<-dcast(Exp, Sample~Exp$cl3,
                   value.var = 'relative_freq')
group=Exp[match(data_wide_d$Sample,Exp$Sample),]
data_wide_d$group=group$IL1B_type

Exp_plot <- data_wide_d
Exp_plot$group <- factor(Exp_plot$group, levels = c( "IL1B_high", "IL1B_low"))
rownames(Exp_plot)=Exp_plot$Sample
Exp_plot=Exp_plot[-1]

cellname=c(colnames(Exp_plot)[-length(Exp_plot)])
# Set colors for each group.
col <- c("#0a6b9d","gray")
comparisons <- list(c("IL1B_low", "IL1B_high"))
plist <- list() # Create an empty list to store plots generated in the loop.
for (i in 1:length(cellname)) {
  bar_tmp <- Exp_plot[, c(cellname[i], "group")] # Extract expression and sample-group data for the current gene.
  colnames(bar_tmp) <- c("Frequency", "group")
  pb1 <- ggboxplot(bar_tmp, # Create a boxplot using ggboxplot.
                   x = "group", # X-axis is for groups.
                   y = "Frequency", # Y-axis is for expression levels.
                   color = "group", # Fill by sample group.
                   fill = NULL,
                   add = "jitter", # Add jitter points.
                   bxp.errorbar.width = 0.8,
                   width = 0.5,
                   size = 0.1,
                   font.label = list(size = 20),
                   palette = col) +
    theme(panel.background = element_blank())
  pb1 <- pb1 + theme(axis.line = element_line(colour = "black")) +
    theme(axis.title.x = element_blank()) # Adjust the axes.

  pb1 <- pb1 + theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1))

  pb1 <- pb1 + theme(axis.text.y = element_text(size = 15)) +
    ggtitle(cellname[i]) +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))

  pb1 <- pb1 + theme(legend.position = "NA") # Remove the redundant legend.

  pb1 <- pb1 + stat_compare_means(method = "wilcox.test", hide.ns = FALSE,
                                  comparisons = comparisons,label = "p.format",vjust=0.02,bracket.size=0.6) # Run the significance test and add pairwise comparisons.

  plist[[i]] <- pb1 # Store the generated plot.
}
pdf("boxplot_IL1B_group.pdf", height = 12, width = 10)
# Align and arrange the plots into a grid.
plot_grid(plist[[1]], plist[[2]], plist[[3]],
          plist[[4]], plist[[5]], plist[[6]],
          plist[[7]], plist[[8]], plist[[9]],
          plist[[10]], plist[[11]], plist[[12]],
          plist[[13]], plist[[14]], plist[[15]],ncol = 5) # ncol = 4 indicates the number of columns in the grid.
dev.off()

###########Roe_CD8T###################
# Install STARTRAC, sscClust, impute, and sscVis during environment setup.
library(sscVis)
library(data.table)
library(Startrac)
library(readr)
library(dplyr)
library(sscClust)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(Startrac)

# Set the output path and file prefix.
out.prefix = "./Roe"
# Define a custom theme.
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
###########hbv_state
# Preprocess the data.
# Preprocess the data.
cd8t=subset(all_filt_immune,cl3=="CD8_T")
meta.tb=cd8t
cellInfo.tb = meta.tb@meta.data
cellInfo.tb$meta.cluster = as.character(cellInfo.tb$cl2)
cellInfo.tb$loc = cellInfo.tb$IL1B_type %>% as.factor()
loc.avai.vec <- levels(cellInfo.tb[["loc"]])
cellInfo.tb$patient=cellInfo.tb$Sample

# Calculate observed-to-expected ratios.
startrac.dist <- unclass(Startrac::calTissueDist(cellInfo.tb,byPatient = F,
                                                 colname.cluster="meta.cluster",
                                                 colname.patient = "patient",
                                                 colname.tissue = "loc"))
startrac.dist <- startrac.dist[,loc.avai.vec]

# Analysis note.
cuts <- c(0.0000000, 0.0000001, 0.8, 1.5, 2, Inf)
startrac.dist.bin.values <- factor(c("-", "+/-", "+", "++", "+++"),levels=c("-", "+/-", "+", "++", "+++"))
startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
                            ncol=ncol(startrac.dist))

colnames(startrac.dist.bin) <- colnames(startrac.dist)
rownames(startrac.dist.bin) <- rownames(startrac.dist)

# Visualize the results.
# Visualize the results.
sscVis::plotMatrix.simple(startrac.dist,
                          out.prefix=sprintf("%s.startrac.dist",out.prefix),
                          show.number=F,
                          clust.row=T,
                          #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                          exp.name=expression(italic(R)[o/e]),
                          z.hi=2,
                          #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                          palatte=viridis::viridis(7),
                          pdf.width = 4.5, pdf.height = 4.5)

# Visualize the results.
startrac.dist1 = ifelse(startrac.dist>2,2,startrac.dist)

plot_2 = pheatmap::pheatmap(startrac.dist1,
                            color = viridis::viridis(7),
                            cluster_rows = T,
                            cluster_cols = F,
                            fontsize = 10,
                            border_color = "white",
                            treeheight_col = 0,
                            treeheight_row = 0
)

# Visualize the results.
plot_data = reshape2::melt(startrac.dist1)
head(plot_data)
colnames(plot_data) = c("meta.cluster", "loc", "Roe")
startrac.dist.bin2 = reshape2::melt(startrac.dist.bin)
plot_data$bin = startrac.dist.bin2$value %>% factor(levels=c("-", "+/-", "+", "++", "+++"))
plot_data$meta.cluster = factor(plot_data$meta.cluster,
                                levels = rev(plot_2$tree_row$labels[plot_2$tree_row$order]))
head(plot_data)

p3 <- ggplot(data = plot_data, aes(x = loc, y = meta.cluster, fill = Roe)) +
  geom_tile() +
  theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  scale_fill_gradientn(name = "Ro/e",
                       colours = viridis::viridis(7)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  mytheme
p3
#plot_data$loc=factor(plot_data$loc,levels = c("HBV_n","HBV_a","HBV_i"))
#plot_data$meta.cluster=factor(plot_data$meta.cluster,levels = c("Neu_08_NFKBIZ","Neu_07_IFIT1","Neu_06_APOA2","Neu_05_CD83","Neu_04_TXNIP","Neu_03_CD163","Neu_02_S100A12","Neu_01_MMP8"))
# Visualize the results.

p4 <- ggplot(data = plot_data, aes(x = loc, y = meta.cluster, fill = bin)) +
  geom_tile() +
  scale_fill_manual(
    name = "Ro/e",
    values = c(
      # "-" = "white",
      "+/-" = "#f2f3f1",
      "+" = "#b3d3ec",
      "++" = "#3f96c1"
    ),
    labels = c(0.8, 1.5, 2, ">2")
  ) +
  geom_text(
    aes(label = bin),
    color = "black",
    size = 3,
    show.legend = T
  ) +
  guides(fill = guide_legend(title = "Ro/e",
                             override.aes = list(label = c(
                               # "-",
                               "+/-", "+", "++"
                             )))) + theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "cm")) + mytheme;p4
#####TCGA_cibersort_meta_article#######
il1b_type=read.csv("ssgsea_TCGALIHC_IL1B_top50gene.csv",header = T)
meta=read.csv("mmc2 (2).csv",header = T)
meta=meta[meta$TCGA.Study=="LIHC",]
il1b_type=il1b_type[il1b_type$X_PATIENT %in% meta$TCGA.Participant.Barcode,]
il1b_type=il1b_type[substr(il1b_type$X,14,15)<11,]
il1b_type=il1b_type[!(duplicated(il1b_type$X_PATIENT)),]
meta=meta[match(il1b_type$X_PATIENT,meta$TCGA.Participant.Barcode),]
out=cbind(il1b_type,meta)
write.csv(out,"TCGA_LIHC_IL1B_ssgsea_meta.csv")

# Analysis note.
cibersort=read.csv("TCGALIHC_cibersort.csv",header = T,row.names = 1)
cibersort$ID=substr(cibersort$ID,1,15)
il1b_type=il1b_type[il1b_type$X %in% cibersort$ID,]
cibersort=cibersort[match(il1b_type$X,cibersort$ID),]
out2=cbind(il1b_type,cibersort)
write.csv(out2,"TCGA_LIHC_IL1B_ssgsea_cibersort_meta.csv")

###########Roe_allcell###################
# Install STARTRAC, sscClust, impute, and sscVis during environment setup.
library(sscVis)
library(data.table)
library(Startrac)
library(readr)
library(dplyr)
library(sscClust)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(Startrac)

# Set the output path and file prefix.
out.prefix = "./Roe"
# Define a custom theme.
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
###########hbv_state
# Preprocess the data.
# Preprocess the data.
meta.tb=all_filt
cellInfo.tb = meta.tb@meta.data
cellInfo.tb$meta.cluster = as.character(cellInfo.tb$cl2)
cellInfo.tb$loc = cellInfo.tb$IL1B_type %>% as.factor()
loc.avai.vec <- levels(cellInfo.tb[["loc"]])
cellInfo.tb$patient=cellInfo.tb$Sample

# Calculate observed-to-expected ratios.
startrac.dist <- unclass(Startrac::calTissueDist(cellInfo.tb,byPatient = F,
                                                 colname.cluster="meta.cluster",
                                                 colname.patient = "patient",
                                                 colname.tissue = "loc"))
startrac.dist <- startrac.dist[,loc.avai.vec]

# Analysis note.
cuts <- c(0.0000000, 0.0000001, 0.8, 1.5, 2, Inf)
startrac.dist.bin.values <- factor(c("-", "+/-", "+", "++", "+++"),levels=c("-", "+/-", "+", "++", "+++"))
startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
                            ncol=ncol(startrac.dist))

colnames(startrac.dist.bin) <- colnames(startrac.dist)
rownames(startrac.dist.bin) <- rownames(startrac.dist)

# Visualize the results.
# Visualize the results.
sscVis::plotMatrix.simple(startrac.dist,
                          out.prefix=sprintf("%s.startrac.dist",out.prefix),
                          show.number=F,
                          clust.row=T,
                          #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                          exp.name=expression(italic(R)[o/e]),
                          z.hi=2,
                          #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                          palatte=viridis::viridis(7),
                          pdf.width = 4.5, pdf.height = 4.5)

# Visualize the results.
startrac.dist1 = ifelse(startrac.dist>2,2,startrac.dist)

plot_2 = pheatmap::pheatmap(startrac.dist1,
                            color = viridis::viridis(7),
                            cluster_rows = T,
                            cluster_cols = F,
                            fontsize = 10,
                            border_color = "white",
                            treeheight_col = 0,
                            treeheight_row = 0
)

# Visualize the results.
plot_data = reshape2::melt(startrac.dist1)
head(plot_data)
colnames(plot_data) = c("meta.cluster", "loc", "Roe")
startrac.dist.bin2 = reshape2::melt(startrac.dist.bin)
plot_data$bin = startrac.dist.bin2$value %>% factor(levels=c("-", "+/-", "+", "++", "+++"))
plot_data$meta.cluster = factor(plot_data$meta.cluster,
                                levels = rev(plot_2$tree_row$labels[plot_2$tree_row$order]))
head(plot_data)

p3 <- ggplot(data = plot_data, aes(x = loc, y = meta.cluster, fill = Roe)) +
  geom_tile() +
  theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  scale_fill_gradientn(name = "Ro/e",
                       colours = viridis::viridis(7)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  mytheme
p3
#plot_data$loc=factor(plot_data$loc,levels = c("HBV_n","HBV_a","HBV_i"))
#plot_data$meta.cluster=factor(plot_data$meta.cluster,levels = c("Neu_08_NFKBIZ","Neu_07_IFIT1","Neu_06_APOA2","Neu_05_CD83","Neu_04_TXNIP","Neu_03_CD163","Neu_02_S100A12","Neu_01_MMP8"))
# Visualize the results.

p4 <- ggplot(data = plot_data, aes(x = loc, y = meta.cluster, fill = bin)) +
  geom_tile() +
  scale_fill_manual(
    name = "Ro/e",
    values = c(
      # "-" = "white",
      "+/-" = "#f2f3f1",
      "+" = "#b3d3ec",
      "++" = "#3f96c1",
      "+++" = "#104c8b"
    ),
    labels = c(0.8, 1.5, 2, ">2")
  ) +
  geom_text(
    aes(label = bin),
    color = "black",
    size = 3,
    show.legend = T
  ) +
  guides(fill = guide_legend(title = "Ro/e",
                             override.aes = list(label = c(
                               # "-",
                               "+/-", "+", "++","+++"
                             )))) + theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "cm")) + mytheme;p4
