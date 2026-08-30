# Use HCC_TME_DATA_DIR for external inputs; no local working directory is embedded.
#####################################################
# Fig2 bh dist
# install.packages("devtools")
# devtools::install_github("arc85/distdimscr")
#library(distdimscr)
data_dir <- Sys.getenv("HCC_TME_DATA_DIR", unset = "data")

source("distdimscr.R")
library(tidyverse)
library(qs)

Sobj<-qread(file.path(data_dir, "all_subtype_anno2.qs"))
Sobj<-subset(Sobj,cl3=="Tumor_Epi",invert=T)

### T vs. AL#####
cellType<-c("B_cell","CD4_T","CD8_T","DC","EC","Fb","Eosinophil","Mast","Mph","Neutrophil","NK","PC","gdT","ILC")

bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(cellType),nrow = 100))
names(bhatt.dist)<-cellType
names(bhatt.dist.rand)<-cellType

for (CT in cellType){

  for (j in 1:10){

    set.seed(j)

    n=length(which(Sobj$Cancer_type=="AL" & Sobj$cl3==CT))
    cells.T<-sample(colnames(Sobj)[which(Sobj$Cancer_type=="HCC" & Sobj$cl3==CT)],n,replace = T)

    cells.Inf <- colnames(Sobj)[Sobj$cl3==CT & Sobj$Cancer_type=="AL"]

    tmp<-Sobj@reductions$harmony@cell.embeddings

    cells.Inf.pca <- tmp[cells.Inf,]
    cells.T.pca <- tmp[cells.T,]

    for (i in 1:10) {

      d<-(j-1)*10+i

      bhatt.dist[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=FALSE)

      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=TRUE)

    }
  }
}
# Analysis note.
library(tidyr)
# Analysis note.
# Analysis note.
data=bhatt.dist.rand
data_long_g<-gather(data, celltype, value, B_cell:ILC)
data_long_g

data=as.data.frame(t(data))
data$median_value=apply(data,1,median)
data=data[order(data$median_value,decreasing = T),]
data_long_g$celltype=factor(data_long_g$celltype,levels = c(rownames(data)))
# Analysis note.
library(ggplot2)
ggplot(data=data_long_g,aes(x=celltype,
                            y=value,
                            color=celltype))+
  geom_jitter(alpha=0.2,
              position=position_jitterdodge(jitter.width = 0.35,
                                            jitter.height = 0,
                                            dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  scale_color_manual(values = SC_color)+
  theme_classic() +
  theme(legend.position="none") +
  theme(text = element_text(size=16)) +
  #ylim(0.0,1.3)+
  ylab("Bhatt distance")+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) # Rotate axis labels by 90 degrees.



### T vs. HBV########
cellType<-c("B_cell","CD4_T","CD8_T","DC","EC","Mph","NK","PC","gdT")

bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(cellType),nrow = 100))
names(bhatt.dist)<-cellType
names(bhatt.dist.rand)<-cellType

for (CT in cellType){

  for (j in 1:10){

    set.seed(j)

    n=length(which(Sobj$Cancer_type=="HBV" & Sobj$cl3==CT))
    cells.T<-sample(colnames(Sobj)[which(Sobj$Cancer_type=="HCC" & Sobj$cl3==CT)],n,replace = T)

    cells.Inf <- colnames(Sobj)[Sobj$cl3==CT & Sobj$Cancer_type=="HBV"]

    tmp<-Sobj@reductions$harmony@cell.embeddings

    cells.Inf.pca <- tmp[cells.Inf,]
    cells.T.pca <- tmp[cells.T,]

    for (i in 1:10) {

      d<-(j-1)*10+i

      bhatt.dist[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=FALSE)

      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=TRUE)

    }
  }
}
# Analysis note.
library(tidyr)
# Analysis note.
# Analysis note.
data=bhatt.dist.rand
data_long_g<-gather(data, celltype, value, B_cell:gdT)
data_long_g

data=as.data.frame(t(data))
data$median_value=apply(data,1,median)
data=data[order(data$median_value,decreasing = T),]
data_long_g$celltype=factor(data_long_g$celltype,levels = c(rownames(data)))
# Analysis note.
library(ggplot2)
ggplot(data=data_long_g,aes(x=celltype,
                            y=value,
                            color=celltype))+
  geom_jitter(alpha=0.2,
              position=position_jitterdodge(jitter.width = 0.35,
                                            jitter.height = 0,
                                            dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  scale_color_manual(values = SC_color)+
  theme_classic() +
  theme(legend.position="none") +
  theme(text = element_text(size=16)) +
  #ylim(0.0,1.3)+
  ylab("Bhatt distance")+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) # Rotate axis labels by 90 degrees.

### T vs. Health#####
cellType<-c("B_cell","CD4_T","CD8_T","EC","Mph","NK","PC","gdT")


bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(cellType),nrow = 100))
names(bhatt.dist)<-cellType
names(bhatt.dist.rand)<-cellType

for (CT in cellType){

  for (j in 1:10){

    set.seed(j)

    n=length(which(Sobj$Cancer_type=="Healthy" & Sobj$cl3==CT))
    cells.T<-sample(colnames(Sobj)[which(Sobj$Cancer_type=="HCC" & Sobj$cl3==CT)],n,replace = T)

    cells.Inf <- colnames(Sobj)[Sobj$cl3==CT & Sobj$Cancer_type=="Healthy"]

    tmp<-Sobj@reductions$harmony@cell.embeddings

    cells.Inf.pca <- tmp[cells.Inf,]
    cells.T.pca <- tmp[cells.T,]

    for (i in 1:10) {

      d<-(j-1)*10+i

      bhatt.dist[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=FALSE)

      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=TRUE)

    }
  }
}
# Analysis note.
library(tidyr)
# Analysis note.
# Analysis note.
data=bhatt.dist.rand
data_long_g<-gather(data, celltype, value, B_cell:gdT)
data_long_g

data=as.data.frame(t(data))
data$median_value=apply(data,1,median)
data=data[order(data$median_value,decreasing = T),]
data_long_g$celltype=factor(data_long_g$celltype,levels = c(rownames(data)))
# Analysis note.
library(ggplot2)
ggplot(data=data_long_g,aes(x=celltype,
                            y=value,
                            color=celltype))+
  geom_jitter(alpha=0.2,
              position=position_jitterdodge(jitter.width = 0.35,
                                            jitter.height = 0,
                                            dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  scale_color_manual(values = SC_color)+
  theme_classic() +
  theme(legend.position="none") +
  theme(text = element_text(size=16)) +
  #ylim(0.0,1.3)+
  ylab("Bhatt distance")+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) # Rotate axis labels by 90 degrees.

#####################################################
# Fig2 bh dist
# install.packages("devtools")
# devtools::install_github("arc85/distdimscr")
library(distdimscr)
library(tidyverse)
library(qs)

Sobj<-qread(file.path(data_dir, "all_downsample.qs"))
Sobj<-subset(Sobj,cl3=="Tumor_Epi",invert=T)

### T vs. AL#####
cellType<-c("B_cell","CD4_T","CD8_T","DC","EC","Fb","Eosinophil","Mast","Mph","Neutrophil","NK","PC","gdT","ILC")

bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(cellType),nrow = 100))
names(bhatt.dist)<-cellType
names(bhatt.dist.rand)<-cellType

for (CT in cellType){

  for (j in 1:10){

    set.seed(j)

    n=length(which(Sobj$Cancer_type=="AL" & Sobj$cl3==CT))
    cells.T<-sample(colnames(Sobj)[which(Sobj$Cancer_type=="HCC" & Sobj$cl3==CT)],n,replace = T)

    cells.Inf <- colnames(Sobj)[Sobj$cl3==CT & Sobj$Cancer_type=="AL"]

    tmp<-Sobj@reductions$harmony@cell.embeddings

    cells.Inf.pca <- tmp[cells.Inf,]
    cells.T.pca <- tmp[cells.T,]

    for (i in 1:10) {

      d<-(j-1)*10+i

      bhatt.dist[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=FALSE)

      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=TRUE)

    }
  }
}
# Analysis note.
library(tidyr)
# Analysis note.
# Analysis note.
data=bhatt.dist.rand
data_long_g<-gather(data, celltype, value, B_cell:ILC)
data_long_g

data=as.data.frame(t(data))
data$median_value=apply(data,1,median)
data=data[order(data$median_value,decreasing = T),]
data_long_g$celltype=factor(data_long_g$celltype,levels = c(rownames(data)))
# Analysis note.
library(ggplot2)
ggplot(data=data_long_g,aes(x=celltype,
                            y=value,
                            color=celltype))+
  geom_jitter(alpha=0.2,
              position=position_jitterdodge(jitter.width = 0.35,
                                            jitter.height = 0,
                                            dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  scale_color_manual(values = SC_color)+
  theme_classic() +
  theme(legend.position="none") +
  theme(text = element_text(size=16)) +
  #ylim(0.0,1.3)+
  ylab("Bhatt distance")+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) # Rotate axis labels by 90 degrees.



### T vs. HBV########
cellType<-c("B_cell","CD4_T","CD8_T","DC","Mph","NK","PC","gdT")

bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(cellType),nrow = 100))
names(bhatt.dist)<-cellType
names(bhatt.dist.rand)<-cellType

for (CT in cellType){

  for (j in 1:10){

    set.seed(j)

    n=length(which(Sobj$Cancer_type=="HBV" & Sobj$cl3==CT))
    cells.T<-sample(colnames(Sobj)[which(Sobj$Cancer_type=="HCC" & Sobj$cl3==CT)],n,replace = T)

    cells.Inf <- colnames(Sobj)[Sobj$cl3==CT & Sobj$Cancer_type=="HBV"]

    tmp<-Sobj@reductions$harmony@cell.embeddings

    cells.Inf.pca <- tmp[cells.Inf,]
    cells.T.pca <- tmp[cells.T,]

    for (i in 1:10) {

      d<-(j-1)*10+i

      bhatt.dist[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=FALSE)

      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=TRUE)

    }
  }
}
# Analysis note.
library(tidyr)
# Analysis note.
# Analysis note.
data=bhatt.dist.rand
data_long_g<-gather(data, celltype, value, B_cell:gdT)
data_long_g

data=as.data.frame(t(data))
data$median_value=apply(data,1,median)
data=data[order(data$median_value,decreasing = T),]
data_long_g$celltype=factor(data_long_g$celltype,levels = c(rownames(data)))
# Analysis note.
library(ggplot2)
ggplot(data=data_long_g,aes(x=celltype,
                            y=value,
                            color=celltype))+
  geom_jitter(alpha=0.2,
              position=position_jitterdodge(jitter.width = 0.35,
                                            jitter.height = 0,
                                            dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  scale_color_manual(values = SC_color)+
  theme_classic() +
  theme(legend.position="none") +
  theme(text = element_text(size=16)) +
  #ylim(0.0,1.3)+
  ylab("Bhatt distance")+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) # Rotate axis labels by 90 degrees.

### T vs. Health#####
cellType<-c("B_cell","CD4_T","CD8_T","Mph","NK","PC","gdT")

bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(cellType),nrow = 100))
names(bhatt.dist)<-cellType
names(bhatt.dist.rand)<-cellType

for (CT in cellType){

  for (j in 1:10){

    set.seed(j)

    n=length(which(Sobj$Cancer_type=="Healthy" & Sobj$cl3==CT))
    cells.T<-sample(colnames(Sobj)[which(Sobj$Cancer_type=="HCC" & Sobj$cl3==CT)],n,replace = T)

    cells.Inf <- colnames(Sobj)[Sobj$cl3==CT & Sobj$Cancer_type=="Healthy"]

    tmp<-Sobj@reductions$harmony@cell.embeddings

    cells.Inf.pca <- tmp[cells.Inf,]
    cells.T.pca <- tmp[cells.T,]

    for (i in 1:10) {

      d<-(j-1)*10+i

      bhatt.dist[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=FALSE)

      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=TRUE)

    }
  }
}
# Analysis note.
library(tidyr)
# Analysis note.
# Analysis note.
data=bhatt.dist.rand
data_long_g<-gather(data, celltype, value, B_cell:ILC)
data_long_g

data=as.data.frame(t(data))
data$median_value=apply(data,1,median)
data=data[order(data$median_value,decreasing = T),]
data_long_g$celltype=factor(data_long_g$celltype,levels = c(rownames(data)))
# Analysis note.
library(ggplot2)
ggplot(data=data_long_g,aes(x=celltype,
                            y=value,
                            color=celltype))+
  geom_jitter(alpha=0.2,
              position=position_jitterdodge(jitter.width = 0.35,
                                            jitter.height = 0,
                                            dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  scale_color_manual(values = SC_color)+
  theme_classic() +
  theme(legend.position="none") +
  theme(text = element_text(size=16)) +
  #ylim(0.0,1.3)+
  ylab("Bhatt distance")+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) # Rotate axis labels by 90 degrees.




#####################################################
#Fig2 roe

ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}




#####################################################
#Fig2 roe

ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}
