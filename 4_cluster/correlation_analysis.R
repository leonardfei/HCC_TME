Idents(all)=all$cl1
all_TME=subset(all,idents=c("Malignant_Epi"),invert=T)
all_hcc=subset(all_TME,Cancer_type=="HCC")
all_al=subset(all_TME,Cancer_type=="AL")

Idents(all_hcc)=all_hcc$cl2
all_hcc=subset(all_hcc,idents = c("EC_T_CD3D","EC_Epi_KRT19","Fb_Epi_KRT19","Stromal_cycling","Myeloid_cycling","B_PC_cycling","T_NK_cycling","Mph_CD68_CD247"),invert=T)
Idents(all_al)=all_al$cl2
all_al=subset(all_al,idents = c("EC_T_CD3D","Stromal_cycling","Myeloid_cycling","B_PC_cycling","T_NK_cycling","Mph_CD68_CD247"),invert=T)

# Analysis note.
meta=all_hcc@meta.data
meta_filt=meta[c(4,28)]
meta_filt$all=paste(meta_filt$Sample,sep = "/",meta_filt$cl2)
meta_filt2=meta_filt[!duplicated(meta_filt$all),]
table(meta_filt2$cl2)
##########cell ratio hcc################################################
library(ggplot2)
library(dplyr)
library(ggalluvial)

# Analysis note.
Ratio <- all_hcc@meta.data %>%
  dplyr::group_by(Sample, cl2) %>% # Group the observations.
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(relative_freq = n/sum(n))

library(reshape2)
data_wide_d<-reshape2::dcast(Ratio, Sample~Ratio$cl2,
                             value.var = 'relative_freq')
data_wide_d[is.na(data_wide_d)]=0
library(ggcorrplot)
library(corrplot)
rownames(data_wide_d)=data_wide_d$Sample
data_wide_d=data_wide_d[-1]
data_wide_d=scale(data_wide_d)

M <- cor(data_wide_d,method = "pearson")
#M=M[,-c(10,17,29,36)]
#M=M[-c(10,17,29,36),]
#color=colorRampPalette(c('#485494','white','#d7232a'))(100)
#color=colorRampPalette(c('#2166ac','white','#b2182b'))(100)
# Analysis note.
library(pheatmap)
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
color=c(colorRampPalette(colors = c("#2166ac","white"))(length(bk)/2),colorRampPalette(colors = c("white","#b2182b"))(length(bk)/2))
pheatmap(M,color = color,breaks = bk,clustering_method = "ward.D2",border_color = NA)
fig1=pheatmap(M,color = color,breaks = bk,clustering_method = "ward.D2",border_color = NA)

id=fig1$tree_col$labels
id=as.data.frame(id)
id$num=1:49
id=id[match(fig1$tree_col$order,id$num),]

##########cell ratio AL################################################
library(ggplot2)
library(dplyr)
library(ggalluvial)

# Analysis note.
Ratio <- all_al@meta.data %>%
  dplyr::group_by(Sample, cl2) %>% # Group the observations.
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(relative_freq = n/sum(n))

library(reshape2)
data_wide_d<-reshape2::dcast(Ratio, Sample~Ratio$cl2,
                             value.var = 'relative_freq')
data_wide_d[is.na(data_wide_d)]=0
library(ggcorrplot)
library(corrplot)
rownames(data_wide_d)=data_wide_d$Sample
data_wide_d=data_wide_d[-1]
data_wide_d=scale(data_wide_d)

M <- cor(data_wide_d,method = "pearson")
M=M[match(id$id,rownames(M)),match(id$id,colnames(M))]
#M=M[,-c(10,17,29,36)]
#M=M[-c(10,17,29,36),]
#color=colorRampPalette(c('#485494','white','#d7232a'))(100)
#color=colorRampPalette(c('#2166ac','white','#b2182b'))(100)
# Analysis note.
library(pheatmap)
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
color=c(colorRampPalette(colors = c("#2166ac","white"))(length(bk)/2),colorRampPalette(colors = c("white","#b2182b"))(length(bk)/2))
pheatmap(M,color = color,breaks = bk,
         #clustering_method = "ward.D2",
         cluster_rows = F,cluster_cols = F,
         border_color = NA)
