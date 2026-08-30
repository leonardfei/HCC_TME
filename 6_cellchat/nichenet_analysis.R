# install.packages("devtools")
data_dir <- Sys.getenv("HCC_TME_DATA_DIR", unset = "data")

# Install nichenetr during environment setup; see the repository README.

library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
# Use HCC_TME_DATA_DIR for external inputs; no local working directory is embedded.
Idents(mph2_filt2)=mph2_filt2$cl3
mph2_filt2$celltype=mph2_filt2$cl3
tumor=qread(file.path(data_dir, "tumor_hcc_noepi.qs"))
tumor$celltype=tumor$cl2
Idents(tumor)=tumor$celltype

Ca_filt=subset(tumor,idents=c("Ca_BICC1","Ca_TOP2A"))
mph_ca=merge(Ca_filt,mph2_filt2)


#fb_03_mph_filt_downsample=subset(fb_03_mph_filt,downsample=2000)
seuratObj <- mph_ca
seuratObj@meta.data %>% head()
seuratObj@meta.data$celltype %>% table() # note that the number of cells of some cell types is very low and should preferably be higher for a real application
seuratObj@meta.data$Cancer_type %>% table()


lr_network = read_rds(file.path(data_dir, "lr_network_human_21122021.rds"))
ligand_target_matrix = read_rds(file.path(data_dir, "ligand_target_matrix_nsga2r_final.rds"))
weighted_networks = read_rds(file.path(data_dir, "weighted_networks_nsga2r_final.rds"))

lr_network = lr_network %>% distinct(from, to)
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

## receiver
Idents(seuratObj) <- 'celltype'
receiver = "Ca_BICC1"
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c( "Mono_like_Mph", "Mph_CXCL9", "Mph_MARCO","Mph_IGF1","Mph_SPP1")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
##############################HCC###############################################
# Define the affected gene set.
#seurat_obj_receiver= subset(seuratObj, idents = receiver)
#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["celltype"]])

condition_oi = "Ca_TOP2A"
condition_reference = "Ca_BICC1"

DE_table_receiver = FindMarkers(object = tumor, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

# Define ligand-receptor candidates.
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()


# Analyze ligand activity.
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

#The different ligand activity measures (auroc, aupr, pearson correlation coefficient) are a measure for how well a ligand can predict the observed differentially expressed genes compared to the background of expressed genes.
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))

best_upstream_ligands = ligand_activities %>% top_n(30, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
DotPlot(seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()


# Identify receptors and target genes.
#Active target gene inference
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
vis_ligand_target=vis_ligand_target[,1:40]
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

# Identify ligand-receptor pairs.
#Receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

# Calculate ligand log-fold changes.
## DE analysis for each sender cell type
## this uses a new nichenetr function - reinstall nichenetr if necessary!
#DE_table_all = Idents(seuratObj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seuratObj, condition_colname = "celltype", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
#DE_table_all[is.na(DE_table_all)] = 0
#
## Combine ligand activities with DE information
#ligand_activities_de = ligand_activities %>% dplyr::select(test_ligand, pearson) %>% dplyr::rename(ligand = test_ligand) %>% left_join(DE_table_all %>% dplyr::rename(ligand = gene))
#ligand_activities_de[is.na(ligand_activities_de)] = 0
#
## make LFC heatmap
#lfc_matrix = ligand_activities_de  %>% dplyr::select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
#rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
#
#order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
#vis_ligand_lfc = lfc_matrix[order_ligands,]
#
#colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
#
#p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
#p_ligand_lfc


###########SUMMARIZE
# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))

# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
rotated_dotplot = DotPlot(seuratObj %>% subset(celltype %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots

figures_without_legend = cowplot::plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  #p_ligand_lfc +
  theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_target)))

legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot

markers=order_targets[1:40]

library(ggplot2)
DotPlot(tumor, features = markers,group.by = "cl2")+coord_flip()+theme_bw()+ # Use a clean theme and rotate the plot.
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+ # Rotate axis labels by 90 degrees.
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#f6fbef","#b2e1ba","#49add0","#094687"))+ # Set the color gradient.
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) +scale_size(range = c(0,4))
