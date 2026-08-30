library(Seurat)
library(qs)

tumor=qread("tumor_hcc_noepi.qs")
mph=qread("../4_mph_analysis/ucell_mph2_filt_final2.qs")
mph$cl3=rep("mph",ncol(mph))
tumor$cl3=tumor$cl2
all=merge(tumor,mph)
set.seed(123) # Set a seed for reproducibility.
cells.use <- sample(
  colnames(all),
  size = 50000,
  replace = FALSE
)
seurat_sub <- subset(all, cells = cells.use)
seurat_sub$samples=seurat_sub$Sample

# Install CellChat during environment setup; see the repository README.
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
# Here we load a scRNA-seq data matrix and its associated cell meta data
#load(url("https://ndownloader.figshare.com/files/25950872")) # This is a combined data from two biological conditions: normal and diseasescc
############cellchat_hbv_a############
data.input = seurat_sub@assays$RNA@data # normalized data matrix

meta = seurat_sub@meta.data # a dataframe with rownames containing cell mata data
#cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data

# Prepare input data for CelChat analysis
#data.input = data.input[, cell.use]
#meta = meta[cell.use, ]
# meta = data.frame(cl3 = meta$cl3[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell cl3
unique(meta$cl3)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cl3")
#> Create a CellChat object from a data matrix
#> Set cell identities for the new CellChat object
#> The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 LC Inflam. DC TC Inflam. TC CD40LG+ TC NKT

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cl3") # set "cl3" as default cell identity
levels(cellchat@idents) # show factor levels of the cell cl3
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
#cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#test=subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)
df.net <- subsetCommunication(cellchat)
#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
#gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions",sources.use = "mph")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength",sources.use = "mph")

A=cellchat@netP$pathways

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pdf("./cxcl.signalingrole.pdf",width = 8,height = 3.5)
netAnalysis_signalingRole_network(cellchat, signaling = A, width = 8, height = 2.5, font.size = 10)
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
pdf("./cxcl.signalingrole2d.pdf",width = 8,height = 2.5)
gg1 + gg2
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf("./outin.cell.pdf",width = 8,height = 6)
ht1 + ht2
dev.off()

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, targets.use = c("Neutrophil"), remove.isolate = FALSE)
#> Comparing communications on a single object

pathways.show <- c("CXCL")
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells
vertex.receiver = seq(4,7) # a numeric vector.
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
qsave(cellchat,"cellchat.qs")
