library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
library(magrittr)
library('MuDataSeurat')
options(stringsAsFactors = FALSE, future.globals.maxSize=20000*1024^2)

od <- '3.3_mergeCellchat_TIL/'
if(!dir.exists(od)){dir.create(od, recursive = TRUE)}
setwd(od)

sc <- MuDataSeurat::ReadH5AD('merged_High_TIL_st_raw.h5ad')
str(sc)
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)

##cell communications in high/low immune infiltration groups(for Fig6j,k)
sc$annotated_cluster < as.factor(sc$annotated_cluster)
table(sc$annotated_cluster)
data.input  <- sc@assays$RNA@data
identity = data.frame(group=sc$annotated_cluster, row.names = names(sc$annotated_cluster)) # create a dataframe consisting of the cell labels
head(identity)
cellchat <- createCellChat(data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
head(cellchat@meta)
str(cellchat)

cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
table(cellchat@idents)
save(cellchat, file="cellchat.RData")

CellChatDB <- readRDS('CellChatDB.human_nichenet.rds')
CellChatDB.use <- CellChatDB
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)
save(cellchat, file = "cellchat_allDB.RData")

## pre-processing
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
save(cellchat, file = "cellchat_preprocess.RData")

## interaction inference
cellchat <- computeCommunProb(cellchat, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)
save(cellchat, file = "cellchat_interaction.RData")

####### after generating the cellchat object of HIGH-TIL and LOW-TIL separately, we merge them and do comparision ########

x <- load('cellchat_interaction.RData')
cellchat_high_TIL <- get(x)
y <- load('cellchat_interaction.RData')
cellchat_low_TIL <- get(y)

cellchat <- mergeCellChat(list(cellchat_high_TIL, cellchat_low_TIL), add.names = c("high_TIL", "low_TIL"), cell.prefix = TRUE)
pdf(paste0(od, "compare_number.pdf"),width=3,height=3)
compareInteractions(cellchat, show.legend = F, group = c(1,2))
compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
dev.off()

pdf(paste0(od, "compare_number_bycell.pdf"),width=6,height=4)
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

pdf(paste0(od, "rankNET.pdf"),width=4,height=7)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
dev.off()

library(ComplexHeatmap)
cellchat_high_TIL <- netAnalysis_computeCentrality(cellchat_high_TIL, slot.name = "netP")
cellchat_low_TIL <- netAnalysis_computeCentrality(cellchat_low_TIL, slot.name = "netP")
object.list <- list(high_TIL = cellchat_high_TIL, low_TIL = cellchat_low_TIL)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16)
pdf(paste0(od, "compare_outgoing_signal_bycell.pdf"),width=8,height=13)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
dev.off()
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16, color.heatmap = "GnBu")
pdf(paste0(od, "compare_incoming_signal_bycell.pdf"),width=8,height=13)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
dev.off()
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16, color.heatmap = "OrRd")
pdf(paste0(od, "compare_all_signal_bycell.pdf"),width=8,height=13)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
dev.off()

##### Compare L-R pairs
#pdf(paste0(od, "L-R_InvEpi_Immune_compare.pdf"),width=6,height=12)
#netVisual_bubble(cellchat, sources.use = c('Epithelial_Invasive_cancer','myCAF'), targets.use = c('B_cell','Macrophage','Plasmocyte','T_cell'),  comparison = c(1, 2), angle.x = 45)
#dev.off()

pdf(paste0(od, "IncreaseorDecrease_IBC_fib_signal.pdf"),width=15,height=11)
gg1 <- netVisual_bubble(cellchat, sources.use = c('IBC cells','Fibroblasts'), targets.use = c('IBC cells','Fibroblasts'),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in low_TIL samples", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c('IBC cells','Fibroblasts'), targets.use = c('IBC cells','Fibroblasts'),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in low_TIL samples", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "high_TIL"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "high_TIL",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "low_TIL",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
print(pairLR.use.up)
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('Fibroblasts','IBC cells'), targets.use = c('Macrophages','T cells','B cells','Plasmocytes'),
                        comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
#gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('IBC cells','Fibroblasts'), targets.use = c('Fibroblasts'), 
                        #comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
pdf(paste0(od, "Up-regulated_IBC&fib_sender_signal.pdf"),width=7,height=33)
gg1
dev.off()

##selected LR
pairLR.use.up.selected = pairLR.use.up %>% filter(interaction_name %in% c('CXCL12_CXCR4', 'CXCL14_CXCR4', 'CXCL2_ACKR1', 'CXCL9_ACKR1', 'CXCL10_ACKR1', 'CXCL12_ITGB1', 'MFGE8_PDGFRB', 'SCGB3A1_NOTCH3', 'GRN_SORT1', 'MDK_NCL', 'MDK_LRP1', 'MDK_LRP2', 'MDK_SDC1'))
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up.selected, sources.use = c('Fibroblasts','IBC cells'), targets.use = c('Macrophages','T cells','B cells','Plasmocytes'),
                        comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", "high_TIL"))
pdf(paste0(od, "Up-regulated_IBC&fib_sender_signal_selected.pdf"),width=6,height=4)
gg1
dev.off()


gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('Fibroblasts','IBC cells'), targets.use = c('Macrophages','T cells','B cells','Plasmocytes'), 
                        comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
pdf(paste0(od, "Down-regulated_IBC&fib_sender_signal.pdf"),width=7,height=30)
gg2
dev.off()

gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('IBC cells'), targets.use = c('IBC cells','Endothelial cells','Fibroblasts'), 
                        comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
pdf(paste0(od, "Up-regulated_IBC_sender_signal.pdf"),width=5,height=28)
gg2
dev.off()

gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('IBC cells'), targets.use = c('IBC cells','Endothelial cells','Fibroblasts'), 
                        comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
pdf(paste0(od, "Down-regulated_IBC_sender_signal.pdf"),width=5,height=34)
gg2
dev.off()

##### plot gene expression
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("high_TIL", "low_TIL")) # set factor level
pdf(paste0(od, "GeneExpression.pdf"),width=13,height=100)
plotGeneExpression(cellchat, signaling = "Other", split.by = "datasets", colors.ggplot = T)
dev.off()