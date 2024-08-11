library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(DoubletFinder)
library(limma)
library(Cairo)
library(car)
library(SingleCellExperiment)
library(ggplot2)
library(data.table)

##Load merged single cell Seurat object
load('saveData/1.1.4_seurat_merge_raw.RData')

##normalize
removeBiasGenes = function(genes){
  RPgenes = genes[intersect(grep("^RP", genes), grep("-", genes))]
  RPgenes2 = genes[grep("^RP[SL]", genes)]
  MTgenes = genes[grep("^MT-", genes)]
  CTCgenes = genes[intersect(grep("^CTC", genes), grep("-", genes))]
  MIRgenes = genes[grep("^MIR", genes)]
  ACgenes = genes[intersect(grep("^AC[0-9]", genes), grep(".", genes))]
  CTgenes = genes[intersect(grep("^CT", genes), grep("-", genes))]
  LINCgenes = genes[grep("^LINC[0-9]", genes)]
  ALgenes = genes[intersect(grep("^AL", genes), grep(".", genes))]
  HEMOgene = c("HBB", "HBA1", "HBA2")

  rmgenes = c(RPgenes, RPgenes2, MTgenes, CTCgenes, MIRgenes, ACgenes, CTgenes, LINCgenes, ALgenes, HEMOgene)
  return(rmgenes)
}
rmgenes=removeBiasGenes(rownames(seurat_merge))
seurat_merge = seurat_merge[setdiff(rownames(seurat_merge), rmgenes),]
seurat_merge = NormalizeData(seurat_merge)
seurat_merge = FindVariableFeatures(seurat_merge, selection.method = "vst", nfeatures = 2000)
seurat_merge = ScaleData(seurat_merge)
seurat_merge = RunPCA(seurat_merge)
seurat_merge = RunUMAP(seurat_merge, dims = 1:30)
seurat_merge = RunTSNE(seurat_merge, dims = 1:30)
save(seurat_merge, file='saveData/1.1.4_seurat_merge.RData')

##Choose resolution
seurat_merge = FindNeighbors(seurat_merge, dims = 1:30)
seurat_merge = FindClusters(seurat_merge, resolution = seq(0.4,1.5,0.1))
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_umap_cluster.pdf",width = 7, height = 6)
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.6", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.7", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.9", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.1", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.3", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.4", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.5", label = TRUE))
dev.off()
save(seurat_merge, file = "saveData/1.1.4_seurat_merge.RData")

##Annotation
res4=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE)
res5=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE)
res6=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.6", label = TRUE)
res7=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.7", label = TRUE)
res8=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE)
res9=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.9", label = TRUE)
res10=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1", label = TRUE)
res11=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.1", label = TRUE)
res12=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE)
res13=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.3", label = TRUE)
res14=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.4", label = TRUE)
res15=DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.5", label = TRUE)
ress = plot_grid(res4,res5,res6,res7,res8,res9,res10,res11,res12,res13,res14,res15, ncol = 2)

for (celltype in names(SigGeneral_all_sort)) {
    print(celltype)
    genemarker = SigGeneral_all_sort[[celltype]]
    genemarker = genemarker[genemarker %in% rownames(seurat_merge)]
    width = 18+ceiling(length(genemarker)/3) * 6
    umap_celltype = FeaturePlot(seurat_merge, features = genemarker, reduction = "umap", ncol = ceiling(length(genemarker)/3))
    pdf(file = paste0(workdir, "Results/1.2.1_cluster_anno/1_resolution_selection/",celltype,".pdf"), width=width, height=18)
    print(plot_grid(ress, umap_celltype))
    dev.off()
    tiff(file = paste0(workdir, "Results/1.2.1_cluster_anno/1_resolution_selection/",celltype,".tiff"), res=300, width=width, height=18, compression="lzw", units="in")
    print(plot_grid(ress, umap_celltype, rel_widths = c(18,ceiling(length(genemarker)/3) * 6)))
    dev.off()
}

Idents(seurat_merge) = as.factor(as.numeric(as.character(seurat_merge$RNA_snn_res.1.5)))
annotation = read_excel('rawData/celltype_annotation.xlsx', sheet='celltype-32samples') %>% as.data.frame(stringsAsFactors=F)
celltypes_new = annotation[,2]
names(celltypes_new) = levels(seurat_merge)
seurat_merge = RenameIdents(seurat_merge, celltypes_new)
seurat_merge$cellTypes_new = Idents(seurat_merge)
write.table(seurat_merge@meta.data, file = '1.2.1_seurat_merge_annot_metadata.txt')
save(seurat_merge, file = paste0(workdir, "saveData/1.2.1_seurat_merge_annot.RData"))

clin = read.xlsx('20221130BRCAST_Pathology_Anonymous.xlsx', sheet = 1, check.names = F)
seurat_meta = read.table('1.2.1_seurat_merge_annot_metadata.txt', sep='\t')
seurat_meta_infercnv = read.table('1.4.12_seurat_merge_infercnv_meta.data.txt', sep='\t') #This file contains the label of normal and tumor epithelial cells annotated using Infercnv.
epi = rownames(seurat_meta %>% filter(cellTypes_new=='Epithelial'))
seurat_meta$cellTypes_new2 = seurat_meta$cellTypes_new
seurat_meta[epi, 'cellTypes_new2'] = seurat_meta_infercnv[epi, 'subtypes_latest']
seurat_meta2 = merge(seurat_meta %>% tibble::rownames_to_column('cell_id'), clin, by.x='orig.ident', by.y='sample')
seurat_merge2 = AddMetaData(seurat_merge2, metadata = seurat_meta2)
write.table(seurat_merge2@meta.data, file = '1.2.1_seurat_merge_annot_clin2_meta.data.txt')
saveRDS(seurat_merge2, file='1.2.1_seurat_merge_annot_clin2.rds')

##Plot
sc = readRDS('1.2.1_seurat_merge_annot_clin2.rds')
Idents(sc) = sc$cellTypes_new2
anno_mannual = c('B cells','Endothelial cells','Fibroblasts','Macrophages','Mast cells','Pericytes','Plasmocytes','T cells','Tumor cells','LEC','IPC','MEC')
names(anno_mannual) = c('Bcell','endothelial','fibroblast','macrophage','mastCells','Pericyte','plasmocyte','Tcell','Tumor','LEC','IPC','MEC')
sc = RenameIdents(sc, anno_mannual)

pdf(file = "fig1b.pdf", width = 11.2, height = 8.5)
DimPlot(sc,label = TRUE,cols=c('B cells'='#00C19F', 'Endothelial cells'='#00B9E3', 'Tumor cells'='#F8766D', 'IPC'='#A25C00', 'LEC'='#FFDE59', 'MEC'='#FF00AAFF', 'Fibroblasts'='#93AA00', 'Macrophages'='#00BA38', 'Mast cells'='#FF61C3', 'Pericytes'='#DA6FFC', 'Plasmocytes'='#619CFF', 'T cells'='#D39200'))
dev.off()

#### marker genes of each cell type
data.anno = data.frame(features.plot=c('EPCAM','KRT18',"ERBB2",'CDH1','KRT8',
                                        'PTPRB','VWF','CDH5','LDB2','DOCK9',
                                        'TPSB2','TPSAB1','MS4A2','SLC18A2','HPGDS',
                                        'CD68','CD14','AIF1','APOC1','CD163',
                                        'CD3D','CD3E','CD3G','CD2','CD8A',
                                        'CD79A','CD79B','CD19','MS4A1','BLNK',
                                        'JCHAIN','IGLC2','IGHG3','IGHG1','IGLC1',
                                        'COL1A1','COL1A2','DCN','LUM','ACTA2',
                                        'NOTCH3','MCAM','HIGD1B','ANGPT2',
                                        'AGR2','ESR1','KRT5','KRT14','TP63',"CNN1"),
                        label = c('Tumor cells','Tumor cells','Tumor cells','Tumor cells','Tumor cells',
                                  'Endothelial cells','Endothelial cells','Endothelial cells','Endothelial cells','Endothelial cells',
                                  'Mast cells','Mast cells','Mast cells','Mast cells','Mast cells',
                                  'Macrophages','Macrophages','Macrophages','Macrophages','Macrophages',
                                  'T cells','T cells','T cells','T cells','T cells',
                                  'B cells','B cells','B cells','B cells','B cells',
                                  'Plasmocytes','Plasmocytes','Plasmocytes','Plasmocytes','Plasmocytes',
                                  'Fibroblasts','Fibroblasts','Fibroblasts','Fibroblasts','Fibroblasts',
                                  'Pericytes','Pericytes','Pericytes','Pericytes',
                                  'IPC','IPC','LEC','LEC','MEC','MEC'))
data.usage = DotPlot(sc, features=data.anno$features.plot)$data
df.plot = plyr::join(data.usage,data.anno)
df.plot$id = factor(df.plot$id, levels = sort(levels(df.plot$id),decreasing = T))
print(df.plot)

p = ggplot(df.plot, aes(x=features.plot, y = as.numeric(id), size = pct.exp, color = avg.exp.scaled)) +
  geom_point() + 
  scale_size("% detected", range = c(0,10)) +
  scale_color_gradientn(colours = viridis::viridis(20),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Relative\nexpression") +
  cowplot::theme_cowplot() +
  ylab("") + xlab("") + theme_bw() +
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id),sec.axis = dup_axis()) +
  facet_grid(~label, scales="free_x",space = "free") + theme_classic() +
  theme(
    axis.text.x = element_text(size=18, angle=90, hjust=0.5, color="black",face="bold"),
    axis.text.y = element_text(size=18, color="black",face="bold"),
    #axis.title.x = element_text(size=14,colour = 'black',vjust = -0.8,hjust = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = 'grey30',size = 0.2),
    panel.spacing=unit(0, "mm"),
    strip.text.x = element_text(size=14, face="bold",color = "#FFFFFF",vjust = 0.5,margin = margin(b=3,t=3)),
    strip.background = element_rect(colour="grey30", fill="grey60",size = 1)
  )

pdf('fig1c.pdf',width=23,height=5)
print(p)
dev.off()
