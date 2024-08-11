library(Seurat)
library(dplyr)
library(ggplot2)
library(MuDataSeurat)
library(infercnv)
options(bitmapType = "cairo")

args = commandArgs(T)
sample = args[1]
st_ids = c('P13'='148BR_D1_13-1', 'P15'='148BR_F2_15-1', 'P161'='148TR_E1_16-1', 'P22'='22-1', 'P28'='28-1', 'P29'='29-1', 'P35'='35-1', 'P36'='36-1', 'P41'='41-1')
st_id = st_ids[sample]
#indir = paste0(st_id, '/Results/8.23_cellbin_merge_cluster_bin100_merge2rds')
indir = paste0(st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds')
indir2 = paste0(st_id, '/Results/8.30_cellbin_bin100_tree')
indir3 = paste0(st_id, '/Results/8.28.2_cellbin_bin100_analyze_infercnv_selectFeatures')
workdir=paste0(st_id, '/Results/8.33_markers_expr')
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

##读取数据
#seurat_cellbin_bin100 = readRDS(paste0(indir, '/adata_merge.rds'))
seurat_cellbin_bin100 = ReadH5AD(paste0(indir, '/merged_adata.h5ad'))
subclone = read.delim(paste0(indir2, '/8.30_labels_subclone.txt'), sep='\t')

##添加subclone label
seurat_cellbin_bin100$subclone = NA
seurat_cellbin_bin100@meta.data[rownames(subclone), 'subclone'] = subclone[,1]
seurat_cellbin_epi = seurat_cellbin_bin100[,rownames(subclone)]

ref_group_names_list = list('P13'=c('B cells', 'Endothelial cells', 'Fibroblasts', 'MEC', 'Macrophages', 'Mast cells', 'Normal epithelial cells', 'Pericytes', 'Plasmocytes', 'T cells'),
                            'P161'=c("Epithelial_normal", 'apCAF', 'CAF_MGP+', 'iCAF', 'myCAF_1', 'myCAF_2', 'myCAF_3', 'Plasmocyte_1', 'Plasmocyte_2', 'Plasmocyte_3', 'Plasmocyte_4'),
                            'P15'=c('Epithelial_normal'),
                            'P22'=c('B cells', 'Endothelial cells', 'Fibroblasts', 'MEC', 'Macrophages', 'Mast cells', 'Normal epithelial cells', 'Pericytes', 'Plasmocytes', 'T cells'),
                            'P28'=c('B cells', 'Endothelial cells', 'Fibroblasts', 'MEC', 'Macrophages', 'Mast cells', 'Normal epithelial cells', 'Pericytes', 'Plasmocytes', 'T cells'),
                            'P29'=c('Fibroblast_1', 'Plasmocyte_1', 'Plasmocyte_2'),
                            'P35'=c('B cells', 'Endothelial cells', 'Fibroblasts', 'MEC', 'Macrophages', 'Mast cells', 'Normal epithelial cells', 'Pericytes', 'Plasmocytes', 'T cells'),
                            'P36'=c('Epithelial_8', 'Epithelial_normal'),
                            'P41'=c('B cells', 'Endothelial cells', 'Fibroblasts', 'MEC', 'Macrophages', 'Mast cells', 'Normal epithelial cells', 'Pericytes', 'Plasmocytes', 'T cells'))
ref_group_names = ref_group_names_list[[sample]]
seurat_cellbin_ref = seurat_cellbin_bin100[,rownames(seurat_cellbin_bin100@meta.data)[seurat_cellbin_bin100$annotated_cluster %in% ref_group_names]]

#violin plot of marker genes expression(Fig3p,q)
markers = c('ESR1', "PGR", 'AR', 'ERBB2', 'MKI67')
markers %in% rownames(seurat_cellbin_bin100)
seurat_cellbin_merge = merge(seurat_cellbin_ref, seurat_cellbin_epi)
seurat_cellbin_merge$subclone[is.na(seurat_cellbin_merge$subclone)] = 'Reference'
#saveRDS(seurat_cellbin_merge, file=paste0(workdir, '/8.33_seurat_cellbin_merge.rds'))
vlnplot = VlnPlot(seurat_cellbin_merge, features = markers, group.by = 'subclone', pt.size = 0, stack=T, flip = T) +
  NoLegend() +
  labs(x=NULL, y='Relative Expression')
pdf(file = paste0(workdir, '/8.33_vlnplot.pdf'), width = 18, height = 9)
print(VlnPlot(seurat_cellbin_merge, features = markers, group.by = 'subclone', pt.size = 0, stack=T, flip = T) +
  NoLegend() +
  labs(x=NULL, y='Relative Expression') +
  theme(axis.text.x=element_text(angle=0)),
        text=element_text(size=7))
dev.off()

genesite <- yaGST::gmt2GO("c1.all.v7.5.1.symbols.gmt")
pq_genes = data.frame()
for (pq in setdiff(names(genesite), 'MT')) {
    pq_i = data.frame(region=stringr::str_extract(pq,"chr[0-9A-Z]*[pq]"), gene=genesite[[pq]])
    pq_genes = rbind(pq_genes, pq_i)
}
markers_pq = pq_genes %>% filter(gene %in% markers)
rownames(markers_pq) = markers_pq$gene
markers_pq = markers_pq[markers,]

#CNV score for ERBB2 gene of each subclone(Fig3r)
infercnv_obj = readRDS(paste0(indir3, '/infercnv_obj2.rds'))
markers2 = c('ESR1', "PGR", 'ERBB2', 'MKI67', 'EGFR', 'KRT5')
markers2 %in% rownames(infercnv_obj@expr.data)
markers2 = intersect(markers2, rownames(infercnv_obj@expr.data))
seurat_cellbin_her2 = seurat_cellbin_epi[markers2, colnames(infercnv_obj@expr.data)]
#seurat_cellbin_epi@assays$Spatial@data[markers2,] = infercnv_obj@expr.data[markers2, colnames(seurat_cellbin_epi@assays$Spatial@data)]
seurat_cellbin_epi@assays$RNA@data[markers2,] = infercnv_obj@expr.data[markers2, colnames(seurat_cellbin_epi@assays$RNA@data)]
vlnplot = VlnPlot(seurat_cellbin_her2, features = markers2, group.by = 'subclone', pt.size = 0, stack=T, flip = T) +
  NoLegend() +
  labs(x=NULL, y='CNV Score')
vlnplot$data = vlnplot$data %>% filter(feature=='ERBB2')
pdf(file = paste0(workdir, '/8.33_vlnplot_her2.pdf'), width = 5, height = 8)
print(vlnplot)
dev.off()
