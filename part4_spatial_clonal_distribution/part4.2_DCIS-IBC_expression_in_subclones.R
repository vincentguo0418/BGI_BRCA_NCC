library('dplyr')
library('limma')
library('MuDataSeurat')
library('Seurat')
library('infercnv')

args = commandArgs(T)
sample = args[1]
st_ids = c('P13'='148BR_D1_13-1', 'P15'='148BR_F2_15-1', 'P161'='148TR_E1_16-1', 'P22'='22-1', 'P28'='28-1', 'P29'='29-1', 'P35'='35-1', 'P36'='36-1', 'P41'='41-1')
samples = names(st_ids)
workdir=paste0('/paper/ST_pca_subclone')
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

##UMAP of DCIS and IBC expression(Fig4c,f,i,l)
for (sample in samples) {
    st_id = st_ids[sample]
    seurat_cellbin_bin100 = ReadH5AD(paste0('/ST/', st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds/merged_adata.h5ad'))
    seurat_i = subset(seurat_cellbin_bin100, subset=annotated_cluster %in% c('DCIS cells', 'IBC cells'))
    seurat_i = NormalizeData(seurat_i)
    seurat_i <- FindVariableFeatures(seurat_i, selection.method = "vst", nfeatures = 2000)
    seurat_i <- ScaleData(seurat_i)
    seurat_i <- RunPCA(seurat_i)
    seurat_i <- RunUMAP(seurat_i, dims = 1:30)

    labels_hclust = read.delim(paste0('/ST/', st_id, '/Results/8.30_cellbin_bin100_tree/8.30_labels_subclone.txt'), sep = '\t')
    cells = intersect(rownames(labels_hclust), colnames(seurat_i))
    seurat_i = seurat_i[,cells]
    seurat_i$subclone = labels_hclust[rownames(seurat_i@meta.data),]
    table(seurat_i$subclone, as.character(seurat_i$annotated_cluster))

    pdf(paste0(workdir, '/', sample, '_umap_expr.pdf'),width=7,height=7)
    print(DimPlot(seurat_i, reduction = "umap", pt.size = 1.2, group.by = "subclone", shape.by='annotated_cluster', label = TRUE) + ggforce::geom_mark_hull())
    dev.off()
    saveRDS(seurat_i, file=paste0(workdir, '/', sample, '.rds'))
}
