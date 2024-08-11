library('Seurat')
library('dplyr')
library('limma')
library('MuDataSeurat')
library('car')

args = commandArgs(T)
sample = args[1]
st_ids = c('P13'='148BR_D1_13-1', 'P15'='148BR_F2_15-1', 'P161'='148TR_E1_16-1', 'P22'='22-1', 'P28'='28-1', 'P29'='29-1', 'P35'='35-1', 'P36'='36-1', 'P41'='41-1')
samples = names(st_ids)
workdir=paste0('/paper/ST_subclone_region_nhood')
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

##Spatial coherence between spatial CNV subclones(Fig4m)
#*P13
sample = 'P13'
st_id = st_ids[sample]
seurat_cellbin = ReadH5AD(paste0('/ST/', st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds/P13_raw_annotated.h5ad'))
cells_map = read.delim(paste0('/ST/', st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds/merged_cluster_labels.txt'), colClasses = "character")
labels_hclust = read.delim(paste0('/ST/', st_id, '/Results/8.30_cellbin_bin100_tree/8.30_labels_subclone.txt'), sep = '\t')

seurat_i = seurat_cellbin
labels_hclust2 = merge(labels_hclust %>% tibble::rownames_to_column('big_cell'), cells_map, by.x='big_cell', by.y='merged_cluster')
seurat_i$'subclone' = 0
seurat_i@meta.data[labels_hclust2$cell_id, 'subclone'] = labels_hclust2$subclone
seurat_i = subset(seurat_i, subset=subclone!=0)
seurat_i@meta.data$subclone2 = paste0(seurat_i@meta.data$subclone, '_', seurat_i@meta.data$annotated_cluster)
`%nin%` <- Negate(`%in%`)
seurat_i = subset(seurat_i, subset=(subclone2 %nin% c('Subclo_A_DCIS cells', 'Subclo_A_IBC cells', 'Subclo_B_Atypical hyperplasia epithelial cells', 'Subclo_C_Atypical hyperplasia epithelial cells')))
seurat_i@meta.data$subclone2 = Recode(seurat_i@meta.data$subclone2, "'Subclo_A_Atypical hyperplasia epithelial cells'='Subclo_A';'Subclo_B_DCIS cells'='Subclo_B_DCIS';'Subclo_B_IBC cells'='Subclo_B_IBC';'Subclo_C_DCIS cells'='Subclo_C_DCIS';'Subclo_C_IBC cells'='Subclo_C_IBC';'Subclo_D_DCIS cells'='Subclo_D_DCIS';'Subclo_D_IBC cells'='Subclo_D_IBC';'Subclo_E_DCIS cells'='Subclo_E_DCIS';'Subclo_E_IBC cells'='Subclo_E_IBC'")
saveRDS(seurat_i, file='P13_subclone.rds')
system(sprintf('Rscript /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/codes/BRCA_NCC_20220729/MuDataSeurat_rds2h5ad.R %s/P13_subclone.rds %s/P13_subclone.h5ad', workdir, workdir))

#*P28
sample = 'P28'
st_id = st_ids[sample]
seurat_cellbin = ReadH5AD(paste0('/ST/', st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds/P28_raw_annotated.h5ad'))
cells_map = read.delim(paste0('/ST/', st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds/merged_cluster_labels.txt'), colClasses = "character")
labels_hclust = read.delim(paste0('/ST/', st_id, '/Results/8.30_cellbin_bin100_tree/8.30_labels_subclone.txt'), sep = '\t')

seurat_i = seurat_cellbin
labels_hclust2 = merge(labels_hclust %>% tibble::rownames_to_column('big_cell'), cells_map, by.x='big_cell', by.y='merged_cluster')
seurat_i$'subclone' = 0
seurat_i@meta.data[labels_hclust2$cell_id, 'subclone'] = labels_hclust2$subclone
seurat_i = subset(seurat_i, subset=subclone!=0)
seurat_i@meta.data$subclone2 = paste0(seurat_i@meta.data$subclone, '_', seurat_i@meta.data$annotated_cluster)
`%nin%` <- Negate(`%in%`)
seurat_i@meta.data$subclone2 = Recode(seurat_i@meta.data$subclone2, "'Subclo_A_DCIS cells'='Subclo_A_DCIS';'Subclo_A_IBC cells'='Subclo_A_IBC';'Subclo_B_DCIS cells'='Subclo_B_DCIS';'Subclo_B_IBC cells'='Subclo_B_IBC';'Subclo_C_DCIS cells'='Subclo_C_DCIS';'Subclo_C_IBC cells'='Subclo_C_IBC'")
saveRDS(seurat_i, file='P28_subclone.rds')
system(sprintf('Rscript MuDataSeurat_rds2h5ad.R %s/P28_subclone.rds %s/P28_subclone.h5ad', workdir, workdir))
