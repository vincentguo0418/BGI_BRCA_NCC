library(dplyr)
library(Seurat)
library(MuDataSeurat)
library(limma)
library(car)

st_ids = c('P13'='148BR_D1_13-1', 'P15'='148BR_F2_15-1', 'P161'='148TR_E1_16-1', 'P22'='22-1', 'P28'='28-1', 'P29'='29-1', 'P35'='35-1', 'P36'='36-1', 'P41'='41-1')
samples = names(st_ids)
workdir=paste0('/paper/figure2F_P13_subcloneBCD_deg')
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

sample = 'P13'
st_id = st_ids[sample]
seurat_cellbin = ReadH5AD(paste0('/ST/', st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds/merged_adata.h5ad'))
labels_hclust = read.delim(paste0('/ST/', st_id, '/Results/8.30_cellbin_bin100_tree/8.30_labels_subclone.txt'), sep = '\t')

seurat_i = seurat_cellbin
seurat_i$subclone = 0
seurat_i@meta.data[rownames(labels_hclust), 'subclone'] = labels_hclust$subclone
seurat_i = subset(seurat_i, subset=subclone!=0)
source('DEG_wilcox.R')

#compare Subclo_B & Subclo_D(for Fig4s&Extended Data Fig8d)
seurat_bd = subset(seurat_i, subset=subclone %in% c('Subclo_B', 'Subclo_D'))
samplesCluster = list()
for (subclo in sort(unique(seurat_bd$subclone))) {
    samplesCluster[[subclo]] = rownames(seurat_bd@meta.data)[seurat_bd@meta.data$subclone==subclo]
}
options(future.globals.maxSize= 89128960000)
degs_raw_bd = DEG_wilcox2(samplesCluster, as.matrix(seurat_bd@assays$RNA@data))
saveRDS(degs_raw_bd, file=paste0(workdir, '/degs_raw_bd.rds'))
degs_raw_bd = readRDS(paste0(workdir, '/degs_raw_bd.rds'))
degs_raw_bd_df = do.call(rbind, degs_raw_bd) %>% as.data.frame()
degs_bd_df_b = degs_raw_bd_df %>% filter(Pval<0.05 & FDR<0.05 & MeanSelect>0.2 & FC>log2(1.5)) %>% mutate('upregulate'='b')
degs_bd_df_d = degs_raw_bd_df %>% filter(Pval<0.05 & FDR<0.05 & MeanOther>0.2 & FC<(-log2(1.5))) %>% mutate('upregulate'='d')
degs_bd_df = rbind(degs_bd_df_b, degs_bd_df_d)
genesite <- yaGST::gmt2GO("MSigDB/c1.all.v7.5.1.symbols.gmt")
pq_genes = data.frame()
for (pq in setdiff(names(genesite), c('MT', "chrYp11", "chrYq11"))) {
    pq_i = data.frame(region=stringr::str_extract(pq,"chr[0-9A-Z]*[pq]"), gene=genesite[[pq]])
    pq_genes = rbind(pq_genes, pq_i)
}
degs_bd_df = merge(degs_bd_df %>% tibble::rownames_to_column('gene'), pq_genes, by='gene')
source('signatureGenes.R')
nonepi_markers = SigGeneral_all[!(names(SigGeneral_all) %in% c('Epithelial', 'Epithe_Lum', 'Epithe_Mam'))] %>% unlist() %>% unique()
epi_markers = SigGeneral_all[(names(SigGeneral_all) %in% c('Epithelial', 'Epithe_Lum', 'Epithe_Mam'))] %>% unlist() %>% unique()
degs_bd_df = degs_bd_df %>% filter(!(gene %in% nonepi_markers))
rmgenes=removeBiasGenes(degs_bd_df$gene)
`%nin%` <- Negate(`%in%`)
degs_bd_df2 = degs_bd_df %>% filter(gene %nin% rmgenes)
degs_bd_df2 = degs_bd_df2[!grepl('-', degs_bd_df2$gene),]
degs_bd_df2_b = degs_bd_df2 %>% filter(FC>log2(1.5) & MeanSelect>0.2)
degs_bd_df2_d = degs_bd_df2 %>% filter(FC<(-2) & MeanOther>0.8)
write.table(degs_bd_df2_b, paste0(workdir, '/degs_bd_df2_b.txt'), sep='\t', quote=F, row.names=F, col.names=T)
write.table(degs_bd_df2_d, paste0(workdir, '/degs_bd_df2_d.txt'), sep='\t', quote=F, row.names=F, col.names=T)
