library('Seurat')
library('infercnv')
library('dendextend')
library('phylogram')
library('miscTools')
library('ggthemes')
library('RColorBrewer')
library('umap')
library('ggplot2')
library('car')
library('limma')
library('tibble')
library('dplyr')
library('stringr')
library('parallel')
library('yaGST')
options(bitmapType = "cairo")

args = commandArgs(T)
sample = args[1]
st_ids = c('P13'='148BR_D1_13-1', 'P15'='148BR_F2_15-1', 'P161'='148TR_E1_16-1', 'P22'='22-1', 'P28'='28-1', 'P29'='29-1', 'P35'='35-1', 'P36'='36-1', 'P41'='41-1')
st_id = st_ids[sample]
indir = paste0(st_id, '/Results/8.27_cellbin_bin100_infercnv_noCAF_v4')
indir2 = paste0(st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds')
indir3 = paste0(st_id, '/Results/8.28.2_cellbin_bin100_analyze_infercnv_selectFeatures')
workdir=paste0(st_id, '/Results/8.30_cellbin_bin100_tree')
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

##load data
infercnv_obj = readRDS(paste0(indir3, '/infercnv_obj2.rds'))
expr = infercnv_obj@expr.data
gene_info = read.delim(paste0(indir, '/gene_ordering_file.txt'), header = F, sep = '\t')
threds = read.delim(paste0(indir3, '/infercnv.heatmap_thresholds.txt'))[,1]
thred_high = threds[threds>1] %>% sort() %>% first()
thred_low = threds[threds<1] %>% sort() %>% last()
genesite <- gmt2GO("MSigDB/c1.all.v7.5.1.symbols.gmt")
pq_genes = data.frame()
for (pq in setdiff(names(genesite), 'MT')) {
    pq_i = data.frame(region=str_extract(pq,"chr[0-9A-Z]*[pq]"), gene=genesite[[pq]])
    pq_genes = rbind(pq_genes, pq_i)
}
dropped_genes = setdiff(rownames(expr), pq_genes$gene)
labels_hclust = read.delim(paste0(indir3, '/8.28.2_labels_subclone.txt'), sep = '\t')
labels_hclust = labels_hclust %>% filter(subclone != 'CNV_Normal')

##check positions of genes with CNV
expr2 = t(expr)[rownames(labels_hclust),]
expr2[expr2>thred_high] = 2 #gain
expr2[expr2<thred_low] = 0 #loss
expr2[(expr2>=thred_low)&(expr2<=thred_high)] = 1 #neutral
cnv_genes = list()
cnv_pqs = list()
cnv_df = data.frame()
thred_pcts = c('P13'=0.2, 'P28'=0.3, 'P35'=0.3, 'P41'=0.3, 'P15'=0.3, 'P161'=0.3, 'P22'=0.2, 'P29'=0.3, 'P36'=0.3)
thred_pct = thred_pcts[sample]
for (subclone in sort(unique(labels_hclust)[,1])) {
    cnv_genes_i = c()
    for (gene in intersect(colnames(expr2), pq_genes$gene)) {
        expr_i = expr2[rownames(labels_hclust)[labels_hclust[,1]==subclone],gene]
        gain_pct = sum(expr_i==2)/length(expr_i)
        loss_pct = sum(expr_i==0)/length(expr_i)
        if (gain_pct>thred_pct) {
            expr_i2 = expr[gene, rownames(labels_hclust)[labels_hclust[,1]==subclone]]
            cnv_df = rbind(cnv_df, c(sample, subclone, gene, 'gain', mean(expr_i2), median(expr_i2), sd(expr_i2), gain_pct))
            cnv_genes_i = c(cnv_genes_i, paste0(gene, '_gain'))
        } else if (loss_pct>thred_pct) {
            expr_i2 = expr[gene, rownames(labels_hclust)[labels_hclust[,1]==subclone]]
            cnv_df = rbind(cnv_df, c(sample, subclone, gene, 'loss', mean(expr_i2), median(expr_i2), sd(expr_i2), loss_pct))
            cnv_genes_i = c(cnv_genes_i, paste0(gene, '_loss'))
        }
    }
    cnv_genes[[subclone]] = cnv_genes_i
    pq_genes_i = intersect(substr(cnv_genes_i, 1, nchar(cnv_genes_i)-5), pq_genes$gene)
    cnv_pqs[[subclone]] = pq_genes %>% filter(gene %in% pq_genes_i) %>% tibble::column_to_rownames('gene') %>% .[pq_genes_i,] %>% paste0(., substring(cnv_genes_i, nchar(cnv_genes_i)-4)) %>% unique()
}
cnv_genes
for (subclo in names(cnv_pqs)) {
    cnv_pqs[[subclo]] = sort(cnv_pqs[[subclo]])
}
cnv_pqs
colnames(cnv_df) = c('sample', 'subclone', 'gene', 'if_gain', 'mean_cnv_score', 'median_cnv_score', 'sd_cnv_score', 'cnv_pct')
cnv_df2 = merge(cnv_df, pq_genes, by='gene')
cnv_df2$cnv_id = cnv_df2 %>% select(gene, if_gain, region) %>% apply(., 1, paste, collapse = "_")
cnv_df2$cnv_pq = cnv_df2 %>% select(region, if_gain) %>% apply(., 1, paste, collapse = "_")
cnv_df2 = cnv_df2 %>% group_by(subclone, cnv_pq) %>% mutate(mean_cnv_pct=mean(as.numeric(cnv_pct))) %>% as.data.frame()
saveRDS(cnv_genes, file = paste0(workdir, '/8.30_cnv_genes.rds'))
saveRDS(cnv_pqs, file = paste0(workdir, '/8.30_cnv_pqs.rds'))
write.table(cnv_df2, paste0(workdir, '/8.30_cnv_df2.txt'), row.names=F, col.names=T, quote=F, sep='\t')
cnv_pqs_df = do.call(cbind, lapply(lapply(cnv_pqs, unlist), 'length<-', max(lengths(cnv_pqs))))
cnv_pqs_df[is.na(cnv_pqs_df)] = ''
write.table(cnv_pqs_df, paste0(workdir, '/8.30_cnv_pqs_df.txt'), quote=F, row.names=F, col.names=T, sep='\t')

