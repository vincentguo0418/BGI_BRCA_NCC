library('yaGST')
library('stringr')
library('dplyr')
library('tidyverse')

##Calculate mean CNV score of DCIS and IBC cells across 9 Cellbin samples(Fig2f)
st_ids = c('P13'='148BR_D1_13-1', 'P15'='148BR_F2_15-1', 'P161'='148TR_E1_16-1', 'P22'='22-1', 'P28'='28-1', 'P29'='29-1', 'P35'='35-1', 'P36'='36-1', 'P41'='41-1')
workdir = 'Figure2_cnv_heatmap'
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)
genesite <- gmt2GO("MSigDB/c1.all.v7.5.1.symbols.gmt")
pq_genes = data.frame()
for (pq in setdiff(names(genesite), 'MT')) {
    pq_i = data.frame(region=str_extract(pq,"chr[0-9A-Z]*[pq]"), gene=genesite[[pq]])
    pq_genes = rbind(pq_genes, pq_i)
}
pq_mean_mtx_all = data.frame()
for (sample in names(st_ids)) {
      print(sample)
      st_id = st_ids[[sample]]
      indir = paste0(st_id, '/Results/8.27_cellbin_bin100_infercnv_noCAF_v4')
      load(paste0(indir, '/infercnv_obj2.RData'))
      expr = infercnv_obj@expr.data
      dropped_genes = setdiff(rownames(expr), pq_genes$gene)
      expr2 = as.data.frame(expr) %>% tibble::rownames_to_column('gene')
      expr2 = merge(expr2, pq_genes, by='gene')
      pq_mean_mtx = expr2 %>% select(-c(gene, region)) %>% split(expr2$region) %>% sapply(., function(x){apply(x, 2, mean)}) %>% as.data.frame()
      chr = str_extract(colnames(pq_mean_mtx), "[0-9A-Z]+")
      df_pqs = data.frame('region'=colnames(pq_mean_mtx), 'chr'=as.numeric(chr))
      df_pqs = df_pqs %>% arrange(chr)
      pq_mean_mtx = pq_mean_mtx[,df_pqs$region]

      indir = paste0('/', st_id, '/Results/8.30_cellbin_bin100_tree')
      indir2 = paste0('/', st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds')
      ids_table <- read.csv(paste0(indir, '/8.30_labels_subclone.txt'),header = TRUE,sep='\t')
      ids_table$merged_cluster <- rownames(ids_table)
      anno_table <- read.csv(paste0(indir2, '/merged_cluster_labels.txt'),header = TRUE,sep='\t')
      anno_table1 <- left_join(anno_table, ids_table, by='merged_cluster') %>% drop_na('subclone')
      if (sample=='P161') {
            sample2 = 'P16'
      } else {
            sample2 = sample
      }
      meta <- read.csv(paste0(indir2, '/', sample2, '_cellbin_metadata_annotated.txt'),header = TRUE,sep='\t')
      meta$cell_id <- meta$X
      anno_table1 <- left_join(anno_table1, meta, by='cell_id')
      anno_table2 <- filter(anno_table1, anno_table1$'annotated_cluster' %in% c('DCIS cells','IBC cells'))  #Atypical hyperplasia epithelial cells, 
      anno_table2$'Cell_type' <- anno_table2$'annotated_cluster'
      anno_table2 <- select(anno_table2, c('merged_cluster','Cell_type')) %>% distinct()
      pq_mean_mtx_anno = pq_mean_mtx %>% rownames_to_column('merged_cluster') %>% merge(., anno_table2, by='merged_cluster')
      pq_mean_mtx_anno$'sample' = sample
      pq_mean_mtx_all = rbind(pq_mean_mtx_all, pq_mean_mtx_anno)
}

pq_mean_mtx_all = pq_mean_mtx_all %>% arrange(sample, Cell_type)
pq_mean_mtx_all2 = pq_mean_mtx_all %>% select(-c('merged_cluster')) %>% column_to_rownames('cell')
pq_mean_mtx_all2 = reshape2::melt(pq_mean_mtx_all2, id.vars=c('Cell_type', 'sample')) %>% dplyr::group_by(sample, Cell_type, variable) %>% mutate(mean_pq=mean(value)) %>% select(-c('value')) %>% ungroup() %>% as.data.frame() %>% distinct() %>% tidyr::spread(., key='variable', value='mean_pq') %>% arrange(sample, Cell_type)
pq_mean_mtx_all2$'id' = paste0('id_', 1:nrow(pq_mean_mtx_all2))
pq_mean_mtx_all3 = pq_mean_mtx_all2 %>% column_to_rownames('id') %>% select(-c(sample, Cell_type))
pq_mean_mtx_all3 = as.data.frame(t(pq_mean_mtx_all3))
ann_col = data.frame(pq_mean_mtx_all2[,c('id', 'Cell_type', 'sample')]) %>% column_to_rownames('id') %>% rename('Cell type'='Cell_type', 'Sample'='sample')
thred_cut = max(1-min(pq_mean_mtx_all3), max(pq_mean_mtx_all3)-1)
pdf(file = paste0(workdir, '/Figure2_cnv_heatmap3.pdf'),width = 7, height = 6, onefile=F)
pq_mean_mtx_all3_ht = pheatmap::pheatmap(as.matrix(pq_mean_mtx_all3)-1, color = colorRampPalette(colors = c("blue","white","red"))(100), cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, scale='none', breaks = unique(c(seq(-thred_cut/2, thred_cut/2, length=100))), clustering_method='ward.D2', annotation_col=ann_col, border_color=NA)
print(pq_mean_mtx_all3_ht)
dev.off()
