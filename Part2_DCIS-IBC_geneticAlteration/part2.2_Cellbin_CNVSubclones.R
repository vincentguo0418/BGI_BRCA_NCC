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
library('MuDataSeurat')
options(bitmapType = "cairo")

args = commandArgs(T)
sample = args[1] #P13
metadata_path = args[2] #SS200000148BR_D1_metadata.txt, P13_cellbin_metadata_annotated.txt
mask_ft_path = args[3] #SS200000148BR_D1_regist_mask_ft.tif
k1 = args[4]
st_ids = c('P13'='148BR_D1_13-1', 'P15'='148BR_F2_15-1', 'P161'='148TR_E1_16-1', 'P22'='22-1', 'P28'='28-1', 'P29'='29-1', 'P35'='35-1', 'P36'='36-1', 'P41'='41-1')
st_id = st_ids[sample]
indir = paste0(st_id, '/Results/8.27_cellbin_bin100_infercnv_noCAF_v4')
indir3 = paste0(st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds')
workdir=paste0(st_id, '/Results/8.28.2_cellbin_bin100_analyze_infercnv_selectFeatures')
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

k1 = 18
##Identify CNV subclones
tree = read.dendrogram(paste0(indir, '/infercnv.observations_dendrogram.txt'))
labels_hclust = cutree(tree, k=k1)
labels_hclust = data.frame(row.names = names(labels_hclust), 'subclone'=paste0('Subclo_', labels_hclust))
subdends = get_subdendrograms(as.dendrogram(tree), k=k1, order_clusters_as_data = T)
load(paste0(indir, '/infercnv_obj2.RData'))
expr = infercnv_obj@expr.data
gene_info = read.delim(paste0(indir, '/gene_ordering_file.txt'), header = F, sep = '\t')

labels_hclust$subclone = Recode(labels_hclust$subclone, "c('Subclo_8', 'Subclo_9', 'Subclo_10')='CNV_Normal';c('Subclo_1', 'Subclo_2', 'Subclo_3', 'Subclo_4', 'Subclo_5', 'Subclo_6')='Subclo_B';c('Subclo_16', 'Subclo_17')='Subclo_C';c('Subclo_11', 'Subclo_12')='Subclo_D';c('Subclo_13', 'Subclo_14', 'Subclo_15', 'Subclo_18')='Subclo_E';'Subclo_7'='Subclo_A'")

selected_genes1 = intersect(gene_info %>% filter(V2=='chr16') %>% pull(V1), rownames(expr))
expr_selected = expr[selected_genes1,rownames(labels_hclust)[labels_hclust$subclone %in% paste0('Subclo_', 'A')]]
out.dist = wordspace::dist.matrix(t(expr_selected), method='euclidean')
expr_selected_hc = hclust(as.dist(out.dist), method='ward.D')
labels_subdend4 = cutree(expr_selected_hc, k=5)
table(labels_subdend4)
labels_hclust[names(labels_subdend4), 'subclone'] = paste0('Subclo_A.', labels_subdend4)

labels_hclust$subclone = Recode(labels_hclust$subclone, "c('Subclo_A.1', 'Subclo_A.3', 'Subclo_A.4')='Subclo_A';c('Subclo_A.2', 'Subclo_A.5')='CNV_Normal'")

table(labels_hclust)
labels_hclust = labels_hclust %>% arrange(subclone)
labels_hclust2 = labels_hclust
labels_hclust_withNormal = labels_hclust
labels_hclust2 = labels_hclust_withNormal %>% filter(subclone!='CNV_Normal')
write.table(labels_hclust2, paste0(workdir, '/8.28.2_labels_subclone.txt'), sep = '\t', row.names = T, col.names = T, quote = F)
infercnv_obj2 = infercnv_obj
ref_ids = setdiff(colnames(infercnv_obj2@expr.data), rownames(labels_hclust_withNormal))
infercnv_obj2@expr.data = infercnv_obj2@expr.data[,union(rownames(labels_hclust2), ref_ids)]
name = list()
index = list()
index_length = list()
hc_list = list()
for (i in unique(labels_hclust2[,1])) {
    name[[i]] <- rownames(labels_hclust2)[labels_hclust2$subclone==i]
    index[[i]] <- which(colnames(infercnv_obj2@expr.data) %in% name[[i]])
    index_length[[i]] <- seq(from=1,to=length(name[[i]]))
    hc_list[[i]] = list(order=seq(from=1,to=length(name[[i]])),labels=name[[i]])
}
infercnv_obj2@observation_grouped_cell_indices <- index_length
infercnv_obj2@tumor_subclusters$hc <- hc_list
infercnv_obj2@reference_grouped_cell_indices = list('ref'=which(colnames(infercnv_obj2@expr.data) %in% ref_ids))
source('SC_1.4.10_infercnv_heatmap_plot.R')
saveRDS(infercnv_obj2, file=paste0(workdir, '/infercnv_obj2.rds'))
##Plot CNV heatmap(Fig2c)
plot_cnv(infercnv_obj2,
        out_dir=workdir,
        title="inferCNV",
        obs_title="Observations (Cells)",
        ref_title="References (Cells)",
        cluster_by_groups=T,
        cluster_references=FALSE,
        plot_chr_scale=FALSE,
        chr_lengths=NULL,
        k_obs_groups = 1,
        contig_cex=1,
        x.center=mean(infercnv_obj2@expr.data),
        x.range="auto", #NA,
        hclust_method='ward.D',
        custom_color_pal=NULL,
        color_safe_pal=FALSE,
        output_filename="infercnv",
        output_format="pdf", #pdf, png, NA
        png_res=300,
        dynamic_resize=0,
        ref_contig = NULL,
        write_expr_matrix=F,
        useRaster=TRUE)

##Plot Cellbin spatial images colored by CNV subclones(Fig2d)
labels_cellID2bin = read.delim(sprintf('ST/%s/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds/merged_cluster_labels.txt', st_id), sep='\t')
labels_cellID2bin2 = read.delim(sprintf('ST/%s/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds/merged_cluster_labels.txt', st_id), sep='\t', colClasses = "character")
labels_cellID2bin[,1] = labels_cellID2bin2[,1]
labels_merge = merge(labels_cellID2bin, labels_hclust2 %>% tibble::rownames_to_column('merged_cluster'), by='merged_cluster')
rownames(labels_merge) = labels_merge$cell_id
cell_ids = rownames(labels_merge)
metadata = read.delim(metadata_path, sep = '\t')
metadata2 = read.delim(metadata_path, sep = '\t', colClasses = "character")
metadata[,1] = metadata2[,1]
rownames(metadata) = metadata[,1]
metadata$annotated_cluster = 'Low_quality'
metadata[cell_ids, 'annotated_cluster'] = labels_merge[cell_ids, 'subclone']
write.table('1\t#ff0000\nLow_quality\t#ffffff', paste0(workdir, '/color.list'), quote=F, row.names=F, col.names=F)
res = metadata[,c("x", "y", 'annotated_cluster')]
write.table(res, paste0(workdir, '/bin1clu.txt'), quote=F, row.names=F, col.names=F, sep='\t')
clusters = setdiff(unique(metadata$annotated_cluster), 'Low_quality')
plot_cellbin = function(i) {
    tmp = res
    tmp[tmp$annotated_cluster != i, 'annotated_cluster'] = 'Low_quality'
    tmp[tmp$annotated_cluster == i, 'annotated_cluster'] = '1'
    write.table(tmp, paste0(workdir, sprintf('/bin1clu_%s.txt', i)), quote=F, row.names=F, col.names=F, sep='\t')
    system(sprintf('cell_bin_plot %s/bin1clu_%s.txt %s %s/color.list %s/cluster_plot_%s.tif', workdir, i, mask_ft_path, workdir, workdir, i))
}
clus <- makeCluster(min(20, length(clusters)))
clusterExport(clus,c('res', 'workdir', 'mask_ft_path'), envir = environment())
clusterEvalQ(clus, library('Cairo'))
parLapply(clus, clusters, fun = plot_cellbin)
stopCluster(clus)
