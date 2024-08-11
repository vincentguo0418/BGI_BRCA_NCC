library('Seurat')
library('infercnv')
library('dplyr')
library('car')
library('AnnoProbe')
library('tidyverse')
library('limma')
library('MuDataSeurat')
options(bitmapType = "cairo")


args = commandArgs(T)
sample = args[1]
st_ids = c('P13'='148BR_D1_13-1', 'P15'='148BR_F2_15-1', 'P161'='148TR_E1_16-1', 'P22'='22-1', 'P28'='28-1', 'P29'='29-1', 'P35'='35-1', 'P36'='36-1', 'P41'='41-1')
st_id = st_ids[sample]
workdir=paste0(st_id, '/Results/8.27_cellbin_bin100_infercnv_noCAF_v4')
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

seurat_cellbin_bin100 = ReadH5AD(paste0(st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds/merged_adata.h5ad'))
ref_group_names_list = list('P13'=c('MEC', 'Normal epithelial cells'),
                            'P161'=c('MEC', 'Normal epithelial cells'),
                            'P15'=c('MEC', 'Normal epithelial cells'),
                            'P22'=c('MEC', 'Normal epithelial cells'),
                            'P28'=c('MEC', 'Normal epithelial cells'),
                            'P29'=c('MEC', 'Normal epithelial cells'),
                            'P35'=c('MEC', 'Normal epithelial cells'),
                            'P36'=c('MEC', 'Normal epithelial cells'),
                            'P41'=c('MEC', 'Normal epithelial cells'))
ref_group_names = ref_group_names_list[[sample]]
obs_group_names_list = list('P13'=c('DCIS cells', 'IBC cells', 'Atypical hyperplasia epithelial cells'),
                            'P161'=c('DCIS cells', 'IBC cells'),
                            'P15'=c('DCIS cells', 'IBC cells'),
                            'P22'=c('DCIS cells', 'IBC cells'),
                            'P28'=c('DCIS cells', 'IBC cells'),
                            'P29'=c('DCIS cells', 'IBC cells'),
                            'P35'=c('DCIS cells', 'IBC cells'),
                            'P36'=c('DCIS cells', 'IBC cells'),
                            'P41'=c('DCIS cells', 'IBC cells', 'Atypical hyperplasia epithelial cells'))
obs_group_names = obs_group_names_list[[sample]]
seurat_i = subset(seurat_cellbin_bin100, subset=annotated_cluster %in% c(ref_group_names, obs_group_names))
selected_clusters = table(seurat_i$annotated_cluster) %>% as.data.frame() %>% filter(Freq>20) %>% pull(Var1) %>% as.character()
ref_group_names = intersect(ref_group_names, selected_clusters)
seurat_i = subset(seurat_cellbin_bin100, subset=annotated_cluster %in% c(ref_group_names, obs_group_names))
#obs_group_names = intersect(obs_group_names, selected_clusters)

mat_comb <- as.matrix(GetAssayData(seurat_i, slot = "counts"))

geneinfo <- annoGene(rownames(mat_comb), "SYMBOL", "human")
geneinfo <- geneinfo[with(geneinfo,order(chr,start)),c(1,4:6)]
geneinfo <- geneinfo[!duplicated(geneinfo[,1]),]

rownames(geneinfo) <- geneinfo$SYMBOL

geneinfo$SYMBOL <- NULL
geneinfo = geneinfo %>% mutate(chr_n=strsplit2(geneinfo$chr, split='chr')[,2]) %>%
 mutate(chr_n = replace(chr_n, chr_n=='X', 23)) %>%
 mutate(chr_n = replace(chr_n, chr_n=='Y', 24)) %>%
 mutate(chr_n = replace(chr_n, chr_n=='M', 25)) %>%
 mutate(chr_n = as.numeric(chr_n)) %>%
 arrange(chr_n, start, end) %>%
 select(-chr_n)

seurat_i$cluster <- seurat_i$annotated_cluster
metada <- seurat_i@meta.data
metadata <- metada[,c("cluster", "cluster")]
colnames(metadata) <- c("cellType", "group")

write.table(metadata, file=paste0(workdir,"/cellAnnotations.txt"), sep="\t", col.names = FALSE,quote = FALSE)
write.table(geneinfo, file=paste0(workdir,"/gene_ordering_file.txt"), sep = "\t", col.names = FALSE, quote = FALSE)

count_mat <- mat_comb[rownames(geneinfo),]

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = count_mat, 
                                     gene_order_file = paste0(workdir,"/gene_ordering_file.txt"), 
                                     delim = "\t", min_max_counts_per_cell = c(10, +Inf),
                                     ref_group_names = ref_group_names,
                                     annotations_file = paste0(workdir,"/cellAnnotations.txt"))
save(infercnv_obj, file = paste0(workdir,"/infercnv_obj1.RData"))
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir='.',  # dir is auto-created for storing outputs
                             cluster_by_groups=F,   # cluster
                             scale_data=FALSE,
                             denoise=T,
                             sd_amplifier=1.5,
                             HMM=F,
                             #analysis_mode='subclusters',
			     num_threads=50
                             )
save(infercnv_obj, file = paste0(workdir,"/infercnv_obj2.RData"))
