library(dplyr)
library(tibble)
library(limma)
library(Seurat)
library(MuDataSeurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggsci)
library(dplyr)
library(tibble)
library(limma)
library(Seurat)
options(stringsAsFactors = FALSE)

args = commandArgs(T)
sample = args[1]
st_ids = c('P13'='148BR_D1_13-1', 'P15'='148BR_F2_15-1', 'P161'='148TR_E1_16-1', 'P22'='22-1', 'P28'='28-1', 'P29'='29-1', 'P35'='35-1', 'P36'='36-1', 'P41'='41-1')
samples = names(st_ids)
st_id = st_ids[sample]
workdir=paste0(topdir, st_id, '/Results/8.41.2_cell_communication')
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

##Cell types interacting with MEC(Fig6f)
sample_file = c('P13'='raw_annotated_bdl_4.h5ad', 'P28'='raw_annotated_bdl_0.h5ad', 'P35'='raw_annotated_bdl_0123.h5ad')
raw_adata_i = ReadH5AD(sprintf('ST/%s/Results/5_cellbin_230105/8.41_borderline/%s', st_id, sample_file[[sample]]))
raw_adata_i_ext = subset(raw_adata_i, subset=borderline_surrounding_500=='external')
raw_adata_i_int = subset(raw_adata_i, subset=borderline_surrounding_500=='internal')
WriteH5AD(raw_adata_i_ext, file='raw_annotated_bdl_i_external.h5ad')
WriteH5AD(raw_adata_i_int, file='raw_annotated_bdl_i_internal.h5ad')

system(sprintf('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/bin/python /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/codes/BRCA_NCC_20220729/8.23.4_cellbin_merge_cluster_bin100_merge2rds.py -a /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/%s/Results/5_cellbin_230105/8.41.2_cell_communication/raw_annotated_bdl_i_external.h5ad -o %s/external_i -n 10', st_id, workdir))
system(sprintf('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/bin/python /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/codes/BRCA_NCC_20220729/8.23.4_cellbin_merge_cluster_bin100_merge2rds.py -a /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/%s/Results/5_cellbin_230105/8.41.2_cell_communication/raw_annotated_bdl_i_internal.h5ad -o %s/internal_i -n 10', st_id, workdir))

##单克隆
raw_adata_others = ReadH5AD(sprintf('/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/%s/Results/5_cellbin_230105/8.41_borderline/raw_annotated_bdl_others.h5ad', st_id))
raw_adata_others_ext = subset(raw_adata_others, subset=borderline_surrounding_500=='external')
raw_adata_others_int = subset(raw_adata_others, subset=borderline_surrounding_500=='internal')
WriteH5AD(raw_adata_others_ext, file='raw_annotated_bdl_others_external.h5ad')
WriteH5AD(raw_adata_others_int, file='raw_annotated_bdl_others_internal.h5ad')

system(sprintf('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/bin/python /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/codes/BRCA_NCC_20220729/8.23.4_cellbin_merge_cluster_bin100_merge2rds.py -a /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/%s/Results/5_cellbin_230105/8.41.2_cell_communication/raw_annotated_bdl_others_external.h5ad -o %s/external_others -n 10', st_id, workdir))
system(sprintf('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/software/miniconda/bin/python /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/software/codes/BRCA_NCC_20220729/8.23.4_cellbin_merge_cluster_bin100_merge2rds.py -a /jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/project/BRCA_NCC_20220729/ST/%s/Results/5_cellbin_230105/8.41.2_cell_communication/raw_annotated_bdl_others_internal.h5ad -o %s/internal_others -n 10', st_id, workdir))

adata_ext_others = ReadH5AD('external_others/merged_adata.h5ad')
adata_int_others = ReadH5AD('internal_others/merged_adata.h5ad')
table(adata_ext_others$annotated_cluster)
table(adata_int_others$annotated_cluster)
adata_ext_others$'borderline_surrounding_500' = 'ext'
adata_int_others$'borderline_surrounding_500' = 'int'
merge_adata_others = merge(adata_ext_others, adata_int_others)
table(merge_adata_others$annotated_cluster)
table(merge_adata_others$borderline_surrounding_500)
merge_adata_others$source = paste(merge_adata_others$borderline_surrounding_500, merge_adata_others$annotated_cluster, sep='_')
merge_adata_others$source = car::Recode(merge_adata_others$source, "c('int_MEC', 'ext_MEC')='MEC'")

##合并
adata_ext_i = ReadH5AD('external_i/merged_adata.h5ad')
adata_int_i = ReadH5AD('internal_i/merged_adata.h5ad')
adata_ext_others = ReadH5AD('external_others/merged_adata.h5ad')
adata_int_others = ReadH5AD('internal_others/merged_adata.h5ad')
# table(adata_ext_i$annotated_cluster)
# table(adata_int_i$annotated_cluster)
adata_ext_i$'borderline_surrounding_500' = 'ext'
adata_int_i$'borderline_surrounding_500' = 'int'
adata_ext_others$'borderline_surrounding_500' = 'ext'
adata_int_others$'borderline_surrounding_500' = 'int'
merge_adata_i = merge(adata_ext_i, adata_int_i)
merge_adata_others = merge(adata_ext_others, adata_int_others)
# table(merge_adata_i$annotated_cluster)
# table(merge_adata_i$borderline_surrounding_500)
merge_adata_all = merge(merge_adata_i, merge_adata_others)
merge_adata_all$source = paste(merge_adata_all$borderline_surrounding_500, merge_adata_all$annotated_cluster, sep='_')
merge_adata_all$source = car::Recode(merge_adata_all$source, "c('int_MEC', 'ext_MEC')='MEC'")

##统计
source_stat = table(merge_adata_all$source) %>% as.data.frame()
exclude_clusters = source_stat %>% filter(Freq<=5) %>% pull(Var1) %>% as.character()
exclude_clusters = union(exclude_clusters, c('ext_DCIS cells', 'ext_low_quality', 'int_low_quality'))
'%nin%' <- Negate('%in%')
merge_adata_all2 = subset(merge_adata_all, subset=(source %nin% exclude_clusters))
merge_adata_all2 = NormalizeData(merge_adata_all2)
cell_nums = table(merge_adata_all2$annotated_cluster, merge_adata_all2$borderline_surrounding_500) %>% as.data.frame.matrix()
write.table(cell_nums, 'cell_nums.txt', quote=F, col.names=T, row.names=T, sep='\t')

##CellChat: merge_adata_i
data.input <- merge_adata_all2@assays$RNA@data # normalized data matrix
labels <- merge_adata_all2$source
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat <- createCellChat(object = merge_adata_all2, group.by = "source", assay = "RNA")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# db.new = readRDS("/jdfsbjcas1/ST_BJ/P21H28400N0232/guojingwen/database/CellphoneDB/cellphonedb-data-5.0.0/CellChatDB.human_user.rds")
# cellchat@DB <- db.new
CellChatDB <- readRDS('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/database/cellchat/human/CellChatDB.human_nichenet.rds')
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 5)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
saveRDS(df.net, 'df.net.rds')
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, 'cellchat.rds')

groupSize <- as.numeric(table(cellchat@idents))
pdf(file = paste0(workdir, '/8.41.2_netVisual_circle.pdf'), width = 10, height = 6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
pdf(file = paste0(workdir, '/8.41.2_netVisual_circle2.pdf'), width = 15, height = 20)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
mat <- cellchat@net$weight
i = 'MEC'
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[,i] <- mat[,i]
groupSize <- as.numeric(table(cellchat@idents))
pdf(file = paste0(workdir, '/8.41.2_netVisual_circle_toMEC.pdf'), width = 7, height = 7)
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
dev.off()

##Cell type distribution within the internal and external regions(Fig6e)
cell_nums_all = data.frame()
for (sample in samples) {
    st_id = st_ids[sample]
    indir = paste0(topdir, st_id, '/Results/8.41.2_cell_communication')
    cell_nums = read.delim(paste0(indir, '/cell_nums.txt'))
    cell_nums2 = reshape2::melt(cell_nums %>% rownames_to_column('Cell Type'), variable.name='Region', value.name='Cell Proportion i')
    cell_nums2$id = paste(cell_nums2$Region, cell_nums2$'Cell Type', sep='_')
    if (dim(cell_nums_all)[1]==0) {
        cell_nums_all = cell_nums2
    } else {
        cell_nums_all = merge(cell_nums_all, cell_nums2[,c('id', 'Cell Proportion i')], by='id', all.x=T)
    }
}
cell_nums_all[is.na(cell_nums_all)] = 0
cell_nums_all$'Cell Proportion all' = apply(cell_nums_all[,4:12], 1, sum)
cell_nums_all = cell_nums_all[,c('Cell Type', 'Region', 'Cell Proportion all')]
cell_nums_all$Region = car::Recode(cell_nums_all$Region, "'int'='Internal';'ext'='External'")
cell_nums_all$'Cell Type'[cell_nums_all$'Cell Type'=='Atypical hyperplasia epithelial cells'] = 'HEC'
cell_nums_all$'Cell Type' = factor(cell_nums_all$'Cell Type', levels=c('HEC', 'B cells', 'DCIS cells', 'Endothelial cells', 'Fibroblasts', 'IBC cells', 'MEC', 'Macrophages', 'Mast cells', 'Plasmocytes', 'Normal epithelial cells', 'Pericytes', 'T cells'), ordered=T)
cell_nums_all = cell_nums_all %>% rename('Cell Proportion'='Cell Proportion all')
cell_nums_all = tidyr::spread(cell_nums_all, key='Region', value='Cell Proportion')
cell_nums_all = cell_nums_all %>% column_to_rownames('Cell Type') %>% apply(., 2, function(x){x/sum(x)}) %>% as.data.frame()
cell_nums_all = reshape2::melt(cell_nums_all %>% rownames_to_column('Cell Type'), variable.name='Region', value.name='Cell Proportion')
pdf(file = paste0(workdir, '/8.41.3_cell_nums_all.pdf'), width = 8, height = 10)
ggplot(cell_nums_all, aes(x=Region, y=`Cell Proportion`)) +
    geom_bar(aes(fill=`Cell Type`), position='stack', stat = 'identity') +
    scale_fill_manual(values=c('HEC'='#8B1B6C', 'B cells'='#84C4B6', 'DCIS cells'='#E6400D', 'Endothelial cells'='#4E3A8A', 'Fibroblasts'='#C6C0A7', 'IBC cells'='#EEE919', 'MEC'='#86BC28', 'Macrophages'='#866236', 'Mast cells'='#415F8E', 'Plasmocytes'='#12B8D7', 'Normal epithelial cells'='#368039', 'Pericytes'='#DDC02C', 'T cells'='#DE3E90')) +
    theme(text = element_text(size = 20),
        panel.background=element_rect(fill="white"),
        plot.title = element_text(hjust=0.5, size = 22),
        legend.position = "right",
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
        labs(x=NULL) +
        guides(fill = guide_legend(title = NULL))
dev.off()
