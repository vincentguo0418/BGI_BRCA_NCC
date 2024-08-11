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
options(bitmapType = "cairo")

sample = args[1] #'P13'
st_ids = c('P13'='148BR_D1_13-1', 'P15'='148BR_F2_15-1', 'P161'='148TR_E1_16-1', 'P22'='22-1', 'P28'='28-1', 'P29'='29-1', 'P35'='35-1', 'P36'='36-1', 'P41'='41-1')
st_id = st_ids[sample]
indir = paste0('/Results/1.4.9_rerun_infercnv_onlyTumor/ref_NormalEpi/', sample)
workdir = paste0('/Results/1.4.13_SC_analyze_inferCNV_selectFeatures/', sample)
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

##Identify CNV subclones
tree = read.dendrogram(paste0(indir, '/infercnv.observations_dendrogram.txt'))
k1 = 10
labels_hclust = cutree(tree, k=k1)
labels_hclust = data.frame(row.names = names(labels_hclust), 'subclone'=paste0('Subclo_', labels_hclust))
load(paste0(indir, '/infercnv_obj2.RData'))
expr = infercnv_obj@expr.data
gene_info = read.delim(paste0(indir, '/gene_ordering_file.txt'), header = F, sep = '\t')

selected_genes1 = intersect(gene_info %>% filter(V2 %in% c('chr16')) %>% pull(V1), rownames(expr))
expr_selected = expr[selected_genes1,rownames(labels_hclust)[labels_hclust$subclone %in% paste0('Subclo_', 2)]]
selected_genes1 = rownames(expr_selected)[apply(expr_selected, 1, function(x){sum(x>1)>sum(x<=1)})]
expr_selected = expr[selected_genes1,rownames(labels_hclust)[labels_hclust$subclone %in% paste0('Subclo_', 2)]]
out.dist = wordspace::dist.matrix(t(expr_selected), method='euclidean')
expr_selected_hc = hclust(as.dist(out.dist), method='ward.D')
labels_subdend4 = cutree(expr_selected_hc, k=2)
table(labels_subdend4)
labels_hclust[names(labels_subdend4), 'subclone'] = paste0('Subclo_2.', labels_subdend4)

selected_genes1 = intersect(gene_info %>% filter(V2 %in% c('chr8')) %>% pull(V1), rownames(expr))
expr_selected = expr[selected_genes1,rownames(labels_hclust)[labels_hclust$subclone %in% paste0('Subclo_', 2.1)]]
selected_genes1 = rownames(expr_selected)[apply(expr_selected, 1, function(x){sum(x>1)>sum(x<=1)})]
expr_selected = expr[selected_genes1,rownames(labels_hclust)[labels_hclust$subclone %in% paste0('Subclo_', 2.1)]]
out.dist = wordspace::dist.matrix(t(expr_selected), method='euclidean')
expr_selected_hc = hclust(as.dist(out.dist), method='ward.D')
labels_subdend4 = cutree(expr_selected_hc, k=2)
table(labels_subdend4)
labels_hclust[names(labels_subdend4), 'subclone'] = paste0('Subclo_2.1.', labels_subdend4)

selected_genes1 = intersect(gene_info %>% filter(V2 %in% c('chr16')) %>% pull(V1), rownames(expr))
expr_selected = expr[selected_genes1,rownames(labels_hclust)[labels_hclust$subclone %in% paste0('Subclo_', 3:5)]]
selected_genes1 = rownames(expr_selected)[apply(expr_selected, 1, function(x){sum(x>1)>sum(x<=1)})]
expr_selected = expr[selected_genes1,rownames(labels_hclust)[labels_hclust$subclone %in% paste0('Subclo_', 3:5)]]
out.dist = wordspace::dist.matrix(t(expr_selected), method='euclidean')
expr_selected_hc = hclust(as.dist(out.dist), method='ward.D')
labels_subdend4 = cutree(expr_selected_hc, k=2)
table(labels_subdend4)
labels_hclust[names(labels_subdend4), 'subclone'] = paste0('Subclo_345.', labels_subdend4)

selected_genes1 = intersect(gene_info %>% filter(V2 %in% c('chr3')) %>% pull(V1), rownames(expr))
expr_selected = expr[selected_genes1,rownames(labels_hclust)[labels_hclust$subclone %in% paste0('Subclo_', 345.1)]]
selected_genes1 = rownames(expr_selected)[apply(expr_selected, 1, function(x){sum(x>1)>sum(x<=1)})]
expr_selected = expr[selected_genes1,rownames(labels_hclust)[labels_hclust$subclone %in% paste0('Subclo_', 345.1)]]
out.dist = wordspace::dist.matrix(t(expr_selected), method='euclidean')
expr_selected_hc = hclust(as.dist(out.dist), method='ward.D')
labels_subdend4 = cutree(expr_selected_hc, k=2)
table(labels_subdend4)
labels_hclust[names(labels_subdend4), 'subclone'] = paste0('Subclo_345.1.', labels_subdend4)

labels_hclust$subclone = Recode(labels_hclust$subclone, "c('Subclo_6', 'Subclo_7', 'Subclo_8', 'Subclo_9', 'Subclo_10')='Subclo_678910'")
labels_hclust$subclone = Recode(labels_hclust$subclone, "c('Subclo_2.1.1', 'Subclo_2.1.2')='Subclo_B';'Subclo_2.2'='Subclo_G';'Subclo_345.1.1'='Subclo_D';'Subclo_345.1.2'='Subclo_F';'Subclo_345.2'='Subclo_E';'Subclo_678910'='Subclo_C';'Subclo_1'='Subclo_H'")
labels_hclust$subclone = Recode(labels_hclust$subclone, "'Subclo_B'='Subclo_A';'Subclo_C'='Subclo_B';'Subclo_D'='Subclo_C';'Subclo_E'='Subclo_D';'Subclo_F'='Subclo_E';'Subclo_G'='Subclo_F';'Subclo_H'='Subclo_G'")

table(labels_hclust)
labels_hclust2 = labels_hclust
labels_hclust2 = labels_hclust2 %>% arrange(subclone)
write.table(labels_hclust2, paste0(workdir, '/1.4.13_labels_hclust2.txt'), row.names = T, col.names = T, sep = '\t', quote = F)
labels_hclust2 = read.delim(paste0(workdir, '/1.4.13_labels_hclust2.txt'), sep='\t')
labels_hclust_withNormal = labels_hclust2
labels_hclust2 = labels_hclust2 %>% filter(!(subclone %in% c('Subclo_F', 'Subclo_G')))

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
##Plot CNV heatmap(Fig2a)
plot_cnv(infercnv_obj2,
        out_dir=workdir,
        title="inferCNV_test",
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

##Plot spatial mapping of SC CNV subclones(Fig2b)
load(paste0('ST/', st_id, '/saveData/8.1.3_seurat_spatialObj_rmMT.RData')) #seurat_spatialObj
load(paste0("/saveData/1.4.3_seurat_merge_infercnv.RData"))
ref = c('Lum_mat', 'Lum_prog', 'myoepithelium')
seurat_i = subset(seurat_merge_infercnv, subset=((subtypes_latest=='Tumor'&Patients==paste0(sample, 'T'))|(subtypes_latest %in% ref)))
seurat_tumor = subset(seurat_i, subset = subtypes_latest=='Tumor')
seurat_normalEpi = subset(seurat_i, subset = subtypes_latest!='Tumor')
seurat_tumor@meta.data$subtypes_latest = labels_hclust2[rownames(seurat_tumor@meta.data),]
seurat_Epithe = merge(seurat_tumor, seurat_normalEpi)
seurat_Epithe$subtypes_latest = Recode(seurat_Epithe$subtypes_latest, "c('Lum_prog', 'Lum_mat', 'myoepithelium')='Epithelial_normal'")
seurat_Epithe <- NormalizeData(seurat_Epithe)
seurat_Epithe <- FindVariableFeatures(seurat_Epithe, selection.method = "vst", nfeatures = 2000)
seurat_Epithe <- ScaleData(seurat_Epithe)
seurat_Epithe <- RunPCA(seurat_Epithe)
seurat_Epithe@meta.data[rownames(seurat_tumor@meta.data), 'subtypes_latest'] = labels_hclust2[rownames(seurat_tumor@meta.data),]
seurat_tumor2 = seurat_tumor
seurat_tumor2 <- NormalizeData(seurat_tumor2)
seurat_tumor2 <- FindVariableFeatures(seurat_tumor2, selection.method = "vst", nfeatures = 2000)
seurat_tumor2 <- ScaleData(seurat_tumor2)
seurat_tumor2 <- RunPCA(seurat_tumor2)
seurat_tumor2@meta.data[rownames(seurat_tumor@meta.data), 'subtypes_latest'] = labels_hclust2[rownames(seurat_tumor@meta.data),]

anchors <- FindTransferAnchors(reference = seurat_Epithe, query = seurat_spatialObj, normalization.method = "LogNormalize")
save(anchors, file = paste0(workdir2, "/1.4.10_anchors.RData"))
singlecellsubTypes.assay.infercnv <- TransferData(anchorset = anchors, refdata = seurat_Epithe$subtypes_latest, prediction.assay = FALSE, weight.reduction = seurat_spatialObj[["pca"]], dims = 1:30)
seurat_spatialObj2_infercnv <- AddMetaData(object = seurat_spatialObj, metadata = singlecellsubTypes.assay.infercnv)
save(singlecellsubTypes.assay.infercnv, file = paste0(workdir2, "/1.4.10_singlecellsubTypes.assay.infercnv.RData"))
Idents(seurat_spatialObj2_infercnv) <- seurat_spatialObj2_infercnv$predicted.id
seurat_spatialObj2_infercnv$subTypes_predicted.id <- seurat_spatialObj2_infercnv$predicted.id
save(seurat_spatialObj2_infercnv, file = paste0(workdir2,"/1.4.10_seurat_spatialObj2_infercnv.RData"))
options(bitmapType = "cairo")
for(i in colnames(singlecellsubTypes.assay.infercnv)[-c(1,length(colnames(singlecellsubTypes.assay.infercnv)))]){
    print(i)
    tiff(file = paste0(workdir, "/1.4.13_seurat2space_infercnv", i, ".tiff"), res=300, width=9, height=8, compression="lzw", units="in")
    tryCatch(print(SpatialFeaturePlot(seurat_spatialObj2_infercnv, features = i, pt.size.factor = 1.3, crop = TRUE, alpha = c(0.1, 1), stroke = 0)), error = function(e){
        print(c(i,e))
    })
    dev.off()
}

saveRDS(seurat_Epithe, file=paste0(workdir, '/1.4.13_seurat_Epithe.rds'))
saveRDS(seurat_spatialObj, file=paste0(workdir, '/1.4.13_seurat_spatialObj.rds'))
system(sprintf('Rscript 4.rdstoh5ad2.R --infile %s/1.4.13_seurat_Epithe.rds --outfile %s/1.4.13_seurat_Epithe.h5ad', workdir, workdir))
system(sprintf('Rscript 4.rdstoh5ad2.R --infile %s/1.4.13_seurat_spatialObj.rds --outfile %s/1.4.13_seurat_spatialObj.h5ad', workdir, workdir))
file.remove(paste0(workdir, '/1.4.13_seurat_Epithe.h5Seurat'))
file.remove(paste0(workdir, '/1.4.13_seurat_spatialObj.h5Seurat'))
system(sprintf('python cellbin_SPACEL_Spoint.py %s/1.4.13_seurat_Epithe.h5ad %s/1.4.13_seurat_spatialObj.h5ad %s subtypes_latest', workdir, workdir, workdir))
spacel = read.delim(paste0(workdir, '/Spoint_metadata.txt')) %>% column_to_rownames('X')
seurat_spatialObj_spacel = seurat_spatialObj
seurat_spatialObj_spacel$'spacel' = 'unmapped'
seurat_spatialObj_spacel@meta.data[rownames(spacel), 'spacel'] = spacel$spacel_pre
spacel_score_ids = colnames(spacel)[grep(colnames(spacel), pattern='Subclo_')]
spacel_score_ids = c(spacel_score_ids, 'Epithelial_normal')
seurat_spatialObj_spacel@meta.data = merge(seurat_spatialObj_spacel@meta.data %>% rownames_to_column('bin_id'), spacel[,c(spacel_score_ids, 'spacel_pre')] %>% rownames_to_column('bin_id'), by='bin_id', all.x=T) %>% column_to_rownames('bin_id')
seurat_spatialObj_spacel@meta.data[,spacel_score_ids][is.na(seurat_spatialObj_spacel@meta.data[,spacel_score_ids])] = 0

Idents(seurat_spatialObj_spacel) <- seurat_spatialObj_spacel$spacel_pre
seurat_spatialObj_spacel$subTypes_predicted.id <- seurat_spatialObj_spacel$spacel
save(seurat_spatialObj2_infercnv, file = paste0(workdir2,"/1.4.10_seurat_spatialObj2_infercnv.RData"))
options(bitmapType = "cairo")
for(i in spacel_score_ids){
    print(i)
    #tiff(file = paste0(workdir, "/1.4.13_spacel2space_infercnv", i, ".tiff"), res=300, width=9, height=8, compression="lzw", units="in")
    pdf(file = paste0(workdir, "/1.4.13_spacel2space_infercnv", i, ".pdf"),width = 9, height = 8)
    tryCatch(print(SpatialFeaturePlot(seurat_spatialObj_spacel, features = i, pt.size.factor = 1.3, crop = TRUE, alpha = c(0.1, 1), stroke = 0)), error = function(e){
        print(c(i,e))
    })
    dev.off()
}
