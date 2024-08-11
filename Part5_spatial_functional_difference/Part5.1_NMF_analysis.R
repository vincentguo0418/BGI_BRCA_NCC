library(dplyr)
library(corrplot)
library(methods)
library(ggplot2)
library(cowplot)
library(msigdbr)
library(aplot)
library(stringr)
library(clusterProfiler)
library(NbClust)
library(magrittr)
library(factoextra)
library(purrr)
setwd("test2/")

##NMF factorization(for Fig5a)
# Read all w_matrices and merge them by taking the intersection of genes
spatiallist <- list.files("BRCA_binbase/", "nmf_w_matrix.xls", full.names = TRUE, recursive = TRUE)
merged_data <- NULL
for (fl in spatiallist) {
    tmp <- read.table(fl, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
    prefix <- str_extract(fl,'Patient\\w+')
    colnames(tmp) <- paste0(prefix, ".", colnames(tmp))
    if (is.null(merged_data)) {
        merged_data <- tmp
    } else {
        merged_data <- merge(merged_data, tmp, by = "row.names")
        rownames(merged_data) <- merged_data[, 1]
        merged_data <- merged_data[, -1]
    }
}
write.table(merged_data,"merged_data.xls",sep='\t')
merged_data_filtered <- merged_data[rowSums(merged_data) > 0,]
write.table(merged_data_filtered,"merged_data_filtered.xls",sep='\t')
col=colorRampPalette(c("blue","white", "red"))
cor_mat = cor(merged_data_filtered, method = "pearson")
# Assign the lower triangular matrix of the correlation coefficient matrix to NA
cor_mat[lower.tri(cor_mat)] <- NA
# Find columns with correlation coefficients greater than 0.8
col_index <- which(cor_mat > 0.8, arr.ind = TRUE)
col_index <- col_index[order(col_index[,2]),]

# Merge columns with high correlation into a new column, where the new column is equal to the PC1 of these columns
# Extract columns with high correlation and merge them into a new table
data_modules <- NULL
single_modules <- NULL
for (i in unique(col_index[, 2])) {
  ind <- col_index[, 2] == i
  if (sum(ind) > 1) {
    new_col_name <- paste0(colnames(merged_data_filtered)[i], "_", 1:sum(ind))
    # new_col_values <- rowMeans(merged_data[, col_index[ind, 1]], na.rm = TRUE)
    new_col_values <- merged_data_filtered[, col_index[ind, 1]]
    if (is.null(data_modules)){
        data_modules <- new_col_values
    }else{
        data_modules <- cbind(data_modules, new_col_values)
    }
    # colnames(data_modules)[ncol(data_modules)] <- new_col_name
  }
}
data_modules <- data_modules[,!duplicated(colnames(data_modules))]
data_modules <- data_modules[rowSums(data_modules) > 0,]

# Calculate correlation
cor_mat_module <- cor(data_modules, method = "pearson")

# Cluster according to correlation and determine the optimal number of clusters based on NbClust
hc_cormat_order <- function(cormat, hc.method = "complete") {
    dd <- stats::as.dist((1 - cormat) / 2)
    hc <- stats::hclust(dd, method = hc.method)
    hc
}
hc <- hc_cormat_order(cor_mat_module)
y <- NbClust(cor_mat_module, distance = "euclidean", min.nc = 2, max.nc = 40 , method = "complete", index = "kl", alphaBeale = 0.1)
k <- y$Best.nc[[1]] # best cluster number
cluster_modules <- cutree(hc, k = k)
length(cluster_modules)
M.clust <- cbind(cor_mat_module, cluster_modules)
write.table(M.clust, file = "heat_cluster.xls", quote = F, sep = "\t")
write.table(cor_mat_module, file = "cor_mat_module.xls", quote = F, sep = "\t",col.names = NA)

# Plot 1: correlation plot
source("ggcorrplot.R")
source("colors.R")
p1 <- ggcorrplot(cor_mat_module, hc.order = FALSE, outline.col = "white",lab = FALSE) +
    scale_y_discrete(position = "right") +
    scale_x_discrete(labels = NULL)
legend_p1 <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")
lim <- xlim2(p1)
lim$limits <- rownames(cor_mat_module)[hc$order]
p1 <- p1 + lim
lim$axis <- 'y'
p1 <- p1 + lim

# Plot 2: the color annotation of samples
data_anno <- gsub("\\.nmf.*", "", rownames(cor_mat_module))
data_anno <- data.frame(data = rownames(cor_mat_module), anno = as.factor(data_anno))
data_anno_color <- getDefaultColors(length(unique(data_anno$anno)),type=1)
p2 <- ggplot(data_anno, aes(x = data, fill = anno, y = 0.5)) +
    geom_tile() +
    theme_void() + scale_fill_manual(values = data_anno_color) +
    xlim2(p1)
legend_p2 <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")
# Plot 3: the color annotation of clusters
annotation_col <- data.frame(data = names(cluster_modules),cluster =as.character(cluster_modules))
annotation_col_color <- getDefaultColors(length(unique(annotation_col$cluster)),type=2)
p22 <- ggplot(annotation_col, aes(x = data, fill = cluster, y = 0.5)) +
    geom_tile() +
    theme_void() + scale_fill_manual(values = annotation_col_color) +
    xlim2(p1)
legend_p22 <- get_legend(p22)
p22 <- p22 + theme(legend.position = "none")

# Plot 4: The HALLMARK enrichment result of NMF modules
m_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol)
m_df$gs_name <- gsub("HALLMARK_", "", m_df$gs_name)
gc_hyper <- list()
for (i in c(1:ncol(data_modules))) {
    tmp <- data_modules[, i]
    names(tmp) <- rownames(data_modules)
    tmp <- tmp[order(tmp, decreasing = TRUE)]
    tmp <- tmp[tmp > 0]
    tmp_n <- tmp[1:100]
    tmp <- list(tmp)
    names(tmp)[1] <- i
    tmp_n <- names(tmp_n)
    tmp_n <- list(tmp_n)
    names(tmp_n)[1] <- colnames(data_modules)[i]
    gc_hyper <- append(gc_hyper, tmp_n)
}
ehm <- compareCluster(gc_hyper, fun = "enricher", TERM2GENE = m_df, pvalueCutoff = 0.05)
write.table(as.data.frame(ehm), "compare_hyper.txt", sep = "\t", row.names = FALSE)
p3 <- dotplot(ehm)
p3$data$Cluster <- gsub("\n.*", "", p3$data$Cluster)
p3 <- p3 + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), ) + xlim2(p1)
legend_p3 <- get_legend(p3)
p3 <- p3 + theme(legend.position = "none")
pdf("hallmarker_enrichment_cluster.pdf", width = 55, height = 14)
print(p3)
dev.off()

library(org.Hs.eg.db)
ego <- compareCluster(gc_hyper, fun = "enrichGO", keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05, OrgDb = org.Hs.eg.db, universe = rownames(data_modules))
p4 <- dotplot(ego)
p4$data$Cluster <- gsub("\n.*", "", p3$data$Cluster)
p4 <- p4 + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), ) + xlim2(p1)
legend_p4 <- get_legend(p4)
p4 <- p4 + theme(legend.position = "none")
pdf("GO_BP_enrichment_cluster.pdf", width = 55, height = 14)
print(p4)
dev.off()


# Merge
pdf("corrplot_pearson_ggplot_HALLMARK2.pdf", width = 40, height = 47)
pv <- plot_grid(p1, p2, p22, p3, nrow = 4, align = "v", rel_heights = c(0.8, 0.01,0.01, 0.3))
ph <- plot_grid(legend_p1, legend_p2, legend_p22, legend_p3, ncol = 1)
p <- plot_grid(pv, ph, rel_widths = c(1, 0.05))
print(p)
dev.off()

# Obtain PCA analysis results for each cluster
gene_module <- NULL
for (cluster in unique(cluster_modules)){
    nmfs <- names(cluster_modules[cluster_modules==cluster])
    tmp <- data_modules[,nmfs]
    pca_result <- princomp(tmp, cor = TRUE)
    pc1 <- pca_result$scores[, 1]
    pc1_norm <- (pc1 + abs(min(pc1)) + 1) / sum(pc1 + abs(min(pc1)) + 1)
    if (is.null(gene_module)){
        gene_module <- pc1_norm
    }
    else{
        gene_module <- cbind(gene_module,pc1_norm)
    }
}
colnames(gene_module) <- unique(cluster_modules)

# Extract the top 100 gene lists for each gene module and perform HALLMARK and GO enrichment
gc_hyper <- list()
gene_list_top100 <- NULL
for (i in c(1:ncol(gene_module))) {
    tmp <- gene_module[, i]
    names(tmp) <- rownames(gene_module)
    tmp <- tmp[order(tmp, decreasing = TRUE)]
    tmp <- tmp[tmp > 0]
    tmp_n <- tmp[1:100]
    gene_list_top100 <- cbind(gene_list_top100,names(tmp_n))
    tmp <- list(tmp)
    names(tmp)[1] <- i
    tmp_n <- names(tmp_n)
    tmp_n <- list(tmp_n)
    names(tmp_n)[1] <- colnames(gene_module)[i]
    gc_hyper <- append(gc_hyper, tmp_n)
}
colnames(gene_list_top100) <- colnames(gene_module)
write.table(gene_module, file = "gene_module.xls", quote = F, sep = "\t",col.names = NA)
write.table(gene_list_top100,"gene_list_top100.xls",sep='\t',col.names=NA)

ehm <- compareCluster(gc_hyper, fun = "enricher", TERM2GENE = m_df, pvalueCutoff = 0.05)
write.table(as.data.frame(ehm), "gene_module_hyper_HALLMARK.txt", sep = "\t", row.names = FALSE)
p <- dotplot(ehm)
p$data$Cluster <- gsub("\n.*", "", p$data$Cluster)
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), )
pdf("gene_module_hallmark_enrichment_cluster.pdf", width = 10, height = 7)
print(p)
dev.off()
ego <- compareCluster(gc_hyper, fun = "enrichGO", keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05, OrgDb = org.Hs.eg.db, universe = rownames(cluster_modules))
write.table(as.data.frame(ego), "gene_module_hyper_GOBP.txt", sep = "\t", row.names = FALSE)
p <- dotplot(ego)
p$data$Cluster <- gsub("\n.*", "", p$data$Cluster)
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), ) + xlim2(p1)
pdf("gene_module_GO_BP_enrichment.pdf", width = 55, height = 14)
print(p)
dev.off()



#### For each cellbin matrix, calculate the cell weights for each gene module

library(Seurat)
library(RcppML)
library(ComplexHeatmap)
gene_module <- read.delim("gene_module.xls",row.names=1)
predict.nmf <- function(w, data, L1 = 0, L2 = 0, mask = NULL, ...){
  m <- new("nmf", w = w, d = rep(1:ncol(w)), h = matrix(0, nrow = ncol(w), 1))
  predict(m, data, L1 = L1, L2 = L2, mask = mask, ...)
}
x <- load("SS200000495BR_A1_metadata.RData")
seurat_obj <- get(x)
od <- "binbase_NMF"
setwd(od)
#seurat_obj <- NormalizeData(seurat_obj, assay = "RNA", verbose = FALSE)
data <- GetAssayData(seurat_obj, "counts")
w <- as.matrix(gene_module)
data_filtered <- data[rownames(data) %in% rownames(w),]
w <- w[rownames(w) %in% rownames(data_filtered),]
h_matrix_new <- predict.nmf(w=w,data=data_filtered)
rownames(h_matrix_new) <- colnames(w)
h_matrix_new_normed <- LogNormalize(data = h_matrix_new)
write.table(h_matrix_new,paste0(od,"/h_matrix_new.txt"),sep = '\t', row.names = TRUE, col.names = NA)
seurat_obj <- AddMetaData(seurat_obj,t(h_matrix_new_normed),col.name = rownames(h_matrix_new_normed))
seurat_obj[["nmf"]] <- CreateDimReducObject(embeddings = as.matrix(t(h_matrix_new_normed)), loadings = w, key = "NMF_", assay = DefaultAssay(seurat_obj))
saveRDS(seurat_obj,file=paste0(od,"/cellbin_cluster.rds"))
  
nmf_markers <- rownames(h_matrix_new_normed)
#pdf("Vlnplot_nmfs.pdf",height=10,width=8)
#print(VlnPlot(seurat_obj,features=paste0('nmf',c(1:4)),group.by='cellsubtype',pt.size =0,flip = TRUE, stack = TRUE, log = FALSE))
#dev.off()
#pdf("Vlnplot_nmfs.pdf",height=10,width=8)
#print(VlnPlot(seurat_obj,features='GM1',group.by='res.1.4',pt.size =0,y.max=5))
#dev.off()
#hei <- ceiling(length(nmf_markers) / 4) * 3
#pdf(file = paste0(od,"/nmf_markers_spatial_featureplots.pdf"), width = 14.5, height = hei)
#print(SpatialFeaturePlot(seurat_obj, features = nmf_markers, ncol = 4, stroke = 0, pt.size.factor = 30))
#dev.off()
metadata <- read.delim("cellbin_metadata_magic.txt",row.names=1)
# cellbin plot
for (nmf in nmf_markers){
  # exp <- ceiling((as.data.frame(data_filtered))*1000000)
  metadata_seurat <- seurat_obj@meta.data
  print(metadata_seurat)
  plotdata <- cbind(metadata[,c("x","y")],scale(metadata_seurat[nmf]))
  colnames(plotdata)=c("bin1x","bin1y","exp")
  # Remove extreme values
  score_quantile <- quantile(plotdata$exp,c(0.95))
  plotdata$exp[plotdata$exp>score_quantile[[1]]]=score_quantile[[1]]
  genefilename=paste0(od,"/",nmf,"_exp.txt")
  write.table(plotdata,genefilename,sep="\t",quote=F,row.names = F)
  collen=length(unique(plotdata$exp))
  colors=colorRampPalette(c("#4B0082","#9E0142","#FFFF00"))(collen)
  sortexp=sort(as.numeric(unique(plotdata$exp)))
  collist=cbind(sortexp,colors)
  colfilename=paste0(od,"/",nmf,"_col.list")
  write.table(collist,file=colfilename,quote=F,sep="\t",row.names=F,col.names=F)
  tifname=paste0(od,"/",nmf,"_exp.tif")
  mask <- "SS200000495BR_A1_regist_mask.tif"
  shell=paste0("cell_bin_plot  ", genefilename, " ", mask, " ", colfilename, " ", tifname)
  grep_out<-system(shell, intern = TRUE)
  ## Generate a legend for each image      
  col_fun = circlize::colorRamp2(c(-1, 0, 1), c("#4B0082", "#9E0142", "#FFFF00"))
  lgd <- Legend(col_fun = col_fun, title = "score", at = c(min(sortexp),ceiling((max(sortexp)-min(sortexp))/2),max(sortexp)))
  pdf(paste0(od,"/",nmf,"_legend.pdf"))
  draw(lgd)
  dev.off()
}
