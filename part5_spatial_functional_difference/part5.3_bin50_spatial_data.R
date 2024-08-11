library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(ggplot2)
library(ggsci)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ggpubr)
library(SingleCellExperiment)
library(RcppML)
library(ggplot2)
library(Cairo)
library(getopt)
options(bitmapType = "cairo")
options(RcppML.threads = 3)

source("scRNA_primary.R")
source("colors.R")
source("signatureGenes.R")

##Create bin50 Seurat object(for Fig5c)
spec <- matrix(c(
    "help", "h", 0, "logical", "help document",
    "infile", "i", 1, "character", "input file,must bt gem file or gem.gz file",
    "outdir", "o", 1, "character", "output dierectory",
    "binsize", "b", 1, "integer", "BIN size for spatial rna-seq,default=50",
    "id", "d", 1, "character", "project id",
    "pc", "p", 1, "integer", "dimension usage in umap and tsne,default=30",
    "res", "r", 1, "double", "resolution of cluster,default=1.0",
    "tumor", "t", "1", "character", "datatype of sample,[PDAC,BRCA]",
    "nGenes", "n", "1", "character", "differential gene number,default=2000"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)
if (!is.null(opt$help) || is.null(opt$infile) || is.null(opt$outdir)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}
infile <- opt$infile
outdir <- opt$outdir
if (is.null(opt$binsize)) {
    opt$binsize <- 50
}


if (is.null(opt$id)) {
    opt$id <- tail(unlist(strsplit(infile, "/", fixed = TRUE)), 1)
    opt$id <- gsub("mat.|.txt|.tsv|.gz|_filtered|.gem", "", opt$id)
}
if (is.null(opt$nGenes)) {
    opt$nGenes <- 2000
}
if (is.null(opt$pc)) {
    opt$pc <- 30
}
if (is.null(opt$res)) {
    opt$res <- 0.6
}
pc <- opt$pc
bs <- opt$binsize
# setwd("/data/public/gongchanghao/script/spatial_pipeline/test/")
organs <- opt$tumor

dir.create(outdir)
setwd(outdir)

pro <- opt$id

############################## 1. bin data  ##############################
dat <- fread(file = infile)
if (grep("MIDCount|MIDCounts", colnames(dat)) > 0) {
    colnames(dat) <- gsub("MIDCount|MIDCounts", "UMICount", colnames(dat))
}
out <- as.data.frame(dat)

dat$x <- trunc((dat$x - min(dat$x)) / bs + 1)
dat$y <- trunc((dat$y - min(dat$y)) / bs + 1)

out <- cbind(dat$y, dat$x, out)
colnames(out)[1:2] <- c(paste0("bin", bs, ".y"), paste0("bin", bs, ".x"))

fwrite(out, paste0(pro, "_bins", bs, "_information.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

dat <- dat[, sum(UMICount), by = .(geneID, x, y)]
dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
bin.coor <- dat[, sum(V1), by = .(x, y)]

out <- as.data.frame(cbind(paste0("BIN.", unique(dat$bin_ID)), bin.coor$y, bin.coor$x))
colnames(out) <- c(paste0("BIN.", bs), paste0("bin", bs, ".y"), paste0("bin", bs, ".x"))
rownames(out) <- out[, 1]
fwrite(out, paste0(pro, "_bin", bs, "_position.txt"),
    col.names = T, row.names = F, sep = "\t", quote = FALSE
)

##
geneID <- seq(length(unique(dat$geneID))) ## 36249 detected genes
hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID)
gen <- hash.G[dat$geneID, "values"]

##
bin_ID <- unique(dat$bin_ID)
hash.B <- data.frame(row.names = sprintf("%d", bin_ID), values = bin_ID)
bin <- hash.B[sprintf("%d", dat$bin_ID), "values"]

##
cnt <- dat$V1

rm(dat)
gc()

##
tissue_lowres_image <- matrix(1, max(bin.coor$y), max(bin.coor$x))
tissue_positions_list <- data.frame(
    row.names = paste("BIN", rownames(hash.B), sep = "."),
    tissue = 1,
    row = bin.coor$y,
    col = bin.coor$x,
    imagerow = bin.coor$y,
    imagecol = bin.coor$x
)



scalefactors_json <- toJSON(list(
    fiducial_diameter_fullres = 1,
    tissue_hires_scalef = 1,
    tissue_lowres_scalef = 1
))




##
mat <- sparseMatrix(i = gen, j = bin, x = cnt)
rownames(mat) <- rownames(hash.G)
colnames(mat) <- paste("BIN", sprintf("%d", seq(max(hash.B[, "values"]))), sep = ".")


############################## 2. creat Spatial Object  ##############################
seurat_spatialObj <- CreateSeuratObject(mat, project = "Spatial", assay = "Spatial", min.cells = 1, min.features = 1)
generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
    if (filter.matrix) {
        tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
    }

    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef

    spot.radius <- unnormalized.radius / max(dim(image))

    return(new(
        Class = "VisiumV1",
        image = image,
        scale.factors = scalefactors(
            spot = scale.factors$tissue_hires_scalef,
            fiducial = scale.factors$fiducial_diameter_fullres,
            hires = scale.factors$tissue_hires_scalef,
            lowres = scale.factors$tissue_lowres_scalef
        ),
        coordinates = tissue.positions,
        spot.radius = spot.radius
    ))
}

spatialObj <- generate_spatialObj(
    image = tissue_lowres_image,
    scale.factors = fromJSON(scalefactors_json),
    tissue.positions = tissue_positions_list
)

##
spatialObj <- spatialObj[Cells(seurat_spatialObj)]
DefaultAssay(spatialObj) <- "Spatial"

seurat_spatialObj[["slice1"]] <- spatialObj

##
rm(mat)
rm(bin.coor)
rm(hash.G)
rm(hash.B)
rm(bin)
rm(gen)
rm(cnt)

dir.create("figures")

##############################  3. Spatial Analyse  ##############################
seurat_spatialObj[["percent.mt"]] <- PercentageFeatureSet(seurat_spatialObj, pattern = "^MT-")

##
Q1 <- quantile(seurat_spatialObj$nFeature_Spatial)[2]
Q3 <- quantile(seurat_spatialObj$nFeature_Spatial)[4]
upper <- as.numeric(Q3 + 1.5 * (Q3 - Q1))
lower <- as.numeric(Q1 - 1.5 * (Q3 - Q1))

save(list = ls(), file = "tmpfiles.RData")

pdf(paste0("figures/", pro, "_bin", bs, "_preQC.pdf"), width = 10, height = 8)
p1 <- VlnPlot(seurat_spatialObj,
    features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), group.by = "orig.ident",
    ncol = 4, pt.size = 0
) +
    theme(axis.text.x = element_text(angle = 0, size = 0), axis.title.x = element_text(angle = 20, size = 8)) +
    labs(x = paste0("nGene:", dim(seurat_spatialObj)[1], "; ", "nBIN:", dim(seurat_spatialObj)[2]))
print(p1)

p2 <- ggplot(seurat_spatialObj@meta.data, aes(x = nFeature_Spatial)) +
    geom_density(colour = "black") +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold.italic"), legend.position = "none",
        axis.title = element_text(size = 15, face = "bold.italic"), axis.text.x = element_text(size = 12), axis.ticks.x = element_blank()
    ) +
    geom_vline(aes(xintercept = 100, colour = "#999999", linetype = "twodash")) +
    geom_vline(aes(xintercept = 200, colour = "#999999", linetype = "twodash")) +
    geom_vline(aes(xintercept = 300, colour = "#999999", linetype = "twodash")) +
    geom_vline(aes(xintercept = 500, colour = "#999999", linetype = "twodash")) +
    geom_vline(aes(xintercept = lower, colour = "#377EB8", linetype = "twodash")) +
    geom_vline(aes(xintercept = upper, colour = "#E41A1C", linetype = "twodash")) +
    xlim(min(seurat_spatialObj@meta.data$nFeature_Spatial), max(seurat_spatialObj@meta.data$nFeature_Spatial)) +
    ggtitle(paste0(pro, ".nBIN_", bs, ":", dim(seurat_spatialObj@meta.data)[1]))
print(p2)

dev.off()

pdf(paste0("figures/", pro, "_bin", bs, "_spatial_dis.pdf"))

SpatialFeaturePlot(seurat_spatialObj, features = "nFeature_Spatial", stroke = 0) +
    theme(legend.position = "right")
SpatialFeaturePlot(seurat_spatialObj, features = "nCount_Spatial", stroke = 0) +
    theme(legend.position = "right")

dev.off()

seurat_spatialObj <- NormalizeData(seurat_spatialObj, assay = "Spatial", verbose = FALSE)
removeBiasGenes <- function(mat){

  RPgenes <- rownames(mat)[intersect(grep("^RP", rownames(mat)), grep("-", rownames(mat)))]
  RPgenes2 <- rownames(mat)[grep("^RP[SL]", rownames(mat))]
  MTgenes <- rownames(mat)[grep("^MT-", rownames(mat))]
  CTCgenes <- rownames(mat)[intersect(grep("^CTC", rownames(mat)), grep("-", rownames(mat)))]
  MIRgenes <- rownames(mat)[grep("^MIR", rownames(mat))]
  ACgenes <- rownames(mat)[intersect(grep("^AC[0-9]", rownames(mat)), grep(".", rownames(mat)))]
  CTgenes <- rownames(mat)[intersect(grep("^CT", rownames(mat)), grep("-", rownames(mat)))]
  LINCgenes <- rownames(mat)[grep("^LINC[0-9]", rownames(mat))]
  ALgenes <- rownames(mat)[intersect(grep("^AL", rownames(mat)), grep(".", rownames(mat)))]

  rmgenes <- c(RPgenes, RPgenes2, MTgenes, CTCgenes, MIRgenes, ACgenes, CTgenes, LINCgenes, ALgenes)

  # datacount <- mat[!rownames(mat)%in%rmgenes,]
  # datacount <- datacount[rowSums(datacount > 0) > 1,]
  return(rmgenes)
}
rmgenes <- removeBiasGenes(seurat_spatialObj)
filetered_genes <- rownames(seurat_spatialObj)[!rownames(seurat_spatialObj)%in%rmgenes]
seurat_spatialObj <- subset(seurat_spatialObj, features = filetered_genes)
seurat_spatialObj <- subset(seurat_spatialObj, subset = nFeature_Spatial > 100)
seurat_spatialObj <- FindVariableFeatures(seurat_spatialObj, selection.method = "vst", nfeatures = opt$nGenes)
seurat_spatialObj <- ScaleData(seurat_spatialObj)
seurat_spatialObj <- RunPCA(seurat_spatialObj, features = VariableFeatures(object = seurat_spatialObj), verbose = F)
save(seurat_spatialObj, file = "seurat_spatialObj.RData")



##############################  4. Gene expression visualization  ##############################
genes <- c('SPARC','TOP2A','ARAF')
output_file <- paste0(workdir, "genes.tif")
tiff(file = output_file, width = 10, height = 10, units = 'in', res = 300)
SpatialFeaturePlot(seurat_spatialObj, features = genes, crop = TRUE, stroke = 0)
dev.off()