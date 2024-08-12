library(tidyverse)
library(Seurat)
library(MuDataSeurat)
library(SeuratDisk)
library(ggplot2)
library(ggpubr)

##Calculate signature score(for Fig6d,e)
CoStimulatory = list(c('CD27', 'CD28', 'CD40LG', 'ICOS', "SLAMF1", 'TNFRSF14', 'TNFRSF18', 'TNFRSF4', 'TNFRSF9', 'ICOSLG', 'TNFRSF8', 'CD226', 'TNFRSF25'))
Cytotoxic = list(c('GNLY', 'GZMA', 'GZMB', 'GZMK', 'IFNG', 'NKG7', 'PRF1', 'CST7', 'TNFSF10', 'CCL4', 'CCL3', 'FASLG', 'CD44'))
Exhaustion = list(c('BTLA', 'CD276', 'CTLA4', 'ENTPD1', 'HAVCR2', 'IDO1', 'KLRC1', 'LAG3', 'LAYN', 'LGALS9', 'LILRB2', 'LILRB4', 'PD1', 'PD-L1', 'PD-L2', 'TDO2', 'TIGIT', 'VSIR', 'PDCD1', 'CXCL13'))
Inhibitory = list(c("BTLA","C10orf54","CD160","CD244","CD274","CTLA4","HAVCR2","LAG3","LAIR1","TIGIT"))

tumor <- ReadH5AD('cellbin_merged_adata.h5ad')
str(tumor)


tumor <- AddModuleScore(tumor,
                          features = CoStimulatory,
                          ctrl = 100,
                          name = "CoStimulatory_score")
tumor <- AddModuleScore(tumor,
                          features = Cytotoxic,
                          ctrl = 100,
                          name = "Cytotoxic_score")
tumor <- AddModuleScore(tumor,
                          features = Exhaustion,
                          ctrl = 100,
                          name = "Exhaustion_score")
tumor <- AddModuleScore(tumor,
                          features = Inhibitory,
                          ctrl = 100,
                          name = "Inhibitory_score")


colnames(tumor@meta.data)[which(colnames(tumor@meta.data) == "CoStimulatory_score1")] <- "CoStimulatory"
colnames(tumor@meta.data)[which(colnames(tumor@meta.data) == "Cytotoxic_score1")] <- "Cytotoxic"
colnames(tumor@meta.data)[which(colnames(tumor@meta.data) == "Exhaustion_score1")] <- "Exhaustion"
colnames(tumor@meta.data)[which(colnames(tumor@meta.data) == "Inhibitory_score1")] <- "Inhibitory"


VlnPlot(tumor, features = c("CoStimulatory"), group.by = "domain", 
            cols = c('#679966','#FEC643','#C643FE','#993333','#FE6943','#43FE69','#437BFE','#B87A3D','#7F6699'), pt.size = 0) +
            xlab("Spatial domain")
VlnPlot(tumor, features = c("Cytotoxic"), group.by = "domain", 
            cols = c('#679966','#FEC643','#C643FE','#993333','#FE6943','#43FE69','#437BFE','#B87A3D','#7F6699'), pt.size = 0) +
            xlab("Spatial domain")
VlnPlot(tumor, features = c("Exhaustion"), group.by = "domain", 
            cols = c('#679966','#FEC643','#C643FE','#993333','#FE6943','#43FE69','#437BFE','#B87A3D','#7F6699'), pt.size = 0) +
            xlab("Spatial domain")
VlnPlot(tumor, features = c("Inhibitory"), group.by = "domain", 
            cols = c('#679966','#FEC643','#C643FE','#993333','#FE6943','#43FE69','#437BFE','#B87A3D','#7F6699'), pt.size = 0) +
            xlab("Spatial domain")
