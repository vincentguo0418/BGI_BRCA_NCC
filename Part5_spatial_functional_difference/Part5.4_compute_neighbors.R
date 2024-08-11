options(stringsAsFactors = FALSE)
library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(Cairo)
library(future)
library(parallel)
library(parallelDist)
options(bitmapType = "cairo")

##
prefix <- "P41"
outdir <- "all_invasive/"
setwd(outdir)
st <- readRDS('SS200000495BR_A1.rds')
newmeta <- read.csv('cellbin_metadata_annotated_v1.txt', header=TRUE, sep='\t', row.names = 1)
print(head(newmeta))
st@meta.data <- newmeta
str(st)

### decide the distance of nerghbors
dis_max <- 50
core=5

print("confirm x and y")
st@meta.data$y <- st@images$slice1@coordinates$imagerow
st@meta.data$x <- st@images$slice1@coordinates$imagecol
print("cell count matrix")
cell_matrix <-  model.matrix(~0+st@meta.data[,"annotated_cluster"]) %>%as.data.frame()
cellnames <- paste0(st@meta.data[,"annotated_cluster"]%>% table() %>% names(),"_bin")
names(cell_matrix) <- cellnames
#cellnames <- cellnames[-8]

#cell_matrix$immune_bin <- cell_matrix$Plasmocyte_bin+cell_matrix$T_cell_bin+cell_matrix$Macrophage_bin+cell_matrix$B_cell_bin
#cellnames <-c(cellnames,"immune_bin")
#cellnames <- setdiff(cellnames,c("Plasmocyte_bin","T_cell_bin","Macrophage_bin",'B_cell_bin'))

rownames(cell_matrix) <- rownames(st@meta.data)
st = AddMetaData(st, cell_matrix[,cellnames], col.name = cellnames)
print(head(cell_matrix[,cellnames]))
print(cellnames)

print("calculate cell percentage")
clnum<-20
cell_id <- rownames(st@meta.data)
cellnames <- c("IBC cells_bin")
system.time({
  for(index in cellnames){
    cl <- makeCluster(getOption("cl.cores", clnum))
    clusterExport(cl, "st")
    clusterExport(cl, "cell_id")
    clusterExport(cl, "index")
    clusterExport(cl, "dis_max")
    coef_nei = parSapply(cl,cell_id, function(spot){
      x_range<-c((st@meta.data[spot,"x"]-dis_max):(st@meta.data[spot,"x"]+dis_max))
      y_range<-c((st@meta.data[spot,"y"]-dis_max):(st@meta.data[spot,"y"]+dis_max))
      index_id <- setdiff(rownames(st@meta.data)[st@meta.data$x %in% x_range & st@meta.data$y %in% y_range],spot)
      y = st@meta.data[index_id,index]
      y[y < 0] = 0
      mean(y, na.rm = TRUE)
    })
    stopCluster(cl);
  } 
})

coef_nei %>% summary()

cell_nei <- paste0(cellnames,"_nei")
coef_nei <- as.data.frame(coef_nei)
colnames(coef_nei) = cell_nei
# rownames(coef_nei) = rownames(st@meta.data)
st = AddMetaData(st, coef_nei, col.name = cell_nei)
print("get cell percentage")
print("output RDS")
saveRDS(st,file=paste0(outdir,"/st_Nei.RDS"))
write.table(st@meta.data, file = paste0(outdir,"/",prefix,"_st_nei_meta.txt"), quote = F, sep = "\t",row.names = T,col.names=NA)

# st <- readRDS(paste0(outdir,"/st_Nei.RDS"))
print("Spatial Feature Plot of cell percentage")
for(pre in cell_nei){
  # st@meta.data[,pre] <-st@meta.data[,pre]  %>% expr_unoutlier_alt(up=0.995,low=0.005)
  h = SpatialFeaturePlot(st, features = pre,stroke = 0,crop = TRUE, alpha = c(0.8, 1),pt.size.factor = 9) +
    theme(legend.position = "right")  
  # pdf(paste0(outdir, "/",pre,"_spatial_Nei.pdf"),width =3, height = 3)
  # print(h)
  # dev.off()
  tiff(file = paste0(outdir, "/",pre,"_spatial_Nei.tiff"), res=300, width=7, height=4, compression="lzw", units="in")
  print(h)
  dev.off()
}