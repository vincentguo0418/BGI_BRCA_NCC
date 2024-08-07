library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(readxl)

data <- read_excel("20221130BRCAST_Pathology_Anonymous.xlsx") %>% data.frame()
rownames(data) <- data$Number
data <- data[, -c(1,6,7)]
print(data)

m = apply(data,1,function(x){length(unique(x))})
col_fun = colorRamp2(breaks = seq(0, 1, length.out = 12),
                     colors = brewer.pal(12,"Paired"))
qz = sort(unique(as.character(data)));length(qz)
colors = col_fun(seq(0, 1, length.out = length(qz)))
names(colors) = qz
head(colors)

colors = c("#aed688", "#cf4c8b","#6b42c8", "#7cd5c8", "#ffd900",
          "#674c2a", "#c8b693", "#798234", "#f6a97a", "#c6a5cc",
           "#feb308", "#cb7c77", "#68d359", "#6a7dc9", "#8F999F",
          "#F08080","#FFF0F5")
names(colors) = c("+","-","0","1+","2+",
                  "3+","T1bN0","T1cN0","T1cN2a","T2N0",
                  "T2N1","T1cN1","T1cN1a","T1miN0","NA",
                  "Yes","No")

ER = Legend(labels = c("+","-","NA"), title = "ER", 
             legend_gp = gpar(fill = c("#aed688", "#cf4c8b","#8F999F")))
PR = Legend(labels = c("+","-","NA"), title = "PR", 
             legend_gp = gpar(fill = c("#aed688", "#cf4c8b","#8F999F")))
HER2 = Legend(labels = c("0","1+","2+","3+","NA"), title = "HER2", 
             legend_gp = gpar(fill = c("#6b42c8", "#7cd5c8", "#ffd900","#674c2a","#8F999F")))
Stage = Legend(labels = c("T1bN0","T1cN0","T1cN2a","T2N0","T2N1","T1cN1","T1cN1a","T1miN0","NA"), title = "Stage", 
             legend_gp = gpar(fill = c("#c8b693", "#798234", "#f6a97a", "#c6a5cc", "#feb308", "#cb7c77", "#68d359", "#6a7dc9", "#8F999F")))
scRNA = Legend(labels = c("Yes","No"), title = "scRNA-seq", 
            legend_gp = gpar(fill = c("#F08080","#FFF0F5")))
Stereo = Legend(labels = c("Yes","No"), title = "Stereo-seq", 
                   legend_gp = gpar(fill = c("#F08080","#FFF0F5")))
p = packLegend(list = c(ER,PR,HER2,Stage,scRNA,Stereo),
               direction = "horizontal",gap = unit(0.6, "cm"))

# Color
col_group = data[, 1]
color_an = brewer.pal(length(unique(col_group)),"Dark2")
names(color_an) = unique(col_group)
top_annotation = HeatmapAnnotation(cluster = col_group,col = list(cluster = color_an),show_legend = F,show_annotation_name = F)
#top_annotation = HeatmapAnnotation(text=anno_text(col_group))

pdf(file = "fig1A.pdf",width=18,height=3)
f = Heatmap(t(data[,c(2,3,4,5,6,7)]), col = colors,
          column_split = col_group,
          top_annotation = top_annotation,
          border = T,show_heatmap_legend = F,row_names_side = "left",
          rect_gp = gpar(col = "white", lwd = 1),
          column_names_rot = 45)
draw(f,heatmap_legend_list = p,heatmap_legend_side = "right")
dev.off()

png("fig1A.png",width=1800,height=300,units='px')
draw(f,heatmap_legend_list = p,heatmap_legend_side = "right")
dev.off()
