library('car')
library('dplyr')
library('ggplot2')

##Calculate cell type distribution across CISR and IBCR of each Cellbin sample(Fig1i)
stat = read.delim('stat.txt', sep='\t')
stat2 = stat
stat2$'Pat_zone' = paste(stat2$Patient, stat2$Zone, sep='_')
stat2$Cell_type = Recode(stat2$Cell_type, "c('DCIS_cells', 'IBC_cells')='Tumor_cells'")
stat2$Cell_type = factor(stat2$Cell_type, levels=c('Atypical_hyperplasia_epithelial_cell', 'B_cells', 'Tumor_cells', 'Endothelial_cells', 'Fibroblasts', 'MEC', 'Macrophages', 'Mast_cells', 'Plasmocytes', 'Normal_epithelial_cells', 'Pericytes', 'T_cells'), ordered=T)
stat2 = stat2 %>% group_by(Pat_zone) %>% mutate(num2=number/sum(number)) %>% ungroup() %>% as.data.frame()
stat2$num2 = ifelse(stat2$Zone=='DCIS', -stat2$num2, stat2$num2)
stat2$Patient = factor(stat2$Patient, levels=sort(unique(stat2$Patient), decreasing=T), ordered=T)
pdf(file = 'Figure1_cellProportion_barplot.pdf',width = 6, height = 3, onefile=F)
stat2_bar = ggplot(stat2, aes(x=Patient, y=num2)) +
            geom_bar(aes(fill=factor(Cell_type)), position='stack', stat = 'identity') +
            coord_flip() +
            scale_fill_manual(values=c("#81357a", "#81d8ae", "#d5492f", "#d3c3a4",
                               "#FFD700", "#6bd155", "#b96f49", "#6d76cb", "#607d3a", 
                               "#ce9d3f", "#00CED1", "#d04dc7")) +
            labs(x=NULL, y='Cell Proportion') +
            theme_bw()
print(stat2_bar)
dev.off()
