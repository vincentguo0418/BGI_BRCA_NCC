library(dplyr)
library(openxlsx)
library(limma)
library(GenomicRanges)
library(magrittr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(ggthemes)
library(parallel)
options(bitmapType = "cairo")

args = commandArgs(T)
sample = args[1]
vcf1_path = args[2] #1_filter_2vcf/P13/1079361-10_67.7_high.annovar.hg38_multianno.txt
vcf2_path = args[3] #1_filter_2vcf/P13/1079361-10_69.6_high.annovar.hg38_multianno.txt
common_vcf_path = args[4] #1_filter_2vcf/P13/common_1079361-10_69.6__1079361-10_67.7.annovar.hg38_multianno.txt
workdir=paste0('/2_pyclone/', sample)
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

##Use Pyclone to identify SNV subclones
id1 = strsplit2(vcf1_path, '/')[1,] %>% last() %>% strsplit2(., '_high.annovar.hg38_multianno.txt') %>% .[1,1]
id2 = strsplit2(vcf2_path, '/')[1,] %>% last() %>% strsplit2(., '_high.annovar.hg38_multianno.txt') %>% .[1,1]
tc1 = read.delim(sprintf('WES/0_sequenza/%s/%s/%s_alternative_solutions.txt', sample, id1, id1))
tc1 = tc1$cellularity[which.max(tc1$SLPP)]
tc2 = read.delim(sprintf('WES/0_sequenza/%s/%s/%s_alternative_solutions.txt', sample, id2, id2))
tc2 = tc2$cellularity[which.max(tc2$SLPP)]
system(sprintf('PyClone run_analysis_pipeline --in_files %s %s --working_dir %s --prior total_copy_number --seed 7 --num_iters 5000 --max_clusters 10 --tumour_contents %s %s', paste0(workdir, '/', id1, '.tsv'), paste0(workdir, '/', id2, '.tsv'), workdir, tc1, tc2))

##Find driver genes in each subclone
根据tree找到各个分枝上的driver gene
cosmic = read.delim('cosmic/Cosmic_CancerGeneCensus_v98_GRCh38.tsv')
cosmic = cosmic %>% filter(!(ROLE_IN_CANCER %in% c('', 'fusion')))
loci3 = read.delim(paste0('WES/2_pyclone/', sample, '/tables/loci.tsv'))
vcf2 %>% filter(mutation_id %in% (loci3 %>% filter(cluster_id %in% c(2,8) & sample_id=='1079897-9_476.9') %>% pull('mutation_id')))
cis_merge_dict = vcf1_merge %>% select('mutation_id', 'gene') %>% distinct()
rownames(cis_merge_dict) = cis_merge_dict$mutation_id
ic_merge_dict = vcf2_merge %>% select('mutation_id', 'gene') %>% distinct()
rownames(ic_merge_dict) = ic_merge_dict$mutation_id
nodes = list('P13'=list('cis_A'=c('1079361-10_67.7', 4), 'ic_A'=c('1079361-10_69.6', 4), 'cis_B'=c('1079361-10_67.7', 1), 'ic_B'=c('1079361-10_69.6', 1), 'cis_C'=c('1079361-10_67.7', 2), 'ic_D'=c('1079361-10_69.6', 0)),
             'P15'=list('cis_A'=c('1079872-3_49.7', 0), 'ic_A'=c('1079872-3_64.5', 0), 'cis_B'=c('1079872-3_49.7', 1), 'ic_C'=c('1079872-3_64.5', 2)),
             'P16'=list('cis_A'=c('1079897-9_7.6', 3), 'ic_A'=c('1079897-9_476.9', 3), 'cis_B'=c('1079897-9_7.6', 1), 'ic_C'=c('1079897-9_476.9', 0), 'cis_D'=c('1079897-9_7.6', 2), 'ic_D'=c('1079897-9_476.9', 2), 'cis_E'=c('1079897-9_7.6', 8), 'ic_E'=c('1079897-9_476.9', 8)),
             'P22'=list('cis_A'=c('P22_CIS', 1), 'ic_A'=c('P22_IC', 1), 'cis_B'=c('P22_CIS', 0), 'ic_C'=c('P22_IC', 3), 'ic_D'=c('P22_IC', 2)),
             'P28'=list('cis_A'=c('1088262-8_18.8', 2), 'ic_A'=c('1088262-8_55.7', 2), 'cis_B'=c('1088262-8_18.8', 1), 'ic_C'=c('1088262-8_55.7', 0)),
             'P29'=list('cis_A'=c('1088515-12_66.9', 6), 'ic_A'=c('1088515-12_73.6', 6), 'cis_D'=c('1088515-12_66.9', 3), 'cis_E'=c('1088515-12_66.9', 4), 'ic_B'=c('1088515-12_73.6', 5), 'ic_C'=c('1088515-12_73.6', 2)),
             'P35'=list('cis_A'=c('1089424-13_33.4', 2), 'ic_A'=c('1089424-12_25.1', 2), 'cis_B'=c('1089424-13_33.4', 1), 'ic_C'=c('1089424-12_25.1', 0)),
             'P36'=list('cis_A'=c('P36_CIS', 4), 'ic_A'=c('P36_IC', 4), 'cis_B'=c('P36_CIS', 0), 'cis_C'=c('P36_CIS', 3), 'ic_C'=c('P36_IC', 3), 'ic_D'=c('P36_IC', 2)),
             'P41'=list('cis_A'=c('1090325-6_20.9', 0), 'ic_A'=c('1090325-6_130', 0), 'cis_B'=c('1090325-6_20.9', 4), 'ic_B'=c('1090325-6_130', 4), 'ic_C'=c('1090325-6_130', 1), 'cis_D'=c('1090325-6_20.9', 2)))
loci_gene_all = data.frame()
for (cluster in names(nodes[[sample]])) {
    node = nodes[[sample]][[cluster]]
    cluster_i = loci3 %>% filter(sample_id==node[1]&cluster_id==node[2]) %>% select(mutation_id, cellular_prevalence) %>% distinct()
    cluster_i = merge(cluster_i, cis_merge_dict, by='mutation_id') %>% mutate('cluster'=cluster)
    loci_gene_all = rbind(loci_gene_all, cluster_i)
}
write.table(loci_gene_all, paste0(workdir, '/', sample, '_loci_gene_all.txt'), row.names=F, col.names=T, quote=F, sep='\t')
loci_gene_all = read.delim(paste0(workdir, '/', sample, '_loci_gene_all.txt'), sep='\t')
drivers = data.frame()
for (i in 1:nrow(loci_gene_all)) {
    gene = loci_gene_all[i, 'gene']
    for (j in 1:nrow(cosmic)) {
        line = cosmic[j,]
        if ((line$'GENE_SYMBOL' %in% strsplit2(gene, ';')[1,]) | length(intersect(strsplit2(gene, ';')[1,], strsplit2(line$'SYNONYMS', ',')[1,]))!=0) {
            drivers = rbind(drivers, cbind(loci_gene_all[i,], cosmic[j,]))
        }
    }
}
write.table(drivers, paste0(workdir, '/', sample, '_drivers.txt'), row.names=F, col.names=T, quote=F, sep='\t')
#*For P36, there're no driver genes. Find genes where the loci have the higher CCF.
nodes_P36 = list('cis_A'=c('P36_CIS', 4), 'ic_A'=c('P36_IC', 4), 'cis_B'=c('P36_CIS', 0), 'cis_C'=c('P36_CIS', 3), 'ic_C'=c('P36_IC', 3), 'ic_D'=c('P36_IC', 2))
top_ccf_P36 = data.frame()
for (cluster in names(nodes_P36)) {
    node = nodes_P36[[cluster]]
    cluster_i = loci3 %>% filter(sample_id==node[1]&cluster_id==node[2]) %>% select(mutation_id, cellular_prevalence) %>% distinct()
    cluster_i = merge(cluster_i, cis_merge_dict, by='mutation_id') %>% mutate('cluster'=cluster)
    top_ccf_P36 = rbind(top_ccf_P36, cluster_i)
}
top_ccf_P36$'cluster0' = strsplit2(top_ccf_P36$cluster, '_')[,2]
write.table(top_ccf_P36, paste0(workdir, '/', sample, '_top_ccf_P36.txt'), row.names=F, col.names=T, quote=F, sep='\t')

##Heatmap colored by mutation types(Fig2g)
samples = c('P13', 'P15', 'P16', 'P22', 'P28', 'P29', 'P35', 'P36', 'P41')
for (sample in samples) {
  print(sample)
  workdir=paste0('/2_pyclone/', sample)
  setwd(workdir)
  loci3 = read.delim(paste0('WES/2_pyclone/', sample, '/tables/loci.tsv'))
  vcf1_merge = read.delim(paste0(workdir, '/', sample, '_vcf1_merge.txt'), sep='\t')
  vcf2_merge = read.delim(paste0(workdir, '/', sample, '_vcf2_merge.txt'), sep='\t')
  cis_merge_dict = vcf1_merge %>% select('mutation_id', 'gene') %>% distinct()
  rownames(cis_merge_dict) = cis_merge_dict$mutation_id
  ic_merge_dict = vcf2_merge %>% select('mutation_id', 'gene') %>% distinct()
  rownames(ic_merge_dict) = ic_merge_dict$mutation_id
  loci_A = loci3 %>% filter(cluster_id==as.numeric(nodes[[sample]][['cis_A']][2]))
  if (sample=='P16') {
    loci_A = loci3 %>% filter(cluster_id %in% c(as.numeric(nodes[[sample]][['cis_A']][2]), as.numeric(nodes[[sample]][['cis_D']][2]), as.numeric(nodes[[sample]][['cis_E']][2]), as.numeric(nodes[[sample]][['cis_B']][2]), as.numeric(nodes[[sample]][['cis_C']][2])))
  }
  loci_A$'gene' = cis_merge_dict[loci_A$mutation_id, 'gene']
  if (sample=='P16') {
    gene_dict = rbind(cis_merge_dict, ic_merge_dict) %>% distinct()
    loci_A$'gene' = gene_dict[loci_A$mutation_id, 'gene']
  }
  loci_A_filter_genes = loci_A %>% group_by(mutation_id) %>% mutate(sum_CCF=sum(cellular_prevalence)) %>% ungroup() %>% as.data.frame() %>% arrange(desc(sum_CCF)) %>% filter(!grepl(';', gene)) %>% pull(gene) %>% unique()
  loci_A_filter = loci_A %>% filter(gene %in% loci_A_filter_genes)
  vcf12_merge_A = rbind(vcf1_merge %>% filter(mutation_id %in% loci_A_filter$mutation_id), vcf2_merge %>% filter(mutation_id %in% loci_A_filter$mutation_id)) %>% arrange(mutation_id)
  top_locis_vcf = vcf12_merge_A %>% .[,c('chr', 'start', 'end',	'ref',	'alt',	'gene', 'sample_info')]
  write.table(top_locis_vcf, file = paste0(workdir, '/top_locis.avinput'), row.names = F, col.names = F, sep = '\t', quote = F)
  system(sprintf('perl annovar/annotate_variation.pl -out %s/%s_top_loci -build hg38 %s/top_locis.avinput annovar/humandb/', workdir, sample, workdir))
  top_locis = read.delim(paste0(workdir, '/', sample, '_top_loci.exonic_variant_function'), header = F, sep = '\t')
  top_locis_type2 = data.frame('gene'=top_locis$V9, 'type'=top_locis$V2, 'sample_info'=top_locis$V10) %>% distinct() %>% mutate('sample'=sample) %>% filter(type!='synonymous SNV')
  loci_A_filter_genes2 = loci_A_filter_genes[loci_A_filter_genes %in% top_locis_type2$gene] %>% .[1:5]
  if (sample=='P16') {
    loci_A_filter_genes2 = loci_A_filter_genes[loci_A_filter_genes %in% top_locis_type2$gene] %>% .[1:3]
  }
  top_locis_type3 = top_locis_type2 %>% filter(gene %in% loci_A_filter_genes2)
  write.table(top_locis_type3, paste0(workdir, '/', sample, '_top_locis_type3.txt'), sep = '\t', row.names = F, col.names = T, quote = F)
}
top_locis_type = data.frame()
for (sample in samples) {
  print(sample)
  workdir=paste0('/2_pyclone/', sample)
  top_locis_type_i = read.delim(paste0(workdir, '/', sample, '_top_locis_type3.txt'))
  top_locis_type = rbind(top_locis_type, top_locis_type_i)
}
top_locis_type2 = top_locis_type %>% select(-c(sample))

get_driver = function(sample) {
#for (sample in setdiff(samples, 'P36')) {
  print(sample)
  workdir=paste0('/2_pyclone/', sample)
  setwd(workdir)
  vcf1_merge = read.delim(paste0(workdir, '/', sample, '_vcf1_merge.txt'))
  vcf1_only = vcf1_merge %>% filter(sample_info==paste0(sample, '_CIS'))
  vcf2_merge = read.delim(paste0(workdir, '/', sample, '_vcf2_merge.txt'))
  vcf2_only = vcf2_merge %>% filter(sample_info==paste0(sample, '_IC'))
  vcf12_only = rbind(vcf1_only, vcf2_only)
  drivers_all = data.frame()
  for (i in 1:nrow(vcf12_only)) {
      gene = vcf12_only[i, 'gene']
      for (j in 1:nrow(cosmic)) {
          line = cosmic[j,]
          if ((line$'GENE_SYMBOL' %in% strsplit2(gene, ';')[1,]) | length(intersect(strsplit2(gene, ';')[1,], strsplit2(line$'SYNONYMS', ',')[1,]))!=0) {
              drivers_all = rbind(drivers_all, cbind(vcf12_only[i,], cosmic[j,]))
          }
      }
  }
  drivers_vcf = vcf12_only %>% filter(mutation_id %in% drivers_all$mutation_id) %>% .[,c('chr', 'start', 'end',	'ref',	'alt',	'gene', 'sample_info')]
  write.table(drivers_vcf, file = paste0(workdir, '/drivers.avinput'), row.names = F, col.names = F, sep = '\t', quote = F)
  system(sprintf('perl annovar/annotate_variation.pl -out %s/%s -build hg38 %s/drivers.avinput annovar/humandb/', workdir, sample, workdir))
  driver_type = read.delim(paste0(workdir, '/', sample, '.exonic_variant_function'), header = F, sep = '\t')
  driver_type2 = data.frame('gene'=driver_type$V9, 'type'=driver_type$V2, 'sample_info'=driver_type$V10) %>% distinct() %>% mutate('sample'=sample)
  write.table(driver_type2, paste0(workdir, '/', sample, '_driver_type2.txt'), sep = '\t', row.names = F, col.names = T, quote = F)
}
driver_type2_list <- mclapply(setdiff(samples, 'P36'), mc.cores=8, function(j){
  get_driver(j)
})
driver_type = data.frame()
for (sample in setdiff(samples, 'P36')) {
  workdir=paste0('/2_pyclone/', sample)
  driver_type_i = read.delim(paste0(workdir, '/', sample, '_driver_type2.txt'))
  driver_type = rbind(driver_type, driver_type_i)
}
loci_gene_all %>% arrange(desc(cellular_prevalence))

driver_type2 = driver_type %>% select(-c(sample)) %>% filter(type!='synonymous SNV')
driver_type2 = rbind(driver_type2, top_locis_type2) %>% distinct() %>% arrange(gene) %>% tidyr::spread(., sample_info, type, fill='') %>% tibble::column_to_rownames('gene')
sample_info_all = apply(expand.grid(samples, c('CIS', 'IC')), 1, paste, collapse = "_") %>% sort()
driver_type3 = driver_type2
for (sample_info in setdiff(sample_info_all, colnames(driver_type2))) {
  driver_type3 = cbind('', driver_type3)
  colnames(driver_type3)[1] = sample_info
}
driver_type3 = driver_type3[,sample_info_all]

col = c('frameshift deletion'='#008000', 'frameshift insertion'='red', 'nonsynonymous SNV'='blue', 'stopgain'='orange')
alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = "gray92", col = NA))
    },
    # big blue
    'frameshift deletion' = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col['frameshift deletion'], col = NA))
    },
    # big red
    'frameshift insertion' = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col['frameshift insertion'], col = NA))
    },
    # small green
    'nonsynonymous SNV' = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col['nonsynonymous SNV'], col = NA))
    },
    # small green
    'stopgain' = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col['stopgain'], col = NA))
    }
)
heatmap_legend_param = list(title = "Alternations", at = c('frameshift deletion', 'frameshift insertion', 'nonsynonymous SNV', 'stopgain'), 
        labels = c('Frameshift Deletion', 'Frameshift Insertion', 'Nonsynonymous SNV', 'Stopgain'))
column_order = order(colnames(driver_type3))
cicle_colors = c("#437BFE", "#FEC643", "#43FE69", "#FE6943", "#C643FE",
                 "#43D9FE", "#B87A3D", "#679966", "#993333")
names(cicle_colors) = unique(strsplit2(colnames(driver_type3), '_')[,1])
bottom_anno = HeatmapAnnotation(df = data.frame(sample = strsplit2(colnames(driver_type3), '_')[,1]),
    col = list(sample = cicle_colors))
P1 = oncoPrint(driver_type3,
    alter_fun = alter_fun, col = col, 
    heatmap_legend_param = heatmap_legend_param, show_column_names = TRUE, column_names_gp = gpar(fontsize = 8), column_names_rot = 45, column_order = column_order, bottom_annotation = bottom_anno)
tiff(file = paste0('/heatmap.tiff'), res=300, width=9, height = 12, compression="lzw", units="in")
print(P1)
dev.off()
pdf(file = paste0('/heatmap.pdf'), width = 9, height = 12)
print(P1)
dev.off()

##Clonal and subclonal mutation statistics(Fig2j)
workdir=paste0('/2_pyclone')
setwd(workdir)
xls_path = 'WES/pyclone_cluster.xlsx'
sheets = loadWorkbook(xls_path)$sheet_names
clone_nums = data.frame()
for (sheet in sheets) {
  print(sheet)
  sheet_i = read.xlsx(xls_path, sheet = sheet)
  sheet_i = sheet_i %>% filter(!is.na(cluster))
  sample = strsplit2(sheet, '_')[1,1]
  sheet_i$'patient' = sample
  sheet_i2 = sheet_i %>% select(cluster, size) %>% distinct() %>% arrange(cluster)
  clone_nums = rbind(clone_nums, c(sample, 'Clonal mutation', sheet_i2 %>% filter(cluster=='A') %>% pull(size)))
  clone_nums = rbind(clone_nums, c(sample, 'Subclonal mutation', sheet_i2 %>% filter(cluster!='A') %>% pull(size) %>% sum()))
}
colnames(clone_nums) = c('patient', 'if_sub', 'num_mut')
write.table(clone_nums, paste0(workdir, '/2_clone_nums.txt'), sep = '\t', quote = F, col.names = T, row.names = F)

clone_nums_cis = read.xlsx(paste0(workdir, '/2_clone_nums_byZone.xlsx'), sheet=1)
P1=ggplot(clone_nums_cis, aes(reorder(patient, num_mut, decreasing=T), num_mut, fill=if_sub)) +
  geom_histogram(position = "stack", stat = "identity") +
  labs(x='Patient', y='Count') +
  theme_tufte() +
  theme(panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(),
        axis.title.y = element_text(),
        legend.title = element_blank(),
        legend.position = c(0.78,0.8)
  ) +
  scale_fill_manual(values = c('Subclonal mutation'="#3C5488", 'Clonal mutation'="#00A087"))
P2=ggplot(clone_nums_cis, aes(reorder(patient, num_mut, decreasing=T), num_mut, fill=if_sub)) +
  geom_histogram(position = "fill", stat = "identity") +
  labs(x='Patient', y='Proportion') +
  theme_tufte() +
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.title.y = element_text(),
        legend.title = element_blank(),
        legend.position = 'none'
  ) +
  scale_fill_manual(values = c('Subclonal mutation'="#3C5488", 'Clonal mutation'="#00A087"))
P12_CIS = cowplot::plot_grid(P1, P2, 
                   ncol = 1, 
                   axis = "lb",
                   align = "v")
pdf(file = paste0(workdir, '/2_clone_mut_num_CIS.pdf'), width = 3, height = 4)
print(P12_CIS)
dev.off()

clone_nums_ic = read.xlsx(paste0(workdir, '/2_clone_nums_byZone.xlsx'), sheet=2)
P1=ggplot(clone_nums_ic, aes(reorder(patient, num_mut, decreasing=T), num_mut, fill=if_sub)) +
  geom_histogram(position = "stack", stat = "identity") +
  labs(x='Patient', y='Count') +
  theme_tufte() +
  theme(panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(),
        axis.title.y = element_text(),
        legend.title = element_blank(),
        legend.position = c(0.78,0.8)
  ) +
  scale_fill_manual(values = c('Subclonal mutation'="#3C5488", 'Clonal mutation'="#00A087"))
P2=ggplot(clone_nums_ic, aes(reorder(patient, num_mut, decreasing=T), num_mut, fill=if_sub)) +
  geom_histogram(position = "fill", stat = "identity") +
  labs(x='Patient', y='Proportion') +
  theme_tufte() +
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.title.y = element_text(),
        legend.title = element_blank(),
        legend.position = 'none'
  ) +
  scale_fill_manual(values = c('Subclonal mutation'="#3C5488", 'Clonal mutation'="#00A087"))
P12_IC = cowplot::plot_grid(P1, P2, 
                   ncol = 1, 
                   axis = "lb",
                   align = "v")
pdf(file = paste0(workdir, '/2_clone_mut_num_IC.pdf'), width = 3, height = 4)
print(P12_IC)
dev.off()
