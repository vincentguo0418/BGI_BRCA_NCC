library('tidyverse')
library('Seurat')
library('parallel')
library('ggplot2')
library('ggsci')
options(bitmapType = "cairo")

args = commandArgs(T)
sample = args[1]
st_ids = c('P13'='148BR_D1_13-1', 'P15'='148BR_F2_15-1', 'P161'='148TR_E1_16-1', 'P22'='22-1', 'P28'='28-1', 'P29'='29-1', 'P35'='35-1', 'P36'='36-1', 'P41'='41-1')
st_id = st_ids[sample]
indir = paste0('/', st_id, '/Results/8.30_cellbin_bin100_tree')
indir2 = paste0('/', st_id, '/Results/8.23.4_cellbin_merge_cluster_bin100_merge2rds')
workdir=paste0('paper_figures/figures_part1.2')
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

ids_table <- read.csv(paste0(indir, '/8.30_labels_subclone.txt'),header = TRUE,sep='\t')
ids_table$merged_cluster <- rownames(ids_table)

anno_table <- read.csv(paste0(indir2, '/merged_cluster_labels.txt'),header = TRUE,sep='\t')
anno_table1 <- left_join(anno_table, ids_table, by='merged_cluster') %>% drop_na('subclone')

meta <- read.csv(paste0(indir2, '/', sample, '_cellbin_metadata_annotated.txt'),header = TRUE,sep='\t')
meta$cell_id <- meta$X #P22
anno_table1 <- left_join(anno_table1, meta, by='cell_id') #P22

anno_table2 <- filter(anno_table1, anno_table1$'annotated_cluster' %in% c('DCIS cells','IBC cells'))
anno_table2$'Cell_type' <- anno_table2$'annotated_cluster'
anno_table2 <- select(anno_table2, c('subclone','Cell_type'))

anno_table3 <- anno_table2 %>%
  group_by(subclone, Cell_type) %>%
  summarise(n = n()) %>% 
  group_by(subclone) %>%
  mutate(percent = n/sum(n))

##plot distribution barplot(Fig4b,e,h,k)
pdf(paste0(workdir, '/', sample, '_stat.pdf'),width=7,height=7)
ggplot(anno_table3, aes(x = subclone, weight = percent, fill = Cell_type))+
  geom_bar(position = "stack")+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_nejm()+
  labs(y='Cell proportion', x=NULL, title="Distribution of cell types in subclones", fill='Cell type')+
  theme(text = element_text(size = 20),
        panel.background=element_rect(fill="white"),
        plot.title = element_text(hjust=0.5, size = 22),
        legend.position = "right",
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))
dev.off()
