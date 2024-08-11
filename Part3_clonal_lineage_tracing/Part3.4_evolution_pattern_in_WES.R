library(dplyr)
library(limma)
library(ggplot2)
library(ggthemes)
options(bitmapType = "cairo")

workdir=paste0('/2_pyclone/', sample)
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

##CCF dotplot(Fig3g,k,o)
loci_gene_all = read.delim(paste0(workdir, '/', sample, '_loci_gene_all.txt'), sep='\t')
muts = unique(loci_gene_all$mutation_id)
spot_df = data.frame('mutation_id'=muts, 'CIS_CCF'=rep(0, length(muts)), 'IC_CCF'=rep(0, length(muts)), 'region'='-', 'cluster'='-')
rownames(spot_df) = spot_df$mutation_id
for (mut in muts) {
  print(mut)
  mut_df = loci_gene_all %>% filter(mutation_id==mut) %>% arrange(cluster)
  clusters = unique(mut_df$cluster)
  spot_df[mut, 'cluster'] = limma::strsplit2(clusters[1], '_')[1,2]
  if (length(clusters)==2) {
    spot_df[mut, 'CIS_CCF'] = mut_df[1, 'cellular_prevalence']
    spot_df[mut, 'IC_CCF'] = mut_df[2, 'cellular_prevalence']
    spot_df[mut, 'region'] = 'CIS&IC shared'
  } else if ((length(clusters)==1)&(substr(clusters, 1, 3)=='cis')) {
    spot_df[mut, 'CIS_CCF'] = mut_df[1, 'cellular_prevalence']
    spot_df[mut, 'region'] = 'CIS private'
  } else if ((length(clusters)==1)&(substr(clusters, 1, 3)=='ic_')) {
    spot_df[mut, 'IC_CCF'] = mut_df[1, 'cellular_prevalence']
    spot_df[mut, 'region'] = 'IC private'
  }
}
colors = c('#f8481c', '#02ccfe', '#ad0afd', '#ffa62b', '#87a922', '#fc86aa')
P1 = ggplot(spot_df, aes(x=CIS_CCF, y=IC_CCF, color=cluster)) +
  geom_point() +
  scale_color_manual(values = colors) +
  theme_bw() +
  labs(x='CIS CCF', y='IC CCF') +
  guides(color = guide_legend(title = "Subclone"))
pdf(paste0(workdir, '/', sample, '_dotplot.pdf'),width=6,height=4)
print(P1)
dev.off()
