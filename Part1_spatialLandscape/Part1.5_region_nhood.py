import os, sys
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
import squidpy as sq

od = 'SPACEL/Splane_cellbin/n_cluster_10_n_neighbors_8_k_2_gnn_dropout_0.8/P13/'
os.chdir(od)

##Calculate cell proportion matrix between domain pairs(Fig1j)
adata = sc.read("SPACEL/Splane_cellbin/n_cluster_10_n_neighbors_8_k_2_gnn_dropout_0.8/cellbin_merged_adata.h5ad")
adata = adata[adata.obs['batch']=='P13']
sc.set_figure_params(facecolor="white", figsize=(4, 4))
sc.settings.verbosity = 3
sc.settings.dpi = 300
sc.settings.figdir = "./figures"

adata.obs['spatial_domain'] = adata.obs['spatial_domain'].astype('str')
adata.obs['spatial_domain'] = adata.obs['spatial_domain'].astype('category')
cluster_annotation = {'0':'MESR','1':'CISR','2':'ICR','3':'IBCR','4':'IBCAR1','5':'UNDR','6':'MIMR','7':'CISAR', '8':'IBCAR2','9':'NMR'}
adata.obs['region'] = adata.obs['spatial_domain'].map(cluster_annotation).astype('category')

sq.gr.spatial_neighbors(adata,coord_type='generic',library_key='batch',radius=30)
sq.gr.nhood_enrichment(adata, cluster_key="region")
score = adata.uns['region_nhood_enrichment']['zscore']
adata.uns['region_nhood_enrichment']['zscore'] = np.nan_to_num(score)

region_number = adata.obs['region'].value_counts()[adata.obs['region'].cat.categories]
nhood_counts = pd.DataFrame(adata.uns['region_nhood_enrichment']['count'], index=adata.obs['region'].cat.categories, columns=adata.obs['region'].cat.categories)
nhood_percents = nhood_counts/region_number
order = ['CISR','CISAR','IBCR','IBCAR1','IBCAR2','MESR','ICR','MIMR','NMR','UNDR']
nhood_percents = nhood_percents[order]
order = [1,7,3,4,8,0,2,6,9,5]
nhood_percents = nhood_percents.iloc[order]

df = nhood_percents
index_list = df.index.tolist()
# Remove the numbers on the diagonal
df = df.mask(np.eye(len(df), dtype=bool))
# the numbers on the diagonal equal to max value in the matrix
df = df.mask(np.eye(len(df), dtype=bool), df.max().max())
print(df)
# zscore df
# df = (df - df.mean()) / df.std()
import seaborn as sns
sns.set_style('whitegrid')
fig,ax = plt.subplots(figsize = (5,4))
sns.heatmap(df, cmap=sns.diverging_palette(50, 250, n=100),center = 0,ax=ax,linewidths=0.3,vmin=0,vmax=0.8) 
ax.set_facecolor('#dcdcdc')
plt.savefig('region_nhood_heatmap.pdf', dpi=300, bbox_inches='tight')           
