import os, sys
del sys.path[4]
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
import squidpy as sq

od = 'ST_subclone_region_nhood/'
os.chdir(od)

sample = 'P28'
adata = sc.read("%s_subclone.h5ad" % sample)
sc.set_figure_params(facecolor="white", figsize=(4, 4))
sc.settings.verbosity = 3
sc.settings.dpi = 300
sc.settings.figdir = "./figures"

adata.obs['subclone2'] = adata.obs['subclone2'].astype('category')

sq.gr.spatial_neighbors(adata,coord_type='generic',radius=50, spatial_key='X_spatial')
sq.gr.nhood_enrichment(adata, cluster_key="subclone2")

subclone_number = adata.obs['subclone2'].value_counts()[adata.obs['subclone2'].cat.categories]
nhood_counts = pd.DataFrame(adata.uns['subclone2_nhood_enrichment']['count'], index=adata.obs['subclone2'].cat.categories, columns=adata.obs['subclone2'].cat.categories)
nhood_percents = nhood_counts/subclone_number

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
#sns.heatmap(df, cmap='YlGn',center = 0,ax=ax,linewidths=0.3,vmin=0,vmax=1)
sns.heatmap(df, cmap=sns.diverging_palette(50, 250, n=100),center = 0,ax=ax,linewidths=0.3,vmin=0,vmax=0.8)
ax.set_facecolor('#dcdcdc')
plt.savefig('%s_region_nhood_heatmap_r50.pdf' % sample, dpi=300, bbox_inches='tight')