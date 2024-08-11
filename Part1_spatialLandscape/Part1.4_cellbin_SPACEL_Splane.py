# %%
import sys
import os
from SPACEL.setting import set_environ_seed

set_environ_seed(42)
from SPACEL import Splane

import scanpy as sc
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.serif'] = ['Arial']
sc.settings.set_figure_params(dpi=100,
                              dpi_save=300,
                              facecolor='white',
                              fontsize=10,
                              vector_friendly=True,
                              figsize=(3, 3))

# %%
args = sys.argv
n_cluster = int(args[1])
n_neighbors = int(args[2])
k = int(args[3])
gnn_dropout = float(args[4])

od = 'SPACEL/Splane_cellbin/n_cluster_{}_n_neighbors_{}_k_{}_gnn_dropout_{}/'.format(
    n_cluster, n_neighbors, k, gnn_dropout)
os.system(f'mkdir -p {od}')
os.chdir(od)
save_path = f'{od}/Splane_models/'

# %%
import yaml

# %%
samplelist = [
    "P13", "P15", "P16",
    "P22", "P28", "P29",
    "P35", "P36", "P41"]

celltypelist = ['Atypical hyperplasia epithelial cells','B cells','DCIS cells','Endothelial cells','Fibroblasts','IBC cells','MEC','Macrophages','Mast cells','Normal epithelial cells','Pericytes','Plasmocytes','T cells']
st_ad_list = []
for sample in samplelist:
    adata = sc.read_h5ad(
        f'SPACEL/Splane_cellbin/{sample}_raw_annotated.h5ad'
    )
    adata.uns['celltypes'] = celltypelist
    for celltype in celltypelist:
        if celltype not in adata.obs.columns:
            adata.obs[celltype] = 0
    st_ad_list.append(adata)


# ## Training Splane model
# In this step, we initialize Splane model with a anndata object list as input. The parameter `n_clusters` controls the number of spatial domains.
splane = Splane.init_model(st_ad_list,
                           n_clusters=n_cluster,
                           use_gpu=False,
                           n_neighbors=n_neighbors,
                           k=k,
                           gnn_dropout=gnn_dropout)

# %% [markdown]
# Here, we train the model to obtain latent feature of each spots/cells. The parameter `d_l` affects the level of batch effect correction between slices. By default, `d_l` is `0.2`.

# %%
splane.train(d_l=0.2, save_path=save_path)

# %% [markdown]
# Then, we can identify the spatial domain to which each spot/cell belongs. By default, the results will be saved in `spatial_domain` column in `.obs`. If the `key` parameter is provided, the results will be saved in `.obs[key]`.

# %%
splane.identify_spatial_domain()
merged_adata = sc.AnnData.concatenate(*st_ad_list,
                                      join='inner',
                                      batch_categories=samplelist,
                                      batch_key='batch')
merged_adata.write("cellbin_merged_adata.h5ad")
# %% [markdown]


# ## Plot spatial domain results(Fig1g)
import matplotlib.pyplot as plt
sys.path.insert(
    0,
    'SC_reanalysis/script'
)
from utils import getDefaultColors, removeBiasGenes
# 11 sample with 3row 4col figure
cluster_number = merged_adata.obs['spatial_domain'].cat.categories.shape[0]
colors = getDefaultColors(cluster_number, type=2)
clusters_colors = dict(zip(range(cluster_number), colors))
fig, axs = plt.subplots(3, 4, figsize=(14, 9))
fig.tight_layout() 
plt.axis('off')
for i, sample in enumerate(samplelist):
    nrow = i // 4
    ncol = i % 4
    ad = merged_adata[merged_adata.obs['batch'] == sample]
    sc.pl.spatial(ad,
                  color='spatial_domain',
                  ax=axs[nrow, ncol],
                  show=False,
                  spot_size=30,
                  title=sample,
                  palette=[
                      v for k, v in clusters_colors.items()
                      if k in ad.obs.spatial_domain.unique().tolist()
                  ])
fig.savefig('spatial_domain.pdf', bbox_inches='tight')


## calculate cell type fraction in each cluster(Fig1h)
adata = merged_adata
ct_list = adata.obs['annotated_cluster'].cat.categories
tissue_leiden_list = adata.obs['spatial_domain'].cat.categories
deconv_mat = defaultdict(dict)
for x in tissue_leiden_list:
    for y in ct_list:
        deconv_mat[x][y] = 0
for row in adata.obs.iterrows():
    deconv_mat[row[1]['spatial_domain']][row[1]['annotated_cluster']] += 1
deconv_mat = pd.DataFrame(deconv_mat)
deconv_mat/=deconv_mat.sum()
order = [1,7,3,4,8,0,2,6,9,5]
deconv_mat = deconv_mat[order]
print(deconv_mat)


plt.figure(dpi=300,figsize=(12,7.5))
label = ['CISR','CISAR','IBCR','IBCAR1','IBCAR2','MESR','ICR','MIMR','NMR','UNDR']
Atypical_hyperplasia_epithelial_cell = deconv_mat.iloc[0,0:].values
B_cells = deconv_mat.iloc[1,0:].values
DCIS_cells = deconv_mat.iloc[2,0:].values
Endothelial_cells = deconv_mat.iloc[3,0:].values
Fibroblasts = deconv_mat.iloc[4,0:].values
IBC_cells = deconv_mat.iloc[5,0:].values
MEC = deconv_mat.iloc[6,0:].values
Macrophages = deconv_mat.iloc[7,0:].values
Mast_cells = deconv_mat.iloc[8,0:].values
Normal_epithelial_cells = deconv_mat.iloc[9,0:].values
Pericytes = deconv_mat.iloc[10,0:].values
Plasmocytes = deconv_mat.iloc[11,0:].values
T_cells = deconv_mat.iloc[12,0:].values

width = 0.8
plt.bar(label, T_cells, width, label='T cells',color='#d04dc7',edgecolor='k',lw=1.4)
plt.bar(label, Plasmocytes, width,  bottom=T_cells, label='Plasmocytes',color='#00CED1',edgecolor='k',lw=1.4)
plt.bar(label, Pericytes, width,  bottom=(Plasmocytes+T_cells), label='Pericytes',color='#ce9d3f',edgecolor='k',lw=1.4)
plt.bar(label, Normal_epithelial_cells, width,  bottom=(Plasmocytes+T_cells+Pericytes), label='Normal epithelial cells',color='#607d3a',edgecolor='k',lw=1.4)
plt.bar(label, Mast_cells, width, bottom=(Plasmocytes+T_cells+Pericytes+Normal_epithelial_cells), label='Mast cells',color='#6d76cb',edgecolor='k',lw=1.4)
plt.bar(label, Macrophages, width,  bottom=(Plasmocytes+T_cells+Pericytes+Normal_epithelial_cells+Mast_cells), label='Macrophages',color='#b96f49',edgecolor='k',lw=1.4)
plt.bar(label, MEC, width,  bottom=(Plasmocytes+T_cells+Pericytes+Normal_epithelial_cells+Mast_cells+Macrophages), label='MEC',color='#6bd155',edgecolor='k',lw=1.4)
plt.bar(label, IBC_cells, width,  bottom=(Plasmocytes+T_cells+Pericytes+Normal_epithelial_cells+Mast_cells+Macrophages+MEC), label='IBC cells',color='#FFD700',edgecolor='k',lw=1.4)
plt.bar(label, Fibroblasts, width, bottom=(Plasmocytes+T_cells+Pericytes+Normal_epithelial_cells+Mast_cells+Macrophages+MEC+IBC_cells), label='Fibroblasts',color='#d3c3a4',edgecolor='k',lw=1.4)
plt.bar(label, Endothelial_cells, width,  bottom=(Plasmocytes+T_cells+Pericytes+Normal_epithelial_cells+Mast_cells+Macrophages+MEC+IBC_cells+Fibroblasts), label='Endothelial cells',color='#683ec2',edgecolor='k',lw=1.4)
plt.bar(label, DCIS_cells, width,  bottom=(Plasmocytes+T_cells+Pericytes+Normal_epithelial_cells+Mast_cells+Macrophages+MEC+IBC_cells+Fibroblasts+Endothelial_cells), label='DCIS cells',color='#d5492f',edgecolor='k',lw=1.4)
plt.bar(label, B_cells, width,  bottom=(Plasmocytes+T_cells+Pericytes+Normal_epithelial_cells+Mast_cells+Macrophages+MEC+IBC_cells+Fibroblasts+Endothelial_cells+DCIS_cells), label='B cells',color='#81d8ae',edgecolor='k',lw=1.4)
plt.bar(label, Atypical_hyperplasia_epithelial_cell, width, bottom=(Plasmocytes+T_cells+Pericytes+Normal_epithelial_cells+Mast_cells+Macrophages+MEC+IBC_cells+Fibroblasts+Endothelial_cells+DCIS_cells+B_cells), label='Atypical hyperplasia epithelial cell',color='#81357a',edgecolor='k',lw=1.4)

#plt.tick_params(axis='x',length=0)
plt.xlabel('Region',fontsize=15)
plt.ylabel('Proportion of cell type',fontsize=15)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.ylim(0,1)
handles,labels = plt.gca().get_legend_handles_labels()
order = [12,11,10,9,8,7,6,5,4,3,2,1,0]
plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],frameon=False,bbox_to_anchor=(1.01,1),fontsize=12)
plt.tight_layout()
plt.savefig('bar.png',
            dpi=200,bbox_inches='tight')


