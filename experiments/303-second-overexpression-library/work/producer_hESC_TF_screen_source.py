import os
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
from tqdm import tqdm

import anndata
import scanpy as sc

from cytofuture_data.gene_name_mapping import GeneNameMapper

# ---- cell ----
# Load the gene name mapper
gene_name_mapper = GeneNameMapper(
    '../../standard_genes/gene_names/human_genes.csv',
    '../../standard_genes/gene_names/mouse_genes.csv',
    '../../standard_genes/gene_names/orthologue_map_human2mouse_best.csv',
    '../../standard_genes/gene_names/orthologue_map_mouse2human_best.csv'
)

# ---- cell ----
tf_screen_df = pd.read_csv(
    '/home/xingjie/Data/data2/cytofuture/datasets/TF_orf_screen/TF_hESC_screen.csv',
    index_col=0)

tf_screen_df.columns = [c.strip() for c in tf_screen_df.columns]

tf_screen_df = 10 ** tf_screen_df
tf_screen_df = 1e4 * tf_screen_df / tf_screen_df.sum(axis=0)

# Remove genes that are not detected in some samples
tf_screen_df = tf_screen_df.loc[tf_screen_df.min(axis=1) > 1e-6].copy()

tf_screen_df.index = gene_name_mapper.map_gene_names(tf_screen_df.index, 
                                'human', 'human', 'name', 'id')
tf_screen_df = tf_screen_df.loc[tf_screen_df.index.str.startswith('ENSG')].copy()

tf_screen_df = tf_screen_df[~tf_screen_df.index.duplicated(keep='first')]
tf_screen_df

# ---- cell ----
for gn in tf_screen_df.columns:
    gid = gene_name_mapper.map_gene_names([gn], 'human', 'human', 'name', 'id')[0]
    if gid == 'na':
        print(gn)

# ---- cell ----
tf_screen_df = tf_screen_df.rename(columns={
    'H2AFZ':'H2AZ1', 'HIST2H3C':'H3C14', 'SSX6':'SSX6P', 'T':'TBXT',
})
tf_screen_df = tf_screen_df.loc[:,~tf_screen_df.columns.duplicated(keep='first')].copy()
tf_screen_df

# ---- cell ----
standard_gene_df = pd.read_csv('standard_genes.csv').set_index('human_id')

obs_dict = {
    'id' : [],
    'condition' : [],
    'perturbed_gene' : [],
    'perturbation_sign' : [],
}
X_perturb = []
X_control = []
X_measure_mask = []
dataset_id = 'hESC_TF_screen'
control_exp = tf_screen_df['Emerald']

# Get the perturbation data
for c in tqdm(tf_screen_df.columns):
    pg = gene_name_mapper.map_gene_names([c], 'human', 'human', 'name', 'id')[0]
    if not pg.startswith('ENSG'):
        continue

    obs_dict['id'].append(f'{dataset_id}_{c}')
    obs_dict['condition'].append(dataset_id)
    obs_dict['perturbed_gene'].append(pg)
    obs_dict['perturbation_sign'].append(1) # -1 for knockdown, 1 for overexpression
    
    standard_gene_df['perturb'] = tf_screen_df[c]
    standard_gene_df['control'] = control_exp
    standard_gene_df['measure_mask'] = ~standard_gene_df['perturb'].isna()

    X_perturb.append(np.nan_to_num(standard_gene_df['perturb'].values, 
                           nan=0).astype(np.float32))
    X_control.append(np.nan_to_num(standard_gene_df['control'].values, 
                           nan=0).astype(np.float32))
    X_measure_mask.append(standard_gene_df['measure_mask'].values.astype(np.float32))

X_perturb = np.array(X_perturb)
X_control = np.array(X_control)
X_perturb = np.log1p(X_perturb / X_perturb.sum(axis=1, keepdims=True) * 1e4)
X_control = np.log1p(X_control / X_control.sum(axis=1, keepdims=True) * 1e4)
X_shift = X_perturb - X_control

X_measure_mask = np.array(X_measure_mask)

# Create the perturbation anndata
adata_perturb = anndata.AnnData(
    X=X_shift, 
    obs=pd.DataFrame(obs_dict).set_index('id'),
    var=standard_gene_df[[]].copy(),
)
adata_perturb.layers['measure_mask'] = X_measure_mask
adata_perturb.layers['control'] = X_control

adata_perturb.obs['perturbed_gene_name'] = gene_name_mapper.map_gene_names(
    adata_perturb.obs['perturbed_gene'], 'human', 'human', 'id', 'name')
adata_perturb.var['gene_name'] = gene_name_mapper.map_gene_names(
    adata_perturb.var.index, 'human', 'human', 'id', 'name')

adata_perturb.write_h5ad(os.path.join(
    '/home/xingjie/Data/data2/cytofuture/datasets/perturbation_data', 
    f'{dataset_id}.h5ad'))

plt.hist(adata_perturb.X.flatten(), bins=100)
adata_perturb

# ---- cell ----

