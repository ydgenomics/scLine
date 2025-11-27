# 251125

import sys
import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import cellrank as cr
import scanpy as sc
import scvelo as scv
import argparse

# 1. 参数定义
parser = argparse.ArgumentParser(description="Deal with integrated h5ad")
parser.add_argument('--input_h5ad',
                    type=str,
                    default='/data/work/scline/OUTPUT/integrate/cotton_harmony_integrated_dealplus.h5ad',
                    help='input h5ad path')
parser.add_argument('--batch_key',
                    type=str,
                    default='biosample',
                    help='batch key in adata.obs')
parser.add_argument('--cluster_key',
                    type=str,
                    default='celltype',
                    help='cluster key in adata.obs')
parser.add_argument('--cluster_value',
                    type=str,
                    default='1|2',
                    help='cluster values to keep (pipe-separated)')
args = parser.parse_args()
input_h5ad = args.input_h5ad
batch_key = args.batch_key
cluster_key = args.cluster_key
cluster_value = args.cluster_value

h5ad2pca = {'_scVI':'X_scVI', '_harmony':'X_pca_harmony', '_rliger.INMF':'X_inmfnorm', '_SCTransform.CCA':'X_pca', '_SCTransform.harmony':'X_pca'}
name = os.path.splitext(os.path.basename(input_h5ad))[0]

# 查找第一个匹配的 key
key = next((k for k in h5ad2pca if k in name), None)

# 若命中则取对应值，否则回退到默认 'X_pca'
pca_value = h5ad2pca[key] if key else 'X_pca'
print(pca_value)

adata2 = sc.read_h5ad(input_h5ad)
print(adata2)

cluster_value = cluster_value.split('|')
print(cluster_value)
adata = adata2[adata2.obs[cluster_key].isin(cluster_value)].copy()

adata.X = adata.layers['counts'].copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key = batch_key)
sc.pp.neighbors(adata, use_rep = pca_value)
sc.tl.umap(adata)
# Using the igraph implementation and a fixed number of iterations can be significantly faster,
# especially for larger datasets
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2)

sc.pl.scatter(adata, basis='umap', color=[batch_key, cluster_key, "leiden"], legend_loc='on data', save=f'{name}_umap.pdf')



adata.layers["spliced"] = adata.X
adata.layers["unspliced"] = adata.X

# scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
scv.pp.moments(adata)

ctk = cr.kernels.CytoTRACEKernel(adata)
ctk.compute_cytotrace()

scv.pl.scatter(
    adata,
    c=["ct_pseudotime", 'leiden'],
    basis="umap",
    legend_loc="right",
    color_map="viridis",
    save=f'{name}_ct_pseudotime.pdf'
)

df = adata.obs[["ct_pseudotime", 'leiden']].copy()
sns.set_style(style="whitegrid")
fig, ax = plt.subplots(figsize=(6, 4))
sns.violinplot(
    data=df,
    x='leiden',
    y="ct_pseudotime",
    scale="width",
    #palette=["#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"],
    ax=ax,
)
ax.tick_params(axis="x", rotation=45)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
plt.savefig(f'{name}_violin_pseudotime.pdf', bbox_inches='tight')
# plt.show()
sns.reset_orig()

adata.write_h5ad(f"{name}_ct.h5ad",compression='gzip')

data = adata.obs
data = data[['ct_pseudotime', "leiden"]]

data.index = data.index.astype(str)
# 1) 建立全量 DataFrame，缺失填 NaN
full_df = data.reindex(adata2.obs.index)
# 2) 把列写回 adata2.obs
for col in full_df.columns:
    adata2.obs[col] = full_df[col]

#sc.pl.scatter(adata2, basis="umap", color=["Ga12g00054", "Ga14g00277"], legend_loc='on data')
sc.pl.scatter(adata2, basis="umap", color=[cluster_key, "ct_pseudotime", "leiden"], legend_loc='on data', save=f'{name}_ct_all.pdf')