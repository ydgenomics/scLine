### Date: 251121
### Image: harmony-py--
### Coder: ydgenomics

import pandas as pd
import scanpy as sc
import anndata as ad
import sys
import os

files_txt_path = sys.argv[1]
anno_csvs = sys.argv[2]
projects_txt_path = sys.argv[3]
species=sys.argv[4]
group_key=sys.argv[5] #"biosample"

indataget = files_txt_path.split(' ')
print(indataget)
anno_csvs = anno_csvs.split(' ')
print(anno_csvs)
projects = projects_txt_path.split(' ')
print(projects)

adatas={}
for i in range(len(indataget)): 
    key = projects[i]
    value = indataget[i]
    value = sc.read_h5ad(value)
    if 'biosample' in value.obs.columns:
        print('The raw key included `biosample` column, value of raw biosample named to biosample0')
        value.obs['biosample0']=value.obs['biosample']
    if 'counts' in value.layers:
        print("counts exist in raw h5ad")
    else:
        print("counts not exist in raw h5ad, .X as counts")
        value.layers["counts"] = value.X.copy()
    value.X = value.layers["counts"] #ensure concat used by raw data
    csv_path = anno_csvs[i]
    if os.path.isfile(csv_path) and csv_path.lower().endswith('.csv'):
        df = pd.read_csv(csv_path)
        df_idx = df.set_index('barcode')
        value.obs['anno'] = df_idx.loc[value.obs_names, 'anno']
        print(value.obs['anno'].unique())
    adatas[key] = value

adata = ad.concat(adatas, label=group_key, join="inner") # 'inner' or 'outer'
if 'celltype' in adata.obs.columns:
    print('The raw key included `celltype` column, value of raw celltype named to celltype0')
    adata.obs['celltype0']=adata.obs['celltype']

adata.obs_names_make_unique()
print(adata.obs[group_key].value_counts())
print(adata.obs.columns)
adata.write_h5ad(filename=species+'.h5ad',compression="gzip")