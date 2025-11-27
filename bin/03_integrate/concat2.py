### Date: 251126
### Image: harmony-py--
### Coder: ydgenomics

import pandas as pd
import scanpy as sc
import anndata as ad
import sys
import os

def sort_files_by_projects(files_txt_path, anno_csvs, projects_txt_path, species, group_key):
    indataget = files_txt_path.split(' ')
    anno_csvs_list = anno_csvs.split(' ')
    projects = projects_txt_path.split(' ')
    
    print("原始indataget:", indataget)
    print("原始anno_csvs:", anno_csvs_list)
    print("projects顺序:", projects)
    
    # 按照projects顺序重新排列indataget
    sorted_indataget = []
    sorted_anno_csvs = []
    
    for project in projects:
        # 查找对应的indataget文件
        indata_file = next((f for f in indataget if project in f), None)
        if indata_file:
            sorted_indataget.append(indata_file)
        else:
            print(f"警告: 未找到项目 {project} 对应的indataget文件")
        
        # 查找对应的anno_csvs文件
        anno_file = next((f for f in anno_csvs_list if project in f), None)
        if anno_file:
            sorted_anno_csvs.append(anno_file)
        else:
            print(f"警告: 未找到项目 {project} 对应的anno_csvs文件")
    
    print("排序后的indataget:", sorted_indataget)
    print("排序后的anno_csvs:", sorted_anno_csvs)
    
    return sorted_indataget, sorted_anno_csvs, projects, species, group_key

# 使用示例
if __name__ == "__main__":
    if len(sys.argv) >= 6:
        files_txt_path = sys.argv[1]
        anno_csvs = sys.argv[2]
        projects_txt_path = sys.argv[3]
        species = sys.argv[4]
        group_key = sys.argv[5]
        
        indataget, anno_csvs, projects, species, group_key = sort_files_by_projects(
            files_txt_path, anno_csvs, projects_txt_path, species, group_key
        )
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
    else:
        print("参数不足")