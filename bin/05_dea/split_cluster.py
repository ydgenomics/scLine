import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description="Split AnnData by cluster and save per-cluster files.")
parser.add_argument("-i", "--input", default="data/processed/05_dea/adata_clustered.h5ad",
                    help="h5ad")
parser.add_argument("-b", "--group_key", default="batch",
                    help="batch key")
parser.add_argument("-c", "--cluster_key", default="leiden",
                    help="cluster key")
args = parser.parse_args()

input = args.input
group_key = args.group_key
cluster_key = args.cluster_key

adata = sc.read_h5ad(input)
clusters = adata.obs[cluster_key].unique().tolist()

for cluster in clusters:
    adata_cluster = adata[adata.obs[cluster_key] == cluster]
    if len(adata_cluster.obs[group_key].unique()) > 1:
        adata_cluster.write_h5ad(f"adata_{cluster}.h5ad")