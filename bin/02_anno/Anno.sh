input_query_rds=$1
cluster_key=$2
rdsORcsv=$3 # "/data/work/group1.hr_sctype.rds" # "/data/work/scline/bin/02_anno/markergene.csv"
umap_name=$4
ref_cluster_key=$5
anno_singler_r=$6
plot_r=$7
anno_sctype_r=$8

ext="${rdsORcsv##*.}"

pwd=$(pwd)
input_query_rds=("${input_query_rds[@]/#/$pwd/}")
rdsORcsv=("${rdsORcsv[@]/#/$pwd/}")
simple_path=$(basename "$input_query_rds" .rds)

# mkdir -p anno
# cd anno

if [ "$ext" == "rds" ]; then
    Rscript $anno_singler_r --input_ref_rds "$rdsORcsv" --ref_cluster_key "$ref_cluster_key" \
    --input_query_rds "$input_query_rds" --query_cluster_key "$cluster_key" --umap_name "$umap_name"
elif [ "$ext" == "csv" ]; then
    Rscript $plot_r --input_rds "$input_query_rds" --markers_csv "$rdsORcsv" --cluster_key "$cluster_key"
    Rscript $anno_sctype_r --input_marker_csv "$rdsORcsv" --input_query_rds "$input_query_rds" \
    --cluster_key "$cluster_key" --umap_name "$umap_name"
    mv report.txt "$simple_path"_report.txt
else
    echo "Unsupported file extension: $ext. Only 'rds' and 'csv' are supported."
    exit 1
fi