rds=$1
ident_1=$2
ident_2=$3
cluster_key=$4
findmarkers_r=$5
single_volcano_r=$6
only_pos="no" # "yes" or "no"

Rscript $findmarkers_r \
--rds $rds --assay RNA --ident_1 "$ident_1" --ident_2 "$ident_2" --cluster_key $cluster_key --only_pos $only_pos

csvs=$(find . -maxdepth 1 -type f -name "*.csv")
echo $csvs

for input_csv in $(echo $csvs | tr ' ' ' '); do
    echo "$input_csv"
    Rscript $single_volcano_r \
    --gene_csv $input_csv --coef_col "avg_log2FC" --pval_col "p_val_adj" \
    --coef_threshold 1 --pval_threshold 0.01 --n_top 15
done