input_h5ad=$1
group_key=$2
control_value=$3
sample_key=$4
cluster_key=$5
split_cluster_py=$6
dea_memento_py=$7
capture_rate=0.07
min_perc=0.7
pval_threshold=0.05
n_cpu=4
top_number=5
perform_2d_test='no' # "yes" means run 2d test(coefficient)


python $split_cluster_py \
--input $input_h5ad \
--group_key $group_key \
--cluster_key $cluster_key

for f in *.h5ad; do
    python $dea_memento_py \
    --input_h5ad $f --group_key $group_key --control_value $control_value --sample_key $sample_key \
    --capture_rate $capture_rate --min_perc $min_perc --pval_threshold $pval_threshold \
    --n_cpu $n_cpu --top_number $top_number --perform_2d_test $perform_2d_test
done