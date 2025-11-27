python /data/work/scline/bin/03_integrate/concat.py \
"/data/work/scline/output/qc/result/group1/group1.h5ad|/data/work/scline/output/qc/result/group2/group2.h5ad" \
"group1|group2" cotton biosample

input_file="/data/work/cotton.h5ad"
prefix=$(basename "${input_file%.*}")   # → cotton

# ----- scVI -----
method="scVI"
output_h5ad="${prefix}_${method}_integrated.h5ad"
output_pdf="${prefix}_${method}_integrated_UMAP.pdf"
python /data/work/scline/bin/03_integrate/scVI_integration.py \
$input_file $output_h5ad $output_pdf --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set 0.5

# ----- harmony -----
method="harmony"
output_h5ad="${prefix}_${method}_integrated.h5ad"
output_pdf="${prefix}_${method}_integrated_UMAP.pdf"
python /data/work/scline/bin/03_integrate/harmony_integration.py \
$input_file $output_h5ad $output_pdf --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set 0.5

# ----- unintegration -----
method="harmony"
output_h5ad="${prefix}_${method}_integrated.h5ad"
output_pdf="${prefix}_${method}_integrated_UMAP.pdf"
python /data/work/scline/bin/03_integrate/unintegration.py \
$input_file $output_h5ad $output_pdf --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set 0.5

whether_sct="~{whether_sct}"
if [[ "$whether_sct" == "yes" ]]; then
    echo "Running sct workflow ..."
    # sct.cca
    method="SCTransform.CCA"
    output_rds="${prefix}_${method}_integrated.rds"
    output_pdf="${prefix}_${method}_integrated_UMAP.pdf"
    Rscript /data/work/scline/bin/03_integrate/SCTransform.CCA_integration.R \
    --input_rds $input_file --out_rds $output_rds --out_UMAP $output_pdf \
    --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set 0.5

    # sct.harmony
    method="SCTransform.harmony"
    output_rds="${prefix}_${method}_integrated.rds"
    output_pdf="${prefix}_${method}_integrated_UMAP.pdf"
    /opt/conda/bin/Rscript /WDL/Integration/v1.0.3/SCTransform.harmony_integration.R \
    --input_rds $input_file --out_rds $output_rds --out_UMAP $output_pdf \
    --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set 0.5
fi

# rliger.inmf
method="rliger.INMF"
output_rds="${prefix}_${method}_integrated.rds"
output_pdf="${prefix}_${method}_integrated_UMAP.pdf"
/opt/conda/bin/Rscript /WDL/Integration/v1.0.3/rliger.INMF_integration.R \
--input_rds $input_file --out_rds $output_rds --out_UMAP $output_pdf \
--batch_key biosample --sample_key sample --cluster_key celltype --resolution_set 0.5


/software/conda/Anaconda/bin/python /WDL/scIB/v1.0.1/scIB.py \
$unintegrated_h5ad $integrated_file $methods_file $pcas_file $deals_file $tests_file \
--batch_key $batch_key --label_key $label_key --n_jobs $n_jobs $prefix