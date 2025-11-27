nextflow run main.nf --mode qc \
--rawPaths '/data/work/scline/example/SRR31330969_raw|/data/work/scline/example/SRR31330970_raw' \
--filterPaths '/data/work/scline/example/SRR31330969_filter|/data/work/scline/example/SRR31330970_filter' \
--biosampleValues 'group1|group2' \
--sampleValues 'sample1|sample2' \
--runsoupx 'yes' \
--tfidfMin 0 \
--outdir './cotton' 

nextflow run main.nf --mode qc \
--rawPaths '/data/work/SRR31330969_raw_seed11_p10|/data/work/SRR31330969_raw_seed12_p10|/data/work/SRR31330969_raw_seed13_p10' \
--filterPaths '/data/work/SRR31330969_filter_seed11_p10|/data/work/SRR31330969_filter_seed12_p10|/data/work/SRR31330969_filter_seed13_p10' \
--biosampleValues 'group1|group1|group1' \
--sampleValues 'sample1|sample2|sample3' \
--runsoupx 'yes' \
--tfidfMin 0 \
--outdir './cotton' 

nextflow run main.nf --mode all \
--rawPaths '/data/work/scline/input/sample1/sample1_raw|/data/work/scline/input/sample2/sample2_raw' \
--filterPaths '/data/work/scline/input/sample1/sample1_filter|/data/work/scline/input/sample2/sample2_filter' \
--biosampleValues 'group1|group1' \
--sampleValues 'sample1|sample2' \
--python_env '/opt/software/miniconda3/envs/scline/bin/python' \
--rdsORcsv '/data/work/scline/input/markergene.csv' \
--cluster_key 'leiden_res_0.50' \
--prefix 'cotton' \
--other1_key 'biosample' \
--other2_key 'anno' \
--ident_1 'celltype2' --ident_2 'celltype3' \
--cluster_value 'celltype2|celltype3' \
--outdir 'results3' \
--runsoupx 'yes' \
--runsct 'yes' \
--emapper_xlsx '/data/work/scline/input/eggNOG_annotation.xlsx'

# ----------------------------
nextflow run main.nf --mode all \
--rawPaths '/data/work/scline/input/sample1_raw_seed11_p10|/data/work/scline/input/sample1_raw_seed12_p10' \
--filterPaths '/data/work/scline/input/sample1_filter_seed11_p10|/data/work/scline/input/sample1_filter_seed12_p10' \
--biosampleValues 'group1|group1' \
--sampleValues 'sample1|sample2' \
--python_env '/opt/software/miniconda3/envs/scline/bin/python' \
--rdsORcsv '/data/work/scline/input/markergene.csv' \
--cluster_key 'leiden_res_0.50' \
--prefix 'example1' \
--other1_key 'biosample' \
--other2_key 'anno' \
--ident_1 'celltype2' --ident_2 'celltype3' \
--cluster_value 'celltype2' \
--outdir './output/example1' \
--runsoupx 'yes' \
--runsct 'yes' \
--emapper_xlsx '/data/work/scline/input/eggNOG_annotation.xlsx' \
--rlst '0.5'

nextflow run main.nf --mode all \
--rawPaths '/data/work/scline/input/sample1_raw_seed11_p10|/data/work/scline/input/sample1_raw_seed12_p10|/data/work/scline/input/sample1_raw_seed21_p10|/data/work/scline/input/sample1_raw_seed22_p10' \
--filterPaths '/data/work/scline/input/sample1_filter_seed11_p10|/data/work/scline/input/sample1_filter_seed12_p10|/data/work/scline/input/sample1_filter_seed21_p10|/data/work/scline/input/sample1_filter_seed22_p10' \
--biosampleValues 'group1|group1|group2|group2' \
--sampleValues 'sample1|sample2|sample3|sample4' \
--python_env '/opt/software/miniconda3/envs/scline/bin/python' \
--rdsORcsv '/data/work/scline/input/markergene.csv' \
--cluster_key 'leiden_res_0.50' \
--prefix 'example2' \
--other1_key 'biosample' \
--other2_key 'anno' \
--ident_1 'celltype2' --ident_2 'celltype3' \
--cluster_value 'celltype2' \
--outdir './output/example2' \
--runsoupx 'yes' \
--runsct 'yes' \
--emapper_xlsx '/data/work/scline/input/eggNOG_annotation.xlsx' \
--rlst '0.5'


nextflow run main.nf --mode all \
--rawPaths '/data/work/scline/input/sample1_raw_seed11_p10|/data/work/scline/input/sample1_raw_seed12_p10|/data/work/scline/input/sample1_raw_seed21_p10|/data/work/scline/input/sample1_raw_seed22_p10|/data/work/scline/input/sample1_raw_seed31_p10|/data/work/scline/input/sample1_raw_seed32_p10' \
--filterPaths '/data/work/scline/input/sample1_filter_seed11_p10|/data/work/scline/input/sample1_filter_seed12_p10|/data/work/scline/input/sample1_filter_seed21_p10|/data/work/scline/input/sample1_filter_seed22_p10|/data/work/scline/input/sample1_filter_seed31_p10|/data/work/scline/input/sample1_filter_seed32_p10' \
--biosampleValues 'time1|time1|time2|time2|time3|time3' \
--sampleValues 'sample1|sample2|sample3|sample4|sample5|sample6' \
--python_env '/opt/software/miniconda3/envs/scline/bin/python' \
--rdsORcsv '/data/work/scline/input/markergene.csv' \
--cluster_key 'leiden_res_0.50' \
--prefix 'example3' \
--other1_key 'biosample' \
--other2_key 'anno' \
--ident_1 'celltype2' --ident_2 'celltype3' \
--cluster_value 'celltype2' \
--outdir './output/example1' \
--runsoupx 'yes' \
--runsct 'yes' \
--emapper_xlsx '/data/work/scline/input/eggNOG_annotation.xlsx' \
--rlst '0.5'