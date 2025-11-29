# 测试将qc部分封装为一个conda

```
source /opt/software/miniconda3/bin/activate

# R
conda create -n scline r-base=4.3 -y
conda install bioconda::nextflow -y
conda activate scline
conda install conda-forge::papermill -y
conda install conda-forge::r-optparse -y
conda install conda-forge::r-irkernel -y
conda install conda-forge::r-biocmanager -y
conda install conda-forge::r-devtools -y
conda install conda-forge::r-remotes
# python
conda install conda-forge::ipykernel -y

# 1_qc
conda install conda-forge::r-seurat -y
conda install bioconda::bioconductor-dropletutils -y
conda install conda-forge::scanpy -y
conda install bioconda::scrublet -y
conda install conda-forge::leidenalg -y
conda install conda-forge::r-soupx -y
```