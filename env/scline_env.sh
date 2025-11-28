source /opt/software/miniconda3/bin/activate
conda install bioconda::nextflow -y
# R
conda create -n scline r-base=4.3 -y
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


# 02_anno
conda install bioconda::bioconductor-singler -y
conda install conda-forge::r-hgnchelper -y
conda install conda-forge::r-tidyverse -y
conda install conda-forge::r-ggraph -y
conda install conda-forge::r-data.tree -y
conda install bioconda::bioconductor-scater -y

# 03_integrate
Rscript -e 'install.packages("bbknnR")'
conda install pwwang::r-seuratdata -y
conda install bioconda::r-liger -y
conda install conda-forge::r-rcppplanc -y
conda install pwwang::r-seuratwrappers -y
conda install bioconda::harmonypy -y
conda install conda-forge::scvi-tools -y
conda install conda-forge::scikit-misc -y
conda install conda-forge::r-rcppplanc -y
conda install bioconda::bioconductor-glmgampoi -y
conda install conda-forge::r-harmony -y
# conda install bioconda::scib -y # python版本冲突
# pip install scib-metrics


# 04_metaneighbor
conda install bioconda::bioconductor-metaneighbor -y
conda install bioconda::bioconductor-complexheatmap -y
conda install conda-forge::r-circlize -y
conda install conda-forge::plotly -y

# 05_dea
Rscript -e "install.packages('presto')"

# 06_enrich
conda install bioconda::bioconductor-clusterprofiler -y
conda install bioconda::bioconductor-keggrest -y
conda install bioconda::bioconductor-annotationforge -y
Rscript -e 'install.packages("openxlsx", repos="https://cloud.r-project.org")' # install failed with conda

# # recall
# conda install bioconda::bioconductor-scdesign3 -y
# conda install r::r-lamw -y
# Rscript -e "install.packages('knockoff')"
# Rscript -e 'install.packages("countsplit")'