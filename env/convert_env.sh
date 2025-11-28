# convert [安装sceasy](https://mp.weixin.qq.com/s/MzMHKjHl-V_XWKYw9DEmCg)
source /opt/software/miniconda3/bin/activate
conda create -n convert r-seurat=4.4 r-devtools python=3.11 -y
conda activate convert
conda install -c conda-forge scanpy -y
conda install bioconda::r-sceasy -y
# Rscript -e 'devtools::install_github("cellgeni/schard")' # fail
wget https://github.com/cellgeni/schard/archive/refs/heads/main.zip
Rscript -e 'devtools::install_local("main.zip")'
rm main.zip
conda install conda-forge::r-reticulate -y
conda install conda-forge::loompy -y
conda install conda-forge::r-optparse -y