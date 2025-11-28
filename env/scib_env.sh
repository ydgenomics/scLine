# scib

# scib-metrics
conda create -n scib python=3.11 -y
conda activate scib
conda install -c conda-forge gcc=12 gxx=12 -y
conda install conda-forge::h5py -y
pip install scib-metrics # Python >=3.11
conda install conda-forge::scanpy -y
pip install --use-pep517 bbknn #pip install bbknn
conda install conda-forge::ipykernel -y

# 07_pseudotime
source /opt/software/miniconda3/bin/activate
conda activate scib
conda install -c conda-forge cellrank -y
conda install conda-forge::moscot -y
# conda install -c conda-forge -c bioconda palantir -y #python要求3.12
# pip install --user magic-impute
conda install conda-forge::certifi -y
conda install conda-forge::ipykernel -y
conda install conda-forge::rpy2 -y
# pip install fa2-modified
conda install conda-forge::fa2 -y #不支持安装在大于等于3.12python
conda install anaconda::click -y
# python -m ipykernel install --user --name cellrank2 --display-name "Python (cellrank2)"