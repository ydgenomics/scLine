## Using conda builds available running environment

- papermill 终端离线运行notebook
- optparse R脚本传参
- irkernel 运行notebook的R核心
- presto R更快运行检验
- biocmanager R包管理
- devtools R包管理
- remotes R包管理
- ipykernel 运行notebook的python核心

|subtask|software|
|-|-|
|1_qc|`lapply(c("Seurat","DropletUtils","SoupX", "optparse"), library, character.only = T)`|
||`import scanpy; import leidenalg; import scrublet`|
|2_anno|`lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)`|
||`lapply(c("dplyr","Seurat","HGNChelper", "optparse"), library, character.only = T)`|
||`lapply(c("scater","SingleR","SingleCellExperiment"), library, character.only = T)`|
|03_integrate|`lapply(c("bbknnR", "patchwork", "SeuratData", "optparse", "magrittr"), library, character.only = T)`|
||`lapply(c("rliger", "harmony", "SeuratWrappers"), library, character.only = T)`|
||`import harmonypy; import scvi`|
||`import numpy as np; import scanpy as sc; import matplotlib.pyplot as plt; import webbrowser; from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection; import click`|
|04_metaneighbor|`lapply(c("MetaNeighbor", "ComplexHeatmap", "circlize"), library, character.only = T)`|
||`import plotly.graph_objects as go`|
|06_enrich|`lapply(c("clusterProfiler", "stringr", "KEGGREST", "AnnotationForge"), library, character.only = T)`|
||`lapply(c("purrr", "RCurl", "data.table", "readxl"), library, character.only = T)`|
|convert|`lapply(c("Seurat", "sceasy", "schard", "reticulate"), library, character.only = T)`|
|07_pseudotime|`import sys; import os; import numpy as np; import pandas as pd; import matplotlib.pyplot as plt; import seaborn as sns; import cellrank as cr; import scanpy as sc; import scvelo as scv`|

<details> <summary> 查看机器 </summary>

```shell
## stereonote-workflow-basic
root@e20515d48d96:/# cat /etc/os-release
PRETTY_NAME="Ubuntu 22.04.5 LTS"
NAME="Ubuntu"
VERSION_ID="22.04"
VERSION="22.04.5 LTS (Jammy Jellyfish)"
VERSION_CODENAME=jammy
ID=ubuntu
ID_LIKE=debian
HOME_URL="https://www.ubuntu.com/"
SUPPORT_URL="https://help.ubuntu.com/"
BUG_REPORT_URL="https://bugs.launchpad.net/ubuntu/"
PRIVACY_POLICY_URL="https://www.ubuntu.com/legal/terms-and-policies/privacy-policy"
UBUNTU_CODENAME=jammy
```

</details>

<details> <summary> 安装wgt并配置miniconda </summary>

```shell
apt update
apt install wget
wget --version
mkdir software; cd software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # install miniconda
bash Miniconda3-latest-Linux-x86_64.sh # will installed on ~/miniconda3
source /software/miniconda/bin/activate
# conda config --set auto_activate_base false
```

</details>

<details> <summary> 配置scline环境安装依赖包 </summary>

# scline
```shell
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

# recall
conda install bioconda::bioconductor-scdesign3 -y
conda install r::r-lamw -y
Rscript -e "install.packages('knockoff')"
Rscript -e 'install.packages("countsplit")'
```

# scib
```shell
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
```

# convert [安装sceasy](https://mp.weixin.qq.com/s/MzMHKjHl-V_XWKYw9DEmCg)
```shell
source /opt/software/miniconda3/bin/activate
conda create -n convert r-seurat=4.4 r-devtools python=3.11 -y
conda activate convert
conda install -c conda-forge scanpy -y
conda install bioconda::r-sceasy -y
# Rscript -e 'devtools::install_github("cellgeni/schard")'
wget https://github.com/cellgeni/schard/archive/refs/heads/main.zip
Rscript -e 'devtools::install_local("main.zip")'
rm main.zip
conda install conda-forge::r-reticulate -y
conda install conda-forge::loompy -y
conda install conda-forge::r-optparse -y
```

</details>