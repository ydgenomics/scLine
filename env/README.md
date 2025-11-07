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
||scanpy|
|2_anno|`lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)`|
||`lapply(c("dplyr","Seurat","HGNChelper", "optparse"), library, character.only = T)`|
||`lapply(c("scater","SingleR","SingleCellExperiment"), library, character.only = T)`|

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

```shell
source /software/miniconda/bin/activate
# R
conda create -n scline r-base=4.3 -y
conda activate scline
conda install conda-forge::papermill -y
conda install conda-forge::r-optparse -y
conda install conda-forge::r-irkernel -y
conda install conda-forge::r-rpresto -y
conda install conda-forge::r-biocmanager -y
conda install conda-forge::r-devtools -y
conda install conda-forge::r-remotes
# python
conda install conda-forge::ipykernel -y

# 1_qc
conda install conda-forge::r-seurat -y
conda install bioconda::bioconductor-dropletutils -y
conda install conda-forge::scanpy -y

# 02_anno
conda install bioconda::bioconductor-singler -y
conda install conda-forge::r-hgnchelper -y
```

</details>