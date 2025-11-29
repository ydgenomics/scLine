# 测试将qc部分封装为一个conda

[制作conda软件包](https://mp.weixin.qq.com/mp/appmsgalbum?__biz=MjM5NTk0Mzg2Nw==&action=getalbum&album_id=4101076046106296324&subscene=7&scenenote=https%3A%2F%2Fmp.weixin.qq.com%2Fs%3Fsearch_click_id%3D2187534805405378034-1764386102813-9123830668%26__biz%3DMjM5NTk0Mzg2Nw%3D%3D%26mid%3D2247489793%26idx%3D1%26sn%3D92fb84cdffedfe3c1273c1bd297789ad%26chksm%3Da7b5585d29ae15808d493d389b12b4695b271e1e28b81e5800b80963ee0f17e5f6fb39a1eb36%26scene%3D7%26key%3Ddaf9bdc5abc4e8d0fe2176adde0406fa22b6dd8a9842aa0f7f4d92e9695455ed0945b80210315f9e9fc841a57e194190a0ca5ca524c3279d05faab9d79ab311120360580816a45bf950090210d1fc9a10d58b52125660188c8133ddbc46deee9092a62e5a508b2173f77c6e76e2e4bfcf8b9ba099d526630a2069b8273c1c67c%26ascene%3D0%26uin%3DNDIxMzk4MTk3%26devicetype%3DUnifiedPCWindows%26version%3Df254140d%26lang%3Dzh_CN%26countrycode%3DCN%26exportkey%3Dn_ChQIAhIQe3X3N3Is6TDtWXUpw%252FDxIRLoAQIE97dBBAEAAAAAAOMoFO0Ggz8AAAAOpnltbLcz9gKNyK89dVj0vBqGF%252FV1S98bTorKeU50uJXC71y%252FZ1rmtqcd%252BvgPRL9Q5MZ2hYi3ar7uqWgrusDyS%252B8%252FtwPVqNWhoPp%252BsFInbLBKD1YBQ6jOtNshlgLfIK5OK5zLhWaqmUhhupyyk3%252BtS4y%252Buw1ClluOkOTYy9aIyNCM%252FJ2hB5ot8Etyi1ePw34MpMviS1w248nEdmdtD32ECt%252BXDlNcHVyoiDMY8dXxh%252Bevl89hTGgblaSQIhpmIj7Sb%252FBiaqUnyufjLvzez%252Fn1uZc%253D%26acctmode%3D0%26pass_ticket%3DhDRA4x%252Fkwg4tBaDQ0XRMTj8aBJBUAksSN%252BmnqKzEk0ID3vkk1xhXgjVbFsm9Cqu9%26wx_header%3D0&nolastread=1&sessionid=&scene=1#wechat_redirect)

## conda build
**配置conda build环境**
```shell
conda create -n conda-build r-base=4.3 python=3.13 -y
conda install conda-forge::conda-build -y
conda install anaconda::anaconda-client -y
anaconda login
```
**提交包到conda**
```shell
# 方法 1：显式指定文件
anaconda upload \
    /home/ydgenomics/miniconda3/envs/conda-build/conda-bld/noarch/scqc-1.0.0-0.conda

# 方法 2：让客户端自动找最新构建
# anaconda upload $(conda-build . --output)
```
安装mamba `conda install conda-forge::mamba -y`

**清除缓存**
```shell
# 清掉所有中间目录与旧 tarball
conda build purge

# 如果还想把已下载的源码包也删掉（通常几百 MB）
conda build purge-all
```

## 测试创建scqc
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