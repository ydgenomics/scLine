# **scLine**: **s**ingle-**c**ell analysis pipe**line**

---

- **目标**
  - 解决两种类型的单细胞数据的端到端的分析需求
    - 两组数据(对照实验)
    - 时序数据
  - 基于conda环境配置，使用已有R和python软件的API，通过git进行部署
  - 操作使用sh参数命令
- **流程设计**
- **测试**
  - 拟南芥对照数据
  - 拟南芥时序数据
  - 参考资料
    - [【Stress Biology】单细胞RNA测序揭示了拟南芥愈伤组织形成的发育轨迹和环境调控](https://mp.weixin.qq.com/s/ZKYanCM-ZalgHteLOaJErw) 
    - [Single-cell RNA sequencing reveals developmental trajectories and environmental regulation of callus formation in Arabidopsis](https://link.springer.com/article/10.1007/s44154-025-00255-4)
- **参考资料**
  |Software|Year|Article|
  |-|-|-|
  |**scDown**|2025|Sun, L., Ma, Q., Cai, C., Labaf, M., Jain, A., Dias, C., Rockowitz, S., & Sliz, P. (2025). scDown: A Pipeline for Single-Cell RNA-Seq Downstream Analysis. International journal of molecular sciences, 26(11), 5297. https://doi.org/10.3390/ijms26115297|
  |SCMeTA|2024|Pan, X., Pan, S., Du, M., Yang, J., Yao, H., Zhang, X., & Zhang, S. (2024). SCMeTA: a pipeline for single-cell metabolic analysis data processing. Bioinformatics (Oxford, England), 40(9), btae545. https://doi.org/10.1093/bioinformatics/btae545|
  |Panpipes|2024|Curion, F., Rich-Griffin, C., Agarwal, D., Ouologuem, S., Rue-Albrecht, K., May, L., Garcia, G. E. L., Heumos, L., Thomas, T., Lason, W., Sims, D., Theis, F. J., & Dendrou, C. A. (2024). Panpipes: a pipeline for multiomic single-cell and spatial transcriptomic data analysis. Genome biology, 25(1), 181. https://doi.org/10.1186/s13059-024-03322-7|
  |MultiSC|2024|Lin, X., Jiang, S., Gao, L., Wei, Z., & Wang, J. (2024). MultiSC: a deep learning pipeline for analyzing multiomics single-cell data. Briefings in bioinformatics, 25(6), bbae492. https://doi.org/10.1093/bib/bbae492|
  ||2023|Heumos, L., Schaar, A. C., Lance, C., Litinetskaya, A., Drost, F., Zappia, L., Lücken, M. D., Strobl, D. C., Henao, J., Curion, F., Single-cell Best Practices Consortium, Schiller, H. B., & Theis, F. J. (2023). Best practices for single-cell analysis across modalities. Nature reviews. Genetics, 24(8), 550–572. https://doi.org/10.1038/s41576-023-00586-w|
  |**single-cell-best-practice**|2023-|https://www.sc-best-practices.org/|
  |IBRAP|2023|Knight, C. H., Khan, F., Patel, A., Gill, U. S., Okosun, J., & Wang, J. (2023). IBRAP: integrated benchmarking single-cell RNA-sequencing analytical pipeline. Briefings in bioinformatics, 24(2), bbad061. https://doi.org/10.1093/bib/bbad061|
  |**SCP**|2023|https://github.com/zhanghao-njmu/SCP|
  ||2022|Bertolini, A., Prummer, M., Tuncel, M. A., Menzel, U., Rosano-González, M. L., Kuipers, J., Stekhoven, D. J., Tumor Profiler consortium, Beerenwinkel, N., & Singer, F. (2022). scAmpi-A versatile pipeline for single-cell RNA-seq analysis from basics to clinics. PLoS computational biology, 18(6), e1010097. https://doi.org/10.1371/journal.pcbi.1010097|
  ||2021|Shaw, R., Tian, X., & Xu, J. (2021). Single-Cell Transcriptome Analysis in Plants: Advances and Challenges. Molecular plant, 14(1), 115–126. https://doi.org/10.1016/j.molp.2020.10.012|
  ||2020|Ding, J., Adiconis, X., Simmons, S. K., Kowalczyk, M. S., Hession, C. C., Marjanovic, N. D., Hughes, T. K., Wadsworth, M. H., Burks, T., Nguyen, L. T., Kwon, J. Y. H., Barak, B., Ge, W., Kedaigle, A. J., Carroll, S., Li, S., Hacohen, N., Rozenblatt-Rosen, O., Shalek, A. K., Villani, A. C., … Levin, J. Z. (2020). Systematic comparison of single-cell and single-nucleus RNA-sequencing methods. Nature biotechnology, 38(6), 737–746. https://doi.org/10.1038/s41587-020-0465-8|
  |**scTyper**|2020|Choi, J. H., In Kim, H., & Woo, H. G. (2020). scTyper: a comprehensive pipeline for the cell typing analysis of single-cell RNA-seq data. BMC bioinformatics, 21(1), 342. https://doi.org/10.1186/s12859-020-03700-5.|
  |**SINCERA**|2015|Guo, M., Wang, H., Potter, S. S., Whitsett, J. A., & Xu, Y. (2015). SINCERA: A Pipeline for Single-Cell RNA-Seq Profiling Analysis. PLoS computational biology, 11(11), e1004575. https://doi.org/10.1371/journal.pcbi.1004575|

---

<details> <summary> <strong> 手册 </strong> </summary>

## 1.引言

<details> <summary> <strong> 引言 </strong> </summary>

### 1.1 编写目的

### 1.2 项目背景

### 1.3 定义

### 1.4 参考资料

</details>

## 2. 软件概述

<details> <summary> <strong> 软件概述 </strong> </summary>

### 2.1 目标

### 2.2 功能

### 2.3 性能

#### 2.3.1 精度

#### 2.3.2 时间特性

#### 2.3.3 灵活性

</details>

## 3. 运行环境

<details> <summary> <strong> 运行环境 </strong> </summary>

### 3.1 硬件

### 3.2 支持软件

</details>


## 4. 使用说明

<details> <summary> <strong> 软件概述 </strong> </summary>

### 4.1 安装和初始化

### 4.2 输入

### 4.3 输出

### 4.4 出错和恢复

</details>

## 5. 运行说明

## 6. 非常规过程

## 7. 操作命令一览表

## 8. 程序文件和数据文件一览表

</details>