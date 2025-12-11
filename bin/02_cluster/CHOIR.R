### 251211

library(Seurat)
library(CHOIR)
library(dplyr)
library(ggplot2)  

library(optparse)

option_list <- list(
  make_option(c("-i", "--input_rds"), type = "character", default = "/data/users/yangdong/yangdong_231e5ea9edf9491ca6183fe80520203c/online/scline/output/example1/example1.hr.rds",
              help = "Path to input Seurat RDS file [default: %default]"),
  make_option(c("-n", "--n_cores"), type = "integer", default = 4,
              help = "Number of cores to use for CHOIR functions [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input_rds <- opt$input_rds
n_cores <- opt$n_cores

# 读取Seurat对象
seurat.obj <- readRDS(input_rds)
library(tools)
prefix <- tools::file_path_sans_ext(basename(input_rds))

# 质量控制：排除细胞中reads数少于100的细胞和在少于5个细胞中出现的基因
seurat.obj <- subset(seurat.obj, subset = nCount_RNA > 100)
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 5)

# 数据标准化
seurat.obj <- NormalizeData(seurat.obj)

# 构建层次聚类树
seurat.obj <- buildTree(seurat.obj, n_cores = 2)

# 修剪层次聚类树
seurat.obj <- pruneTree(seurat.obj, n_cores = 2)

# 运行UMAP降维
seurat.obj <- runCHOIRumap(seurat.obj, reduction = "P0_reduction")

# 检查降维结果是否被正确保存
print("Available reductions in the object:")
print(Reductions(seurat.obj))

# 设置当前的聚类结果为 CHOIR_clusters_0.05
Idents(seurat.obj) <- "CHOIR_clusters_0.05"

# 可视化聚类结果
pdf(file = paste0(prefix, "_choir_clusters.pdf"), width = 8, height = 8)
DimPlot(seurat.obj, 
        reduction = "CHOIR_P0_reduction_UMAP",  
        group.by = "CHOIR_clusters_0.05", 
        label = TRUE, 
        pt.size = 0.5)
dev.off()

# 差异表达分析
markers <- FindAllMarkers(seurat.obj, only.pos = TRUE)

# 筛选出平均对数变化大于1的基因
markers <- markers %>% 
    group_by(cluster) %>% 
    filter(avg_log2FC > 1)

# 保存差异表达基因结果
write.csv(markers, file = paste0(prefix, "_choir_markers.csv"), row.names = FALSE)  #对应修改名称

# 提取差异表达基因的特征，每个聚类取平均对数变化最大的5个基因
features <- markers %>% 
    group_by(cluster) %>% 
    top_n(5, avg_log2FC) %>%  
    ungroup() %>% 
    select(gene) %>% 
    unique() %>% 
    pull()

# 绘制点图 - 使用umap降维结果
num <- unique(seurat.obj$CHOIR_clusters_0.05)
pdf(file = paste0(prefix, "_choir_dotplot.pdf"), width = num*2, height = num)
p <- DotPlot(seurat.obj, features = features, dot.scale = 6) + 
    RotatedAxis() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
print(p)
dev.off()

# 用hclust绘制以簇为对象的树状图，并生成所有层级的标记基因列表
# 1. 提取聚类中心
cluster_means <- Seurat::AggregateExpression(
  seurat.obj,
  group.by = "CHOIR_clusters_0.05",
  slot = "data",
  assays = "RNA",
  return.seurat = FALSE
)$RNA

# 2. 计算距离矩阵
dist_matrix <- dist(t(cluster_means))

# 3. 层次聚类
hc <- hclust(dist_matrix, method = "ward.D2")

# 4. 绘制树状图
pdf(paste0(prefix, "_dendrogram.pdf"), width = 8, height = 6)
plot(hc, main = "Cluster Dendrogram", xlab = "CHOIR Clusters")
dev.off()

# 5.从full_tree中提取所有L层级的列名
all_levels <- grep("^L\\d+$", colnames(seurat.obj@misc$CHOIR$clusters$full_tree), value = TRUE)

# length(unique(seurat.obj$L1))  # 检查L1的簇数量

# 6.按数字顺序排序（确保L1 < L2 < L3...）
all_levels <- all_levels[order(as.numeric(gsub("L", "", all_levels)))]  # 这里补全了括号
print(paste("检测到", length(all_levels), "个层级:", paste(all_levels, collapse = ", ")))

# 7.设置输出目录
output_dir <- paste0(prefix, "_CHOIR_markers")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 8.为每个层级生成标记基因列表
for (level in all_levels) {
  # 将当前层级的分群结果添加到Seurat对象
  seurat.obj[[level]] <- seurat.obj@misc$CHOIR$clusters$full_tree[[level]]
  
  # 9.运行差异表达分析
  markers <- FindAllMarkers(
    seurat.obj,
    group.by = level,
    only.pos = TRUE
  )
  
  # 10.保存结果到指定目录
  write.csv(
    markers,
    file = file.path(output_dir, paste0("markers_", level, ".csv")),
    row.names = FALSE
  )
  
  message("已完成层级 ", level, " 的标记基因分析，文件保存至：", 
          file.path(output_dir, paste0("markers_", level, ".csv")))
}

# 保存处理后的RDS
print("Saving processed Seurat object to:")
print(paste0(prefix, "_after_choir.rds"))
saveRDS(seurat.obj, file = paste0(prefix, "_after_choir.rds"))  #对应修改名称