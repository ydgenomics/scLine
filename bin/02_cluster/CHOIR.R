# 加载包
library(Seurat)
library(CHOIR)
library(dplyr)
library(ggplot2)  

# 读取Seurat对象
seurat.obj <- readRDS("/data/work/path/")

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
pdf(file = "/data/work/new_D10_merged_choir_clusters.pdf", width = 8, height = 8)
DimPlot(seurat.obj, 
        reduction = "CHOIR_P0_reduction_UMAP",  
        group.by = "CHOIR_clusters_0.05", 
        label = TRUE, 
        pt.size = 0.5)
dev.off()

# 保存处理后的RDS
print("Saving processed Seurat object to:")
print("/data/work/new_D10_merged_obj_after_choir.Rds")
saveRDS(seurat.obj, file = "/data/work/name_obj_after_choir.Rds")  #对应修改名称

# 差异表达分析
markers <- FindAllMarkers(seurat.obj, only.pos = TRUE)

# 筛选出平均对数变化大于1的基因
markers <- markers %>% 
    group_by(cluster) %>% 
    filter(avg_log2FC > 1)

# 保存差异表达基因结果
print("Saving marker genes to:")
print("/data/work/7112.2_markers.csv")
write.csv(markers, file = "/data/work/markers.csv", row.names = FALSE)  #对应修改名称

# 提取差异表达基因的特征，每个聚类取平均对数变化最大的5个基因
features <- markers %>% 
    group_by(cluster) %>% 
    top_n(5, avg_log2FC) %>%  
    ungroup() %>% 
    select(gene) %>% 
    unique() %>% 
    pull()

# 绘制点图 - 使用umap降维结果
print("Saving dot plot to:")
print("/data/work/7112.2_dotplot_adjusted.pdf")
pdf(file = "/data/work/7112.2_dotplot_adjusted.pdf", width = 35, height = 20)  
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
pdf("/data/work/7112.2_cluster_dendrogram.pdf", width = 8, height = 6)
plot(hc, main = "Cluster Dendrogram", xlab = "CHOIR Clusters")
dev.off()

# 5.从full_tree中提取所有L层级的列名
all_levels <- grep("^L\\d+$", colnames(seurat.obj@misc$CHOIR$clusters$full_tree), value = TRUE)

# length(unique(seurat.obj$L1))  # 检查L1的簇数量

# 6.按数字顺序排序（确保L1 < L2 < L3...）
all_levels <- all_levels[order(as.numeric(gsub("L", "", all_levels)))]  # 这里补全了括号
print(paste("检测到", length(all_levels), "个层级:", paste(all_levels, collapse = ", ")))

# 7.设置输出目录
output_dir <- "/data/work/7113_CHOIR_markers"
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