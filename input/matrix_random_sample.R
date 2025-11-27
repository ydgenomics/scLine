# 随机取10%的数据作为测试数据

library(Matrix)
library(DropletUtils)
library(Seurat)

#' 随机采样并输出 10X V3 格式
#' @param raw_dir    原始 10X 路径（RNA）
#' @param filter_dir 过滤后 10X 路径（可选，NULL 时不处理）
#' @param seeds      随机种子向量，长度 = 要生成的组数
#' @param frac       采样比例（0-1），可向量或单值；向量时与 seeds 一一对应
#' @param out_prefix 输出文件夹前缀，函数会自动加 "_seed{seed}_p{frac*100}"
sample_write_10x <- function(raw_dir, filter_dir = NULL, seeds = 42, frac = 0.1) {
  if (length(frac) == 1)  frac <- rep(frac, length(seeds))
  if (length(frac) != length(seeds))
    stop("seeds 与 frac 长度不一致")
  
  # 读入
  raw   <- Read10X(raw_dir)
  filter <- if (!is.null(filter_dir)) Read10X(filter_dir) else NULL
  
  for (i in seq_along(seeds)) {
    set.seed(seeds[i])
    f <- frac[i]
    
    # 1. 在 raw 里采样
    n_cell   <- ncol(raw)
    keep_idx <- sample(n_cell, max(1, round(n_cell * f)))
    cells_keep <- colnames(raw)[keep_idx]
    raw_sub  <- raw[, cells_keep, drop = FALSE]
    
    # 2. 同步筛选 filter（若提供）
    if (!is.null(filter)) {
      common_cells <- intersect(cells_keep, colnames(filter))
      filter_sub   <- filter[, common_cells, drop = FALSE]
    } else {
      filter_sub   <- NULL
    }
    
    # 3. 文件夹命名：{prefix}_seed{seed}_p{percent}
    pct <- round(f * 100)
    raw_out   <- sprintf("%s_seed%d_p%d", basename(raw_dir), seeds[i], pct)
    filter_out <- sprintf("%s_seed%d_p%d", basename(filter_dir), seeds[i], pct)
    
    # 4. 写 10X V3
    write10xCounts(raw_out, raw_sub, version = "3")
    if (!is.null(filter_sub))
      write10xCounts(filter_out, filter_sub, version = "3")
    
    message("Seed ", seeds[i], " p=", f, " 完成")
  }
}

# 用法示例：3 组不同比例 + 3 个种子
sample_write_10x(
  raw_dir   = "/data/work/scline/input/human/sample1_raw",
  filter_dir = "/data/work/scline/input/human/sample1_filter",
  seeds     = c(11, 12),
  frac      = c(0.1, 0.1)
)

sample_write_10x(
  raw_dir   = "/data/work/scline/input/human/sample1_raw",
  filter_dir = "/data/work/scline/input/human/sample1_filter",
  seeds     = c(21, 22),
  frac      = c(0.1, 0.1)
)

sample_write_10x(
  raw_dir   = "/data/work/scline/input/human/sample1_raw",
  filter_dir = "/data/work/scline/input/human/sample1_filter",
  seeds     = c(31, 32),
  frac      = c(0.1, 0.1)
)
