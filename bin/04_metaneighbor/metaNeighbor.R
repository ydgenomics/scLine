# Date: 251120
# Image: metaNeighbor-R--04 /opt/conda/bin/R
# Reference: https://mp.weixin.qq.com/s/tVxalBWsxLn58RJkpb-PaQ
# 基于RNA@counts做分析

library(MetaNeighbor)
library(SummarizedExperiment)
library(Seurat)
library(SingleCellExperiment)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_file"),
    type = "character", default = "/data/input/Files/songyouran/resingle-8/dataget/W202509040023043/01_dataget/lettuce.leaf/lettuce.leaf.rds",
    help = "Path to input file"
  ),
  make_option(c("-o", "--output_name"),
    type = "character", default = "cotton",
    help = "Output file prefix name"
  ),
  make_option(c("-b", "--batch_key"),
    type = "character", default = "biosample",
    help = "Batch key for integration"
  ),
  make_option(c("-c", "--cluster_key"),
    type = "character", default = "leiden_res_0.50",
    help = "Cluster key for integration"
  ),
  make_option(c("-t", "--threshold_value"),
    type = "numeric", default = 0.95,
    help = "Threshold value for top hits"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
input_file <- opt$input_file
output_name <- opt$output_name
batch_key <- opt$batch_key
cluster_key <- opt$cluster_key
threshold_value <- opt$threshold_value


sdata <- readRDS(input_file)
sdata
colnames(sdata@meta.data)
sdata <- as.SingleCellExperiment(sdata, assay = "RNA", slot = "counts")
# print(assay(sdata, "counts")[1:5, 1:5])
head(colData(sdata))

var_genes = variableGenes(dat = sdata, exp_labels = sdata@colData[[batch_key]])

celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = sdata,
                             study_id = sdata@colData[[batch_key]],
                             cell_type = sdata@colData[[cluster_key]],
                             fast_version = TRUE)

cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

# in work directory output a pdf
pdf(paste0(output_name,"_metaNeighbor.pdf"), width=5+0.1*length(rownames(celltype_NV)), height=5+0.1*length(rownames(celltype_NV)))
# Using heatmap do heatmap
gplots::heatmap.2(celltype_NV,
                  col = cols,
                  breaks = breaks,
                  key.xlab = "AUROC",
                  margins = c(8, 8),
                  trace = "none",
                  density.info = "none",
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.7,
                  cexCol = 0.7)
dev.off()
write.csv(file=paste0(output_name,"_metaNeighbor.csv"),celltype_NV,quote=FALSE,row.names=TRUE)
top_hits = topHits(cell_NV = celltype_NV,
                   dat = sdata,
                   study_id = sdata@colData[[batch_key]],
                   cell_type = sdata@colData[[cluster_key]],
                   threshold = threshold_value)

write.csv(file=paste0(output_name,"_metaNeighbor_tophits.csv"),top_hits,quote=FALSE,row.names=FALSE)