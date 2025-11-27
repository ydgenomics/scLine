# Date: 250819
# Image: plantphone-R-02
# Findmarkers使用于两分组找差异基因 `assay` `ident.1` `ident.2`
library(Seurat)
library(dplyr)
library(optparse)

option_list <- list(
    make_option(c("-r", "--rds"), type = "character", default = "/data/users/yangdong/yangdong_faff775391984da0a355d4bd70217714/online/cotton/output/merge/cotton_merged.rds", help = "Path to RDS file"),
    make_option(c("-a", "--assay"), type = "character", default = "RNA", help = "Assay to use"),
    make_option(c("-b", "--ident_1"), type = "character", default = "K2", help = "Batch variable"),
    make_option(c("-c", "--ident_2"), type = "character", default = "C1", help = "Batch variable"),
    make_option(c("-g", "--cluster_key"), type = "character", default = "sample", help = "Group variable"),
    make_option(c("-o", "--only_pos"), type = "character", default = "yes", help = "Whether is only focusing positive genes")
    #make_option(c("-p", "--p_threshold"), type = "numeric", default = 0.01, help = "P value of maxiusm")
)
opt <- parse_args(OptionParser(option_list = option_list))
rds <- opt$rds
assay <- opt$assay
ident_1 <- opt$ident_1
ident_2 <- opt$ident_2
cluster_key <- opt$cluster_key
only_pos <- opt$only_pos
#p_threshold <- opt$p_threshold

seu <- readRDS(rds); DefaultAssay(seu) <- assay
seu
only_pos <- only_pos == "yes"

# if (assay == "RNA") {
#     seu <- NormalizeData(seu)
# }

name <- paste0(basename(rds), "_", ident_1, "_vs_", ident_2)

# FindMarkers
# "p_val" "avg_log2FC" "pct.1" "pct.2" "p_val_adj" "cluster" "gene"
Idents(seu) <- seu@meta.data[[cluster_key]]
markers <- FindMarkers(seu,ident.1 = ident_1,ident.2 = ident_2, assay = assay, only.pos = only_pos)
markers$cluster <- paste0(ident_1, "_vs_", ident_2)
markers$gene <- rownames(markers)
#markers <- markers %>% filter(p_val_adj < p_threshold)
print(head(markers))
write.csv(markers,paste0("markers_", name, ".csv"),row.names = FALSE) # nolint