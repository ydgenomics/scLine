library(readr)
library(openxlsx)

# # 1. 中英文注释（按需增删）
# col_comment <- c(
#   "## 查询基因ID",
#   "## 种子直系同源蛋白",
#   "## E值",
#   "## 得分",
#   "## eggNOG直系同源组",
#   "## 最高注释级别",
#   "## COG功能分类",
#   "## 功能描述",
#   "## 推荐基因名",
#   "## GO注释",
#   "## 酶编号",
#   "## KEGG Orthology",
#   "## KEGG通路",
#   "## KEGG模块",
#   "## KEGG反应",
#   "## KEGG反应类",
#   "## BRITE层级",
#   "## KEGG转运分类",
#   "## 碳水化合物酶",
#   "## BiGG反应",
#   "## PFAM结构域"
# )

# 2. 读入文件（tab 分隔，带表头）
df <- read_tsv("/data/work/scline/input/Galaxy6.tabular")
names(df)[1] <- "query"

# # 3. 构造两行注释头（## 开头）
# comment_df <- matrix("##", nrow = 2, ncol = ncol(df), byrow = TRUE)
# colnames(comment_df) <- names(df)

# 4. 写出 Excel（注释头 + 原表头 + 数据）
wb <- createWorkbook()
addWorksheet(wb, "Annotation")

# ① 第 1-2 行：空注释
for (col in 1:ncol(df)) {
  writeData(wb, "Annotation", x = "##", startCol = col, startRow = 1)
  writeData(wb, "Annotation", x = "##", startCol = col, startRow = 2)
}

# ② 第 3 行：原始列名（分散写）
for (col in 1:ncol(df)) {
  writeData(wb, "Annotation", x = names(df)[col], startCol = col, startRow = 3)
}

# ③ 第 4 行起：数据
writeData(wb, "Annotation", x = df, startRow = 4, colNames = FALSE)

# 5. 保存
saveWorkbook(wb, "eggNOG_annotation.xlsx", overwrite = TRUE)
print("已生成：eggNOG_annotation.xlsx（前两行为##注释）")

# emapper <- read_excel("eggNOG_annotation.xlsx", skip = 2) # check it could be loaded