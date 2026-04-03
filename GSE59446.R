library(affy)
library(limma)
library(genefilter)
library(WGCNA)
library(hgu133a.db)
library(hgu133plus2.db)
library(hgu95av2.db)
library(hgu219.db)
library(annotate)
library(ggplot2)
library(sva)
library(GEOquery)
library(oligo)
library(huex10sttranscriptcluster.db)
library(dplyr)
library(tibble)
library(stringr)
library(dplyr)
library(tidyr)      
library(CIBERSORT)
library(dplyr); library(matrixStats); library(WGCNA)
library(ggplot2); library(pheatmap); library(edgeR)
rm(list=ls())
setwd("D:\\课题/KBD\\data_raw\\GSE59446/")
# lines <- readLines("GSE59446_series_matrix.txt")
# meta_lines <- grep("^!Sample_", lines, value = TRUE)
# # 所有样本名称（假设每个字段都应该有相同样本顺序）
# sample_ids <- strsplit(grep("^!Sample_title", meta_lines, value = TRUE), "\t")[[1]][-1]
# # 提取所有字段名称
# get_field_name <- function(line) {
#   strsplit(sub("^!Sample_", "", line), "\t")[[1]][1]
# }
# fields <- unique(sapply(meta_lines, get_field_name))
# # 初始化存储列表
# sample_info <- list()
# for (field in fields) {
#   lines_for_field <- grep(paste0("^!Sample_", field, "\t"), meta_lines, value = TRUE)
#   if (length(lines_for_field) == 0) next
#   # 拆分每行值并合并（注意缺失补 NA）
#   value_matrix <- do.call(rbind, lapply(lines_for_field, function(line) {
#     vals <- strsplit(line, "\t")[[1]][-1]
#     length(vals) <- length(sample_ids)  # 补齐为标准长度
#     return(vals)
#   }))
#   values <- apply(value_matrix, 2, function(x) {
#     paste(na.omit(x), collapse = "; ")
#   })
#   values <- gsub("^\"|\"$", "", values)
#   sample_info[[field]] <- values
# }
# # 转为 DataFrame
# sample_df <- as.data.frame(sample_info, stringsAsFactors = FALSE)
# rownames(sample_df) <- sample_ids
# write.csv(sample_df,"GSE77298_meta.csv",quote =FALSE)

rm(list=ls())
count1=read.delim(file = "GSE59446_raw_data_KBD.txt/GSE59446_raw_data_KBD.txt")
count1=count1[,-1]
count2=read.delim(file = "GSE59446_raw_data_control.txt/GSE59446_raw_data_control.txt")
count2=count2[,-1]
merged <- merge(count1, count2, by="name", all=TRUE)
rownames(merged) <- merged$name
merged=merged[,-1]
meta=read.csv("GSE59446_series_matrix.txt/GSE56446_meta.csv",row.names = 1)
GSE59446_count=merged
GSE59446_meta=meta
save(GSE59446_count,GSE59446_meta, file = "GSE59446_.RData")
expr <- log2(merged + 1)
rownames(meta)=gsub("-",".",rownames(meta))
meta <- meta[colnames(expr), ]
design <- model.matrix(~ 0 + Group + Age + Sex, data=meta)
fit <- lmFit(expr, design)
contr <- makeContrasts(KBDvsControl = GroupKBD - GroupControl, levels=design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
res <- topTable(fit2, coef="KBDvsControl", number=Inf)

padj_cut <- 0.05
lfc_cut  <- 0.58  # 若想更严格就用 0.58 (~1.5x) 或 1.0 (~2x)
library(dplyr); library(tibble)
get_deg_sets <- function(res, padj = padj_cut, lfc = lfc_cut) {
  res <- res %>% as.data.frame() %>% rownames_to_column("Gene")
  up   <- res %>% filter(adj.P.Val < padj, logFC >=  lfc) %>% arrange(adj.P.Val)
  down <- res %>% filter(adj.P.Val < padj, logFC <= -lfc) %>% arrange(adj.P.Val)
  list(up = up, down = down, all = res %>% filter(adj.P.Val < padj, abs(logFC) >= lfc))
}
deg_blood <- get_deg_sets(res, padj_cut, lfc_cut)
deg_blood$up$Gene
deg=res
up5 <- deg[deg$logFC > 0, ]
up5 <- up5[order(up5$P.Value), ]   
up5 <- head(up5, 5)

down5 <- deg[deg$logFC < 0, ]
down5 <- down5[order(down5$P.Value), ]
down5 <- head(down5, 5)

label_genes <- c(rownames(up5), rownames(down5))

## 标签对象：非 top 基因全部 NA
labels <- ifelse(rownames(deg) %in% label_genes, rownames(deg), NA)


## -------------------------
## 2. 发表级 EnhancedVolcano 参数
## -------------------------

png("Fig2_GSE32127_Volcano_Top5_Publication.png",
    width = 2500, height = 2000, res = 300)
EnhancedVolcano(
  deg,
  lab              = labels,
  x                = "logFC",
  y                = "adj.P.Val",       # 或 adj.P.Val
  pCutoff          = 0.05,
  FCcutoff         = 0.58,
  
  ## ---- 关键：发表级配色 ----
  col              = c("grey70", "#1565C0", "#D32F2F", "#8E24AA"),
  colAlpha         = 0.9,
  shape            = 16,              # 实心点
  border           = "partial",
  
  ## ---- 标签与连接线优化（发表级） ----
  drawConnectors   = TRUE,
  colConnectors    = "black",
  widthConnectors  = 0.6,
  max.overlaps     = Inf,             # 必须，否则 top5 不一定显示
  
  labSize          = 5.2,
  labCol           = "black",
  boxedLabels      = TRUE,           # 如果想用白底，可以改成 TRUE
  
  ## --- 标题与坐标轴：更适合期刊风格 ----
  title            = "Differential Expression in KBD (GSE32127)",
  subtitle         = "limma(adj.P.Val < 0.05,|log2FC|≥0.58)",
  xlab             = bquote(Log[2]~"Fold Change"),
  ylab             = bquote(-Log[10]~italic(adj.P.Val)),
  
  titleLabSize     = 18,
  subtitleLabSize  = 13,
  axisLabSize      = 14,
  captionLabSize   = 10,
  
  ## ---- 图例优化 ----
  legendPosition   = "top",
  legendLabSize    = 12,
  legendIconSize   = 4.5,
  
  ## ---- 网格线关闭（更简洁） ----
  gridlines.major  = FALSE,
  gridlines.minor  = FALSE
)

dev.off()











top50 <- res[order(res$adj.P.Val), ]
top50 <- head(top50, 50)

genes_top50 <- rownames(top50)
mat <- expr[genes_top50, , drop = FALSE]
ann <- data.frame(Group = meta$Group)
rownames(ann) <- colnames(mat)

ann_colors <- list(
  Group = c(
    KBD     = "#FF7043",   # 橙色
    Control = "#29B6F6"    # 天蓝色
  )
)
library(pheatmap)

png("GSE59446_Top50_heatmap.png", width = 2400, height = 3000, res = 320)

pheatmap(
  mat,
  scale = "row",                     # 按行 z-score，更易观察差异
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  fontsize_row = 7,
  fontsize_col = 10,
  
  annotation_col = ann,
  annotation_colors = ann_colors,
  
  color = colorRampPalette(c("#4575B4", "white", "#D73027"))(200),  # 蓝〜白〜红（Nature 风）
  
  border_color = NA,
  main = "Top 50 Differential Genes (GSE59446 KBD vs Control)"
)

dev.off()
