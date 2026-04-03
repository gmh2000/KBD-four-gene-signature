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
setwd("D:\\project/KBD\\data_raw\\GSE186593/")
############
# lines <- readLines("GSE186593_series_matrix (1).txt/GSE186593_series_matrix (1).txt")
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
count=read.delim(file = "GSE186593_raw_counts_GRCh38.p13_NCBI.tsv/GSE186593_raw_counts_GRCh38.p13_NCBI.tsv",row.names = 1)
anno=read.delim(file = "Human.GRCh38.p13.annot.tsv/Human.GRCh38.p13.annot.tsv",row.names = 1)
meta=read.csv(file="GSE77298_meta.csv",row.names = 1)
# 保留 symbol 非 NA 的行
count$symbol <- anno$Symbol[match(rownames(count), rownames(anno))]
count$symbol_unique <- make.unique(count$symbol)
rownames(count) <- count$symbol_unique
count=count[,c(-11,-12)]
colnames(count) == rownames(meta)
GSE186593_count=count
GSE186593_meta=meta
save(GSE186593_count,GSE186593_meta, file = "GSE186593_raw.RData")
library(DESeq2)
library(tidyverse)
x <- DGEList(counts=count,group=factor(meta$Group))
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE, prior.count=2)
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
keep.exprs <- filterByExpr(x, group=factor(meta$Group))
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
lcpm <-lcpm[keep.exprs,]
library(pheatmap)
df=count[keep.exprs,]
annotation_col <- data.frame(
  type = meta$Group
)
rownames(annotation_col) <- colnames(df)
pdf("PCCB_pheatmap.pdf")
pheatmap(
  df,
  scale = "row",
  color = rainbow(7),
  annotation_col = annotation_col,
  cluster_rows = FALSE,
  show_rownames = FALSE
)
dev.off()

dds <- DESeqDataSetFromMatrix(
  countData = count, 
  colData   = meta, 
  design    = ~ Sex + Group
)

# 过滤低表达基因
dds <- dds[rowSums(counts(dds)) > 15, ]

library(sva)

# vst 变换用于 SVA
vsd <- vst(dds, blind=TRUE)
vsd_mat <- assay(vsd)
mod  <- model.matrix(~ Sex + Group, data=meta)   # full model
mod0 <- model.matrix(~ Sex, data=meta)           # null model
n.sv  <- num.sv(vsd_mat, mod, method="leek")
svobj <- sva(vsd_mat, mod, mod0, n.sv = n.sv)
SV <- svobj$sv
colnames(SV) <- paste0("SV", 1:ncol(SV))
meta_with_sv <- cbind(meta, SV)
dds_sv <- DESeqDataSetFromMatrix(
  countData = count,
  colData   = meta_with_sv,
  design    = as.formula(
    paste("~ Sex +", paste(colnames(SV), collapse = " + "), "+ Group")
  )
)

dds_sv <- dds_sv[rowSums(counts(dds_sv)) > 15, ]
dds_sv <- DESeq(dds_sv)
res <- results(dds_sv, contrast = c("Group", "KBD", "Normal"))
res <- res[order(res$padj), ]
res <- na.omit(res)
head(res)

padj_cut <- 0.05
lfc_cut  <- 1.5  # log2 fold change ≥ 1.5

deg_all <- res %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  filter(padj < padj_cut, abs(log2FoldChange) >= lfc_cut)

deg_up   <- deg_all %>% filter(log2FoldChange >= lfc_cut)
deg_down <- deg_all %>% filter(log2FoldChange <= -lfc_cut)
save.image(file="GSE186593_DEG_all.RData")
save(dds_sv,res,deg_all, file="GSE186593_DEG1.5.RData")

top_up   <- head(deg_up, 5)
top_down <- head(deg_down, 5)
label_genes <- c(top_up$Gene, top_down$Gene)

labels <- ifelse(rownames(res) %in% label_genes, rownames(res), NA)
## -------------------------
## 2. 发表级 EnhancedVolcano 参数
## -------------------------


df <- as.data.frame(res)

keyvals <- rep("grey80", nrow(df))
names(keyvals) <- rep("NS", nrow(df))

# 显著上调：padj<0.05 且 log2FC >= 1
idx_up <- df$padj < padj_cut & df$log2FoldChange >=  lfc_cut
keyvals[idx_up] <- "#D81B60"
names(keyvals)[idx_up] <- "Up"
library(EnhancedVolcano)
# 显著下调：padj<0.05 且 log2FC <= -1
idx_down <- df$padj < padj_cut & df$log2FoldChange <= -lfc_cut
keyvals[idx_down] <- "#1E88E5"
names(keyvals)[idx_down] <- "Down"
png("Fig2_GSE186593_Volcano_Top5_Publication1.5.png",
    width = 4000, height = 4000, res = 600)
EnhancedVolcano(
  df,
  lab              = labels,          # 你之前定义好的标签，可用 NA 控制只标 top 基因
  x                = "log2FoldChange",
  y                = "padj",
  pCutoff          = padj_cut,
  FCcutoff         = lfc_cut,
  
  ## 关键：自定义每个点的颜色（只三类：NS/Down/Up）
  colCustom        = keyvals,
  
  xlim             = c(-10, 10),
  shape            = 16,
  colAlpha         = 0.9,
  
  drawConnectors   = TRUE,
  colConnectors    = "black",
  widthConnectors  = 0.5,
  max.overlaps     = Inf,
  boxedLabels      = TRUE,
  caption          = NULL,   
  title            = "Differential Expression",
  subtitle         = "DESeq2, padj < 0.05, |log2FC| ≥ 1",
  xlab             = bquote(Log[2]~Fold~Change),
  ylab             = bquote(-Log[10]~padj),
  
  legendPosition   = "top",
  legendLabSize    = 12,
  legendIconSize   = 4.5,
  gridlines.major  = FALSE,
  gridlines.minor  = FALSE
)

dev.off()

pdf("Fig2_GSE186593_Volcano_Top5_Publication1.5.pdf",
    width = 8, height = 8)  # 将像素转换为英寸
EnhancedVolcano(
  df,
  lab              = labels,          # 你之前定义好的标签，可用 NA 控制只标 top 基因
  x                = "log2FoldChange",
  y                = "padj",
  pCutoff          = padj_cut,
  FCcutoff         = lfc_cut,
  
  ## 关键：自定义每个点的颜色（只三类：NS/Down/Up）
  colCustom        = keyvals,
  
  xlim             = c(-10, 10),
  shape            = 16,
  colAlpha         = 0.9,
  
  drawConnectors   = TRUE,
  colConnectors    = "black",
  widthConnectors  = 0.5,
  max.overlaps     = Inf,
  boxedLabels      = TRUE,
  
  title            = "Differential Expression",
  subtitle         = "DESeq2, padj < 0.05, |log2FC| ≥ 1",
  xlab             = bquote(Log[2]~Fold~Change),
  ylab             = bquote(-Log[10]~padj),
  caption          = NULL,          # 隐藏所有标题
  captionLabSize   = 0    ,          # 将标签大小设为0
  legendPosition   = "top",
  legendLabSize    = 12,
  legendIconSize   = 4.5,
  gridlines.major  = FALSE,
  gridlines.minor  = FALSE
)

dev.off()




library(pheatmap)
vst_mat <- assay(vsd)
top50 <- deg_all$Gene[1:50]
mat50 <- vst_mat[top50, ]
ann <- data.frame(Group = meta$Group)
rownames(ann) <- rownames(meta)
ann_colors <- list(Group = c(KBD="#FF7043", Normal="#29B6F6"))
png("Top50_heatmap_GSE186593_1.5.png", width = 2500, height = 2500, res = 300)
pheatmap(
  mat50,
  scale="row",
  clustering_method="complete",
  annotation_col = ann,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#4575b4","white","#D81B60"))(256),
  fontsize_row = 7
)
dev.off()

pdf("Top50_heatmap_GSE186593_1.5.pdf",
    width = 8, height = 7)  # 将像素转换为英寸
pheatmap(
  mat50,
  scale="row",
  clustering_method="complete",
  annotation_col = ann,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#4575b4","white","#D81B60"))(256),
  fontsize_row = 7
)
dev.off()


save(dds,file="GSE186593_dds.RData")
save.image(file="GSE186593_all.RData")
