rm(list=ls())
setwd("D:\\课题/KBD/data_raw\\GSE32127\\GSE32127_RAW")
pkgs <- c("limma", "stringr", "pheatmap", "EnhancedVolcano", "dplyr")
for(p in pkgs){ if(!requireNamespace(p, quietly=TRUE)) install.packages(p) }
library(limma); library(stringr); library(pheatmap); library(EnhancedVolcano); library(dplyr)
options(stringsAsFactors = FALSE)
## 1. 列出 GPR 文件
gpr_dir <- "D:\\课题/KBD/data_raw\\GSE32127\\GSE32127_RAW"
files <- list.files(gpr_dir, pattern="\\.gpr(\\.gz)?$", full.names = TRUE)
stopifnot(length(files) > 0)
## （可选）如果你能从文件名看出 Z=KBD, Y=CTRL，则构建 targets/isSwap：
## 这里视情况启用，否则暂不使用
 parse_target <- function(f){
   nm <- basename(f)
   left  <- toupper(stringr::str_match(nm, "^([A-Z])\\d*_vs_")[,2])
   right <- toupper(stringr::str_match(nm, "_vs_([A-Z])\\d*")[,2])
   data.frame(FileName=f, Left=left, Right=right, stringsAsFactors = FALSE)
 }
 targets <- do.call(rbind, lapply(files, parse_target))
 grp_map <- c(Z="KBD", Y="CTRL")
 targets$R.Group <- grp_map[targets$Left]
 targets$G.Group <- grp_map[targets$Right]
 targets$isSwap <- with(targets, R.Group=="CTRL" & G.Group=="KBD")
## 2. 预清洗 GPR（去掉带中文路径的 GalFile/FileName 行）
clean_gpr <- function(f, outdir="data_raw/GSE32127_GPR_clean"){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  con <- if (grepl("\\.gz$", f, ignore.case=TRUE)) gzfile(f, open="rt") else file(f, open="rt")
  txt <- readLines(con, warn = FALSE, encoding = "Latin1")
  close(con)
  txt <- iconv(txt, from = "Latin1", to = "UTF-8", sub = "byte")
  drop <- grepl("^(GalFile|FileName)\\s*=", txt, ignore.case = TRUE)
  txt <- txt[!drop]
  out <- file.path(outdir, paste0(basename(f), ".clean.gpr"))
  writeLines(txt, out, useBytes = TRUE)
  out
}
clean_files <- vapply(files, clean_gpr, FUN.VALUE = character(1))

## 3. 读入两色数据 + 归一化
cols <- list(R="F635 Median", G="F532 Median", Rb="B635 Median", Gb="B532 Median")  # 如用 Mean 就改成 Mean
RG <- read.maimages(files = clean_files, source = "genepix", columns = cols)
RGb <- backgroundCorrect(RG, method="normexp", offset=20)
MA  <- normalizeWithinArrays(RGb, method="loess")
MA  <- normalizeBetweenArrays(MA, method="Aquantile")

## （可选）若有染料互换，则在规范化后翻转 M 值
## MA$M[, targets$isSwap] <- -MA$M[, targets$isSwap]

## 4. 先把探针 ID 写入行名，再做过滤 / MAf

# 这里提前设行名（关键）probe_id_candidates <- c("ProbeName","ID","Name","SystematicName","ProbeID","Probe Id")
pick_probe_col <- function(genes_df){
  cand <- probe_id_candidates[probe_id_candidates %in% colnames(genes_df)]
  if(length(cand)) cand[1] else NA
}

probe_col_MA <- pick_probe_col(MA$genes)
if(is.na(probe_col_MA)) {
  probe_col_RG <- pick_probe_col(RG$genes)
  if(is.na(probe_col_RG)) stop("未在 MA$genes 或 RG$genes 找到探针ID列；请检查GPR列名。")
  MA$genes[[probe_col_RG]] <- RG$genes[[probe_col_RG]]
  probe_col_MA <- probe_col_RG
}
rownames(MA$M) <- make.unique(as.character(MA$genes[[probe_col_MA]]))
rownames(MA$A) <- rownames(MA$M)
## 5. 过滤低变异探针（IQR 前 50%）
iqr_all <- apply(MA$M, 1, IQR, na.rm=TRUE)

## 计算 30% 分位数作为过滤阈值（保留前 70%）
th <- quantile(iqr_all, 0.40, na.rm = TRUE)

keep <- iqr_all > th
sum(keep)  # 查看保留多少探针

MAf <- MA
MAf$M     <- MA$M[keep, , drop=FALSE]
MAf$A     <- MA$A[keep, , drop=FALSE]
MAf$genes <- MA$genes[keep, , drop=FALSE]
probe_ids <- as.character(MAf$genes[[probe_col_MA]])
probe_ids[is.na(probe_ids) | probe_ids==""] <- paste0("NA_", seq_len(sum(is.na(probe_ids) | probe_ids=="")))
probe_ids <- make.unique(probe_ids)
rownames(MAf$M) <- probe_ids
rownames(MAf$A) <- probe_ids
MAf$genes[[probe_col_MA]] <- probe_ids
## 6. 保障列名/设计矩阵
if (is.null(colnames(MAf$M))) {
  colnames(MAf$M) <- paste0("Array", seq_len(ncol(MAf$M)))
}
if (is.null(colnames(MAf$A))) {
  colnames(MAf$A) <- colnames(MAf$M)
}
if (any(duplicated(colnames(MAf$M)))) {
  colnames(MAf$M) <- make.unique(colnames(MAf$M))
  colnames(MAf$A) <- colnames(MAf$M)
}
design <- matrix(1, ncol=1, nrow=ncol(MAf$M))
colnames(design) <- "KBDvsCTRL"
rownames(design) <- colnames(MAf$M)
## 7. arrayWeights + limma 
aw <- arrayWeights(as.matrix(MAf$M), design = design)
fitP  <- lmFit(MAf$M, design, weights = aw)
fitP  <- eBayes(fitP, trend = TRUE, robust = TRUE)
tt1   <- topTable(fitP, coef=1, number=Inf)

# 让行名也等于探针 ID
tt1$Probe= rownames(tt1)  
# 此时：rownames(tt) 和 tt$Probe 都是探针 ID
head(tt1[, c("Probe","logFC","P.Value","adj.P.Val")])
# 2) 读取本地 GPL 注释并建立 Probe→Symbol 映射
gpl_file <- "GPL7264-9589.txt"              # 路径改为你的
annot <- read.delim(gpl_file, header = TRUE, quote = "", comment.char = "#",
                    stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
id_col  <- c("ID","ProbeName","SystematicName","Probe ID","PROBE_ID")
sym_col <- c("Gene Symbol","GENE_SYMBOL","Gene symbol","Symbol","SYMBOL","GENE_SYMBOLS")
id_col  <- id_col[id_col %in% colnames(annot)][1]
sym_col <- sym_col[sym_col %in% colnames(annot)][1]
stopifnot(!is.na(id_col), !is.na(sym_col))
map <- annot[, c(id_col, sym_col)]
colnames(map) <- c("Probe","Symbol_raw")
# 3) 清理多符号与空符号
split_symbol <- function(x){
  if (is.na(x) || x=="") return(NA_character_)
  parts <- unlist(strsplit(x, "\\s*///\\s*|\\s*;\\s*|\\s*,\\s*|\\s*\\|\\s*"))
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(NA_character_)
  parts[1]
}
map$Symbol <- vapply(map$Symbol_raw, split_symbol, FUN.VALUE = character(1))
map <- subset(map, !is.na(Symbol) & Symbol != "")

# 4) 合并注释到探针层结果
tt1$Symbol <- map$Symbol[match(tt1$Probe, map$Probe)]
tt1 <- tt1[!is.na(tt1$Symbol) & tt1$Symbol!="", ]

# 5) 每个基因选择“最显著探针”（P 值最小；并保留方向/效应量）
tt1 <- tt1[order(tt1$P.Value), ]
best_by_gene <- tt1[!duplicated(tt1$Symbol), ]          # 去重后每基因一行
rownames(best_by_gene) <- best_by_gene$Symbol
deg = best_by_gene
# --- Fig2: 火山图 ---
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


padj_cut <- 0.05
lfc_cut  <- 0.58  # 若想更严格就用 0.58 (~1.5x) 或 1.0 (~2x)
library(dplyr); library(tibble)
get_deg_sets <- function(res, padj = padj_cut, lfc = lfc_cut) {
  res <- res %>% as.data.frame() %>% rownames_to_column("Gene")
  up   <- res %>% filter(adj.P.Val < padj, logFC >=  lfc) %>% arrange(adj.P.Val)
  down <- res %>% filter(adj.P.Val < padj, logFC <= -lfc) %>% arrange(adj.P.Val)
  list(up = up, down = down, all = res %>% filter(adj.P.Val < padj, abs(logFC) >= lfc))
}
deg_blood <- get_deg_sets(best_by_gene, padj_cut, lfc_cut)
deg_blood$all$Gene


# -- Fig2: 热图（Top 50） ---

top50 <- deg[order(deg$adj.P.Val), ]
top50 <- head(top50, 50)
genes_top50 <- rownames(top50)
# MAf$M 行名是 probe，需要映射到 symbol
probe2symbol <- tt1[, c("Probe", "Symbol")]  # 或者你已有 map

# 取 Top50 所属的 probe
top50_probe <- probe2symbol$Probe[match(genes_top50, probe2symbol$Symbol)]
top50_probe <- na.omit(top50_probe)

# 提取表达矩阵
mat <- MAf$M[top50_probe, , drop=FALSE]

# 行名改成基因名
rownames(mat) <- probe2symbol$Symbol[match(rownames(mat), probe2symbol$Probe)]
ann <- data.frame(Group = targets$R.Group)   # 或者 meta$Group
rownames(ann) <- colnames(mat)

ann_colors <- list(
  Group = c(KBD="#FF7043", CTRL="#29B6F6")
)

library(pheatmap)

png("Fig3_Heatmap_Top50_TwoColor_KBD.png",
    width = 2500, height = 3200, res = 320)

pheatmap(
  mat,
  scale = "row",                 # 按行z-score，凸显差异
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7.0,
  fontsize_col = 10,
  
  annotation_col = ann,
  annotation_colors = ann_colors,
  
  ## ---- 发表级配色（蓝-白-红） ----
  color = colorRampPalette(c("#4575B4", "white", "#D73027"))(200),
  
  border_color = NA,
  main = "Top 50 Differentially Expressed Genes (Two-colour microarray)"
)

dev.off()

save.image()



deg_blood$down


message("✓ GPR 两色差异分析完成。结果：results/tables/Fig2_DEG_GSE32127_gpr.csv；图：results/fig/")


