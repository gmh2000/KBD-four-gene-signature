##==================================================
## GSE32127：2-color → 单样本表达矩阵（基因 × 样本）
## 通道：G = one sample, R = another sample
##==================================================

work.path <- "D:/project/KBD/data_proc/ML/InputData/"
dir.create(work.path, recursive = TRUE, showWarnings = FALSE)

rm(list = ls())

##--------------------------------------------------
## 0. 环境与基本设置
##--------------------------------------------------
raw_dir <- "D:/project/KBD/data_raw/GSE32127/GSE32127_RAW"
setwd(raw_dir)

pkgs <- c("limma", "stringr", "dplyr")
for(p in pkgs){
  if(!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(limma)
library(stringr)
library(dplyr)

options(stringsAsFactors = FALSE)

##--------------------------------------------------
## 1. 定义一个基因层面聚合函数（按行对齐的 map_vec）
##--------------------------------------------------
aggregate_to_gene_vec <- function(expr_mat, gene_vec, fun = stats::median){
  # expr_mat: probe × sample
  # gene_vec: 与 probe 行一一对应的基因名（向量）
  stopifnot(length(gene_vec) == nrow(expr_mat))
  
  keep <- !is.na(gene_vec) & gene_vec != ""
  expr_keep <- expr_mat[keep, , drop = FALSE]
  gene_keep <- gene_vec[keep]
  
  # 用 rowsum 聚合：每个基因一行（默认求和，后面改成 median）
  sum_mat <- rowsum(expr_keep, group = gene_keep)
  
  if(identical(fun, base::sum)) {
    return(sum_mat)
  } else if(identical(fun, stats::median)) {
    # rowsum 得到的是“和”，这里再按基因计算 median
    df <- as.data.frame(expr_keep)
    df$Symbol <- gene_keep
    agg <- aggregate(. ~ Symbol, data = df, FUN = stats::median)
    rn <- agg$Symbol
    agg$Symbol <- NULL
    mat <- as.matrix(agg)
    rownames(mat) <- rn
    return(mat)
  } else {
    # 通用写法：按基因拆分，再用 fun
    df <- as.data.frame(expr_keep)
    df$Symbol <- gene_keep
    agg <- aggregate(. ~ Symbol, data = df, FUN = fun)
    rn <- agg$Symbol
    agg$Symbol <- NULL
    mat <- as.matrix(agg)
    rownames(mat) <- rn
    return(mat)
  }
}

##--------------------------------------------------
## 2. 列出 GPR 文件 & 解析文件名里的组别信息
##--------------------------------------------------
files <- list.files(raw_dir, pattern = "\\.gpr(\\.gz)?$", full.names = TRUE)
stopifnot(length(files) > 0)

# 从文件名解析 Left/Right（例如 Z1_vs_Y2.gpr）
parse_target <- function(f){
  nm <- basename(f)
  left  <- toupper(str_match(nm, "^([A-Z])\\d*_vs_")[,2])
  right <- toupper(str_match(nm, "_vs_([A-Z])\\d*")[,2])
  data.frame(FileName = f, Left = left, Right = right, stringsAsFactors = FALSE)
}
targets <- do.call(rbind, lapply(files, parse_target))

# 这里按你之前的约定：Z=KBD, Y=CTRL
grp_map <- c(Z = "KBD", Y = "CTRL")
targets$R.Group <- grp_map[targets$Left]   # 红通道对应 Left
targets$G.Group <- grp_map[targets$Right]  # 绿通道对应 Right

##--------------------------------------------------
## 3. 预清洗 GPR（去掉含中文路径的 GalFile/FileName 行）
##--------------------------------------------------
clean_gpr <- function(f, outdir = "GSE32127_GPR_clean"){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  con <- if (grepl("\\.gz$", f, ignore.case = TRUE)) gzfile(f, open = "rt") else file(f, open = "rt")
  txt <- readLines(con, warn = FALSE, encoding = "Latin1")
  close(con)
  txt <- iconv(txt, from = "Latin1", to = "UTF-8", sub = "byte")
  drop <- grepl("^(GalFile|FileName)\\s*=", txt, ignore.case = TRUE)
  txt  <- txt[!drop]
  out  <- file.path(outdir, paste0(basename(f), ".clean.gpr"))
  writeLines(txt, out, useBytes = TRUE)
  out
}
clean_files <- vapply(files, clean_gpr, FUN.VALUE = character(1))

##--------------------------------------------------
## 4. 读入两色数据 + 背景校正 + 通道归一化
##--------------------------------------------------
# 按 Genepix 标准列名，这里用 Median，你之前也这样设
cols <- list(R = "F635 Median", G = "F532 Median",
             Rb = "B635 Median", Gb = "B532 Median")

RG <- read.maimages(files = clean_files, source = "genepix", columns = cols)

# 背景校正
RG.bc <- limma::backgroundCorrect(RG, method = "normexp")  # 如果 offset 报错就先不加 offset

# 取校正后的 R/G
R.bc <- RG.bc$R   # probe × array
G.bc <- RG.bc$G

# log2 + 对所有通道一起做分位数归一化（8 个样本）
R.log <- log2(R.bc + 1)
G.log <- log2(G.bc + 1)

expr_all <- normalizeBetweenArrays(cbind(G.log, R.log), method = "quantile")
# 前 ncol(G.log) 列是 G，后 ncol(R.log) 列是 R
expr_G <- expr_all[, 1:ncol(G.log), drop = FALSE]
expr_R <- expr_all[, (ncol(G.log) + 1):ncol(expr_all), drop = FALSE]

##--------------------------------------------------
## 5. 统一行名为探针 ID
##--------------------------------------------------
probe_id_candidates <- c("ProbeName", "ID", "Name", "SystematicName", "ProbeID", "Probe Id")
pick_probe_col <- function(genes_df){
  cand <- probe_id_candidates[probe_id_candidates %in% colnames(genes_df)]
  if (length(cand)) cand[1] else NA
}
probe_col <- pick_probe_col(RG.bc$genes)
if (is.na(probe_col)) stop("在 RG$genes 中未找到探针 ID 列，请检查 GPR 的基因注释列名。")

probe_ids <- make.unique(as.character(RG.bc$genes[[probe_col]]))
stopifnot(nrow(expr_G) == length(probe_ids))

rownames(expr_G) <- probe_ids
rownames(expr_R) <- probe_ids

# 合并成 probe × (2 * arrays)
expr_probe <- cbind(expr_G, expr_R)

# 样本名（来自 GPR 文件名）
raw_names   <- basename(files)
array_names <- gsub(".gpr.*$", "", raw_names)

colnames(expr_G) <- paste0(array_names, "_G")
colnames(expr_R) <- paste0(array_names, "_R")
colnames(expr_probe) <- c(colnames(expr_G), colnames(expr_R))

##--------------------------------------------------
## 6. 读 GPL 注释，构建 Probe→Symbol 映射并按行对齐
##--------------------------------------------------
gpl_file <- "GPL7264-9589.txt"   # 路径按实际修改（建议放在 raw_dir 或绝对路径）
annot <- read.delim(gpl_file, header = TRUE, quote = "", comment.char = "#",
                    stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)

id_col  <- c("ID", "ProbeName", "SystematicName", "Probe ID", "PROBE_ID")
sym_col <- c("Gene Symbol", "GENE_SYMBOL", "Gene symbol", "Symbol", "SYMBOL", "GENE_SYMBOLS")

id_col  <- id_col[id_col %in% colnames(annot)][1]
sym_col <- sym_col[sym_col %in% colnames(annot)][1]
stopifnot(!is.na(id_col), !is.na(sym_col))

map <- annot[, c(id_col, sym_col)]
colnames(map) <- c("Probe", "Symbol_raw")

# 清理多基因 / 空 Symbol
split_symbol <- function(x){
  if (is.na(x) || x == "") return(NA_character_)
  parts <- unlist(strsplit(x, "\\s*///\\s*|\\s*;\\s*|\\s*,\\s*|\\s*\\|\\s*"))
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(NA_character_)
  parts[1]
}
map$Symbol <- vapply(map$Symbol_raw, split_symbol, FUN.VALUE = character(1))

# 构造一个与 expr_probe 行一一对应的 gene_vec
# 用 rownames(expr_probe) 去匹配 map$Probe
gene_vec <- map$Symbol[match(rownames(expr_probe), map$Probe)]
length(gene_vec); nrow(expr_probe)
stopifnot(length(gene_vec) == nrow(expr_probe))

##--------------------------------------------------
## 7. 按基因聚合为 gene × sample 矩阵
##--------------------------------------------------
expr_gene <- aggregate_to_gene_vec(expr_probe, gene_vec, fun = stats::median)
dim(expr_gene)
head(rownames(expr_gene))

##--------------------------------------------------
## 8. 构建 phenotype：每个通道一个样本（0/1 outcome）
##--------------------------------------------------
# G 通道（右样本）：CTRL / KBD
outcome_G <- ifelse(targets$G.Group == "KBD", 1L,
                    ifelse(targets$G.Group == "CTRL", 0L, NA_integer_))
# R 通道（左样本）
outcome_R <- ifelse(targets$R.Group == "KBD", 1L,
                    ifelse(targets$R.Group == "CTRL", 0L, NA_integer_))

sample_G <- colnames(expr_G)   # array_names_G
sample_R <- colnames(expr_R)   # array_names_R

Testing_class_32127 <- data.frame(
  sample  = c(sample_G, sample_R),
  outcome = c(outcome_G, outcome_R),
  Cohort  = "GSE32127",
  stringsAsFactors = FALSE
)
rownames(Testing_class_32127) <- Testing_class_32127$sample

table(Testing_class_32127$outcome, useNA = "ifany")
save(expr_gene,Testing_class_32127,file="GSE32127.RData")
##--------------------------------------------------
## 9. 导出到 InputData 目录
##--------------------------------------------------
setwd(work.path)

