library(Biobase)
library(Biobase)
library(limma)
rm(list=ls())
#########################################E-MEXP-3196##########################################
# 1）读入 eSet 文件
load("E-MEXP-3196.eSet.r")   # 路径按你自己的位置改
# 2）看看加载了什么对象
class(study)                  # 应该是 "NChannelSet"
fdata <- fData(study)
## 取出四个矩阵：前景 + 背景
G  <- assayData(study)$G   # green foreground
Gb <- assayData(study)$Gb  # green background
R  <- assayData(study)$R   # red foreground
Rb <- assayData(study)$Rb  # red background
## 组 RGList 做背景校正
RG <- new("RGList", list(G = G, Gb = Gb, R = R, Rb = Rb))
# normexp + offset 是 Agilent 很常用的
RG.bc <- backgroundCorrect(RG, method = "normexp", offset = 50)
# 取校正后的 G/R
G.bc <- RG.bc$G
R.bc <- RG.bc$R
# log2 + 芯片间分位数归一化（按通道看成 8 个样本）
G.log <- log2(G.bc + 1)
R.log <- log2(R.bc + 1)
# 对 G 和 R 一起做 between-array normalization
expr_all <- normalizeBetweenArrays(cbind(G.log, R.log), method = "quantile")
# 前 4 列是 Cy3(G)，后 4 列是 Cy5(R)
expr_G <- expr_all[, 1:ncol(G.log), drop = FALSE]
expr_R <- expr_all[, (ncol(G.log)+1):ncol(expr_all), drop = FALSE]
locus_raw <- as.character(fdata$Composite.Element.Database.Entry.locus.)

# 清洗 locus：处理 "GENE1 /// GENE2" 这种情况
locus_clean <- sapply(strsplit(locus_raw, "///|;|,"), function(x) trimws(x[1]))
locus_clean[locus_clean == ""] <- NA

anno_spot <- as.character(fdata$Composite.Element.Name)

anno_df <- data.frame(
  spot  = anno_spot,
  locus = locus_clean,
  stringsAsFactors = FALSE
)

# 去掉没有 spot 或 locus 的行
anno_df <- anno_df[!is.na(anno_df$spot) & anno_df$spot != "" &
                     !is.na(anno_df$locus) & anno_df$locus != "", ]

# 如果同一个 spot 有多行，先保留第一行即可
anno_df <- anno_df[!duplicated(anno_df$spot), ]

# 用 spot 作为 rownames，方便后面按行名索引
rownames(anno_df) <- anno_df$spot
spot_exp <- rownames(expr_G)   # "1","2","3",...,"22575"

# 这一步就是“按 exp 行名去 anno_df 里找”，返回的是向量
gene_symbol_full <- anno_df[spot_exp, "locus"]

length(gene_symbol_full); nrow(expr_G)
# 这里应该是相等
stopifnot(length(gene_symbol_full) == nrow(expr_G))
keep <- !is.na(gene_symbol_full) & gene_symbol_full != ""

expr_G_keep <- expr_G[keep, , drop = FALSE]
expr_R_keep <- expr_R[keep, , drop = FALSE]
symbol_keep <- gene_symbol_full[keep]

stopifnot(length(symbol_keep) == nrow(expr_G_keep))

G_gene <- rowsum(expr_G_keep, group = symbol_keep)
R_gene <- rowsum(expr_R_keep, group = symbol_keep)

expr_by_gene <- cbind(G_gene, R_gene)
dim(expr_by_gene)
head(rownames(expr_by_gene))

sample_names <- colnames(expr_by_gene)
pheno <- pData(study)
sample_Cy3 <- make.names(pheno$Source.Name.Cy3)
sample_Cy5 <- make.names(pheno$Source.Name.Cy5)
colnames(expr_by_gene) <- c(sample_Cy3, sample_Cy5)
status_Cy3 <- pheno$Characteristics.DiseaseState..Cy3
status_Cy5 <- pheno$Characteristics.DiseaseState..Cy5

outcome_Cy3 <- ifelse(status_Cy3 %in% c("normal", "Normal"), 0L, 1L)
outcome_Cy5 <- ifelse(status_Cy5 %in% c("normal", "Normal"), 0L, 1L)
Testing_class_3196 <- data.frame(
  sample  = c(sample_Cy3, sample_Cy5),
  outcome = c(outcome_Cy3, outcome_Cy5),
  Cohort  = "E-MEXP-3196",
  stringsAsFactors = FALSE
)
rownames(Testing_class_3196) <- Testing_class_3196$sample
table(Testing_class_3196$outcome)

save(expr_by_gene,Testing_class_3474,file="E-MEXP-3196_raw")
######################################################################################
###################################################GSE32127##########################################################
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
######################################################################################
###################################################GSE59446##########################################################
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
