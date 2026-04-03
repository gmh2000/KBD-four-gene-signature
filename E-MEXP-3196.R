library(Biobase)
library(Biobase)
library(limma)
rm(list=ls())
setwd("D:\\project\\KBD\\data_raw\\E-MEXP-3196/")
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
Testing_class_3474 <- data.frame(
  sample  = c(sample_Cy3, sample_Cy5),
  outcome = c(outcome_Cy3, outcome_Cy5),
  Cohort  = "E-MEXP-3474",
  stringsAsFactors = FALSE
)
rownames(Testing_class_3474) <- Testing_class_3474$sample
table(Testing_class_3474$outcome)

save(expr_by_gene,Testing_class_3474,file="E-MEXP-3196_raw")

