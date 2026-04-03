

############### 输入：你的数据对象 ################
# exp_59446    # GSE59446 表达矩阵（行=基因，列=样本）
# meta_59446   # GSE59446 分组信息（含 Group 列）
# exp_32127    # GSE32127 表达矩阵（行=基因，列=样本）
# meta_32127   # GSE32127 分组信息（含 Group 列）
# genes25      # 25 个基因组成的字符向量

############### STEP 1：筛选 25 基因 ################
getwd()
rm(list=ls())
load("../../../data_raw/GSE186593/GSE186593_DEG1.5.RData")
load("../../../data_raw/GSE186593/WGCNA_hub_genes.RData")
load("../../../data_raw/GSE59446/GSE59446_raw.RData")
load("../../../data_proc/ML/InputData/E-MEXP-3196_raw.RData")
load("../../../data_proc/ML/InputData/GSE32127.RData")
gene26=read.table(file="../../../data_proc/gene28.txt",header=T)
genes_use <- intersect(deg_all$Gene, genes_in_module)
gene24=intersect(gene26$gene, rownames(GSE59446_count))
exp_59446=GSE59446_count
meta_59446=GSE59446_meta
genes_use <- intersect(gene24, rownames(exp_59446))
exp_59446 <- exp_59446[genes_use, , drop = FALSE]
expr_by_gene=as.data.frame(expr_by_gene)
expr_gene=as.data.frame(expr_gene)
expr_by_gene <- expr_by_gene[genes_use, , drop = FALSE]
expr_gene <- expr_gene[genes_use, , drop = FALSE]
train_exp  <- exp_59446
genes_use <- intersect(rownames(expr_by_gene), rownames(expr_gene))
genes_use <- genes_use[!genes_use %in% c("NA", "NA.1")]
expr_by_gene <- expr_by_gene[genes_use, , drop = FALSE]
expr_gene <- expr_gene[genes_use, , drop = FALSE]
combined_expr <- cbind(expr_by_gene, expr_gene)
test_exp  =   combined_expr
cat("训练集基因数：", nrow(train_exp), "\n")
cat("测试集基因数：", nrow(test_exp), "\n")

############### STEP 2：确保基因顺序一致 ################
common_genes <- intersect(rownames(train_exp), rownames(test_exp))
train_exp <- train_exp[common_genes, ]
test_exp  <- test_exp[common_genes, ]

############### STEP 3：整理 Training_expr（行=gene，列=sample） ################
# 对 3.0：必须行=基因，列=样本
Training_expr <- train_exp
write.table(Training_expr,
            file = "Training_expr.txt",
            sep = "\t", quote = FALSE, col.names = NA)

############### STEP 4：准备 Training_class ################
Training_class <- data.frame(
  outcome = ifelse(meta_59446$Group == "KBD", 1, 0)
)
rownames(Training_class) <- rownames(meta_59446)  # 样本名必须匹配
rownames(Training_class)=gsub("-",".",rownames(Training_class))
Training_class <- Training_class[colnames(Training_expr), , drop = FALSE]

write.table(Training_class,
            file = "Training_class.txt",
            sep = "\t", quote = FALSE, col.names = NA)

############### STEP 5：整理 Testing_expr ################
Testing_expr <- test_exp
write.table(Testing_expr,
            file = "Testing_expr.txt",
            sep = "\t", quote = FALSE, col.names = NA)

############### STEP 6：准备 Testing_class（必须含 Cohort 列） ################
Testing_class <- data.frame(
  outcome = ifelse(meta_32127$Group == "KBD", 1, 0),
  Cohort  = "GSE32127"
)
rownames(Testing_class) <- rownames(meta_32127)
Testing_class <- Testing_class[colnames(Testing_expr), , drop = FALSE]

write.table(Testing_class,
            file = "Testing_class.txt",
            sep = "\t", quote = FALSE, col.names = NA)

############### STEP 7：输出检查 ################
cat("==== 文件生成完成 ====\n")
cat("Training_expr.txt：", nrow(Training_expr), "genes ×", ncol(Training_expr), "samples\n")
cat("Training_class.txt：", nrow(Training_class), "samples\n")
cat("Testing_expr.txt：", nrow(Testing_expr), "genes ×", ncol(Testing_expr), "samples\n")
cat("Testing_class.txt：", nrow(Testing_class), "samples\n")
