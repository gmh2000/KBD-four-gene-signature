###############WGCNA#########
rm(list=ls())
load("GSE186593_dds.RData")
load("GSE186593_raw.RData")
vsd <- vst(dds, blind=TRUE)
vst_mat <- assay(vsd)
datExpr0  <- t(vst_mat) 
var_cut <- apply(datExpr0, 2, var)
keep_genes <- names(sort(var_cut, decreasing = TRUE))[1:12000]  # 例如选择 top 1w 基因
datExpr <- datExpr0[, keep_genes]

# 选择前75%变异性的基因
#var_threshold <- quantile(var_cut, probs = 0.75)
#keep_genes <- names(var_cut)[var_cut >= var_threshold]
#datExpr <- datExpr0[, keep_genes]

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
powers = c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf("SoftThreshold_GSE186593.pdf", width=8, height=6)
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     col="red")
abline(h=0.8, col="red")
dev.off()
###########构建邻接矩阵 + TOM##############
softPower <- 10 # 根据上一步调整
adjacency <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

###########层次聚类 + 模块识别###########
geneTree <- hclust(as.dist(dissTOM), method = "average")
minModuleSize <- 80
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit =1,
  pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
)
moduleColors <- labels2colors(dynamicMods)
length(unique(moduleColors))
table(moduleColors)

merge_threshold=0.1
# 自动合并模块
merged <- mergeCloseModules(datExpr, moduleColors, 
                            cutHeight = merge_threshold,
                            verbose = 3)
mergedColors <- merged$colors
mergedMEs <- merged$newMEs
length(unique(mergedColors))
table(mergedColors)
moduleColors=mergedColors


pdf("geneTree12000_80_1_0.1.pdf", width=8, height=6)
plotDendroAndColors(
  geneTree,
  moduleColors,
  groupLabels = "Module Colors",
  main = "Gene Dendrogram and Module Colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()
length(unique(moduleColors))
table(moduleColors)

##########计算模块特征基因（ME）并与表型关联########
meta=GSE186593_meta
datTraits <- data.frame(
  Control = as.numeric(meta$Group == "Normal"),
  KBD     = as.numeric(meta$Group == "KBD")
)
rownames(datTraits) <- rownames(datExpr0)   # datExpr0 = t(vst_mat)

# 模块特征值
MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs)

MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method="average")

pdf("MEs_dendrogram.pdf", width=8, height=6)
plot(METree, main="Clustering of Module Eigengenes", xlab="", sub="")
dev.off()
# 相关性矩阵 & P 值矩阵
moduleTraitCor     <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue  <- corPvalueStudent(moduleTraitCor, nrow(datExpr))


## 展开模块 × 性状
df <- expand.grid(
  Module = rownames(moduleTraitCor),
  Trait  = colnames(moduleTraitCor)
)

df$Correlation <- as.vector(moduleTraitCor)
df$Pvalue      <- as.vector(moduleTraitPvalue)

## P 值 → 星号
df$Signif <- cut(
  df$Pvalue,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)

## 标签：r + 星号（两行）
df$Label <- paste0(
  sprintf("r=%.2f", df$Correlation),
  "\n",
  df$Signif
)

## 模块按 |cor 与 KBD| 大小排序
order_modules <- df %>%
  filter(Trait == "KBD") %>%
  arrange(desc(abs(Correlation))) %>%
  pull(Module)

df$Module <- factor(df$Module, levels = order_modules)

## Trait 顺序：Control → KBD
df$Trait <- factor(df$Trait, levels = c("Control", "KBD"))
# 这里要得到 “每个模块的颜色”
module_levels <- levels(df$Module)         # MEturquoise, MEblue ...
module_cols   <- gsub("^ME", "", module_levels)  # 取真正的颜色名：turquoise, blue ...

# 转成 HEX 颜色
module_hex <- sapply(
  module_cols,
  function(col) {
    rgb(t(col2rgb(col))/255)
  }
)
names(module_hex) <- module_levels

# 给每一行贴上对应模块的 HEX 色
df$ColorHex <- module_hex[as.character(df$Module)]

pdf("Module_Trait_Heatmap_Publication.pdf", width=9, height=10)

ggplot(df, aes(x = Trait, y = Module)) +
  
  ## ---- 第一套 fill：相关性渐变色 ----
geom_tile(aes(fill = Correlation), color = "white", size = 0.4) +
  
  geom_text(aes(label = Label), size = 4.0, fontface = "bold") +
  
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  
  ## ---- 新的 fill 标度，用于模块颜色条 ----
ggnewscale::new_scale_fill() +
  
  geom_tile(
    aes(x = -0.15, y = Module, fill = ColorHex),
    width = 0.15,
    inherit.aes = FALSE,
    color = "white"
  ) +
  
  scale_fill_identity(guide = "none") +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 12, face = "bold"),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  
  ggtitle("Module–Trait Relationships (Control vs KBD)")

dev.off()
###############
pdf("Module_Trait_Heatmap_Publication.pdf", width=9, height=10)

ggplot(df, aes(x = Trait, y = Module)) +
  
  ## ---- 第一套 fill：相关性渐变色 ----
geom_tile(aes(fill = Correlation), color = "white", size = 0.4) +
  
  ## ---- 加大加粗数字标签 ----
geom_text(aes(label = Label), 
          size = 6.0,        # 从 4.0 增加到 5.0
          fontface = "bold") + 
  
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D81B60",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  
  ## ---- 新的 fill 标度，用于模块颜色条 ----
ggnewscale::new_scale_fill() +
  
  geom_tile(
    aes(x = -0.1, y = Module, fill = ColorHex),
    width = 0.15,
    inherit.aes = FALSE,
    color = "white"
  ) +
  
  scale_fill_identity(guide = "none") +
  
  theme_minimal(base_size = 17) +  # 从 14 增加到 16
  
  theme(
    ## ---- 加大加粗横坐标 ----
    axis.text.x = element_text(
      size = 16,        # 从 14 增加到 16
      face = "bold",
      color = "black",
      margin = margin(t = 10)  # 增加上边距
    ),
    
    ## ---- 加大加粗纵坐标 ----
    axis.text.y = element_text(
      size = 14,        # 从 12 增加到 14
      face = "bold",
      color = "black",
      margin = margin(r = 5)  # 增加右边距
    ),
    
    axis.title = element_blank(),
    panel.grid = element_blank(),
    
    ## ---- 标题也适当加大 ----
    plot.title = element_text(
      size = 19,        # 从 16 增加到 18
      face = "bold", 
      hjust = 0.5,
      margin = margin(b = 15)  # 增加下边距
    )
  ) +
  
  ggtitle("Module–Trait Relationships")

dev.off()
################提取与 KBD 最显著相关的模块#########
module_of_interest <- rownames(moduleTraitCor)[which.max(abs(moduleTraitCor[, "KBD"]))]
module_of_interest
module_color <- gsub("^ME", "", module_of_interest)
module_color
genes_in_module <- colnames(datExpr)[moduleColors == module_color]
length(genes_in_module)
unique(moduleColors)
rownames(moduleTraitCor)
ncol(datExpr) == length(moduleColors)  
#################计算 Module Membership (MM) 与 Gene Significance (GS)#########
geneModuleMembership <- cor(datExpr, MEs, use="p")
geneTraitSignificance <- cor(datExpr, datTraits$KBD, use="p")
colnames(geneModuleMembership) <- substring(colnames(geneModuleMembership), 3) 
module_color <- substring(module_of_interest, 3)
genes_in_module <- colnames(datExpr)[moduleColors == module_color]

MM <- geneModuleMembership[ genes_in_module , module_color ]
GS <- geneTraitSignificance[ genes_in_module ,]
summary(datTraits$KBD)
summary(geneTraitSignificance)
cor(MM, GS)

module_genes_ranked <- data.frame(
  gene = genes_in_module,
  MM   = MM,
  GS   = GS
)

head(module_genes_ranked[order(-abs(module_genes_ranked$GS)), ], 20)
hub_genes <- subset(module_genes_ranked, abs(MM) > 0.9 & abs(GS) > 0.4)
save.image(file="WGCNA_12000")
save(hub_genes,genes_in_module,file="WGCNA_hub_genes.RData")
dim(deg_all)
overlap=intersect(deg_all$Gene,hub_genes$gene)
