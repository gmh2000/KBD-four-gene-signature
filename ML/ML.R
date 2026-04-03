# 设置系统编码
Sys.setlocale("LC_ALL", "Chinese")
options(encoding = "UTF-8")
#项目描述：12种机器学习算法进行特征筛选和预测模型构建
#算法包括：Lasso、Ridge、Enet、Stepglm、SVM、glmBoost、LDA、plsRglm、RandomForest、GBM、XGBoost、NaiveBayes
#主要流程：
#在交叉验证的框架下使用一种算法进行特征筛选，再使用另一种算法构建预测模型，
#最后使用得到的模型在外部数据集(除训练集外)计算ROC曲线下面积(AUC)，并通过热图将模型的性能可视化
# 算法输入：
# 训练集表达谱：行为样本(如SYMBOL/ENSEMBL等，需要与验证集一致)，列为特征的文本文件，格式见InputData文件夹中的Train_expr.txt
# 训练集临床信息：行为样本，包含需要预测的二元结局变量(支持[0,1]格式，请指定对应的列名)，格式见InputData文件夹中的Train_class.txt
# 验证集表达谱：行为样本(如SYMBOL/ENSEMBL等，需要与训练集一致)，列为特征的文本文件，格式见InputData文件夹中的Test_expr.txt
# 验证集临床信息：行为样本，包含需要预测的二元结局变量(支持[0,1]格式)以及一个队列指示信息的列名，格式见InputData文件夹中的Test_class.txt
# 注意：如果特征数量未经预筛选，包含数十或数千个感兴趣的特征，模型筛选和构建时间将相应延长

# 主要优势：
# 基于多种算法灵活进行特征筛选，并可视化针对不同队列的预测性能
# 可根据输出结果，选择最稳健的特征及对应的模型
# ...

rm(list=ls())
## 1 配置项目路径 -------------------------------------------------------
setwd("")
work.path <- getwd()
code.path <- file.path(work.path, "Codes")
data.path <- file.path(work.path, "InputData")
res.path  <- file.path(work.path, "Results")
fig.path  <- file.path(work.path, "Figures")

if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)

## 2 加载 R 包 ------------------------------------------------------------
library(randomForest)
library(survival)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
library(openxlsx)

# 加载模型训练以及模型评估的脚本
source(file.path(code.path, "ML.R"))
# 选择最终生成的模型类型：panML表示生成不同算法对应的模型，multiLogistic表示提取特征后构建多因素logistic模型
FinalModel <- c("panML", "multiLogistic")[2]

## 训练队列数据读取 ---------------------------------------------------------
# 训练集表达谱：行为特征(感兴趣的基因集合)，列为样本的表达矩阵(表达谱需要与验证集一致类型，同为SYMBOL或ENSEMBL等)
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"),
                         header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE,
                         stringsAsFactors = FALSE)
# 列为样本，包含需要预测的二元结局变量(支持[0,1]格式)
Train_class <- read.table(file.path(data.path, "Training_class.txt"),
                          header = TRUE, sep = "\t",
                          row.names = 1, check.names = FALSE,
                          stringsAsFactors = FALSE)
# 读取训练集的共同样本
comsam <- intersect(rownames(Train_class), colnames(Train_expr))
Train_expr  <- Train_expr[, comsam]
Train_class <- Train_class[comsam, , drop = FALSE]

## 验证队列数据读取 -------------------------------------------------------
# 验证集表达谱：行为特征(感兴趣的基因集合)，列为样本的表达矩阵(表达谱需要与训练集一致类型，同为SYMBOL或ENSEMBL等)
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"),
                        header = TRUE, sep = "\t",
                        row.names = 1, check.names = FALSE,
                        stringsAsFactors = FALSE)
# 列为样本，包含需要预测的二元结局变量(支持[0,1]格式)以及一个队列指示信息的列名
Test_class <- read.delim(file.path(data.path, "Testing_class.txt"),
                         header = TRUE, row.names = 1,
                         check.names = FALSE, stringsAsFactors = FALSE,
                         fill = TRUE)

# 读取验证集的共同样本
comsam <- intersect(rownames(Test_class), colnames(Test_expr))
Test_expr  <- Test_expr[, comsam]
Test_class <- Test_class[comsam, , drop = FALSE]

# 读取共同基因
comgene    <- intersect(rownames(Train_expr), rownames(Test_expr))
Train_expr <- Train_expr[comgene, ]
Test_expr  <- Test_expr[comgene, ]
colnames(Train_expr) <- make.names(colnames(Train_expr))
colnames(Test_expr)  <- make.names(colnames(Test_expr))
# 转置矩阵，使行=样本，列=特征
Train_expr <- t(Train_expr)
Test_expr  <- t(Test_expr)

## 数据标准化 -------------------------------------------------------------
scaleData <- function(data, cohort = NULL, centerFlags = TRUE, scaleFlags = TRUE) {
  # 若未指定cohort，则视为单一队列
  if (is.null(cohort)) {
    data <- scale(data, center = centerFlags, scale = scaleFlags)
    return(as.data.frame(data))
  } else {
    cohort <- as.factor(cohort)
    data.sc <- lapply(levels(cohort), function(co) {
      idx <- which(cohort == co)
      tmp <- data[idx, , drop = FALSE]
      tmp <- scale(tmp, center = centerFlags, scale = scaleFlags)
      return(tmp)
    })
    data.sc <- do.call(rbind, data.sc)
    data.sc <- data.sc[match(rownames(data), rownames(data.sc)), ]
    return(as.data.frame(data.sc))
  }
}

## 根据cohort信息对训练集和验证集数据进行标准化 ----------------------------
Train_set <- scaleData(data = Train_expr,
                       centerFlags = TRUE, scaleFlags = TRUE)

Test_set <- scaleData(data = Test_expr,
                      cohort = Test_class$Cohort,
                      centerFlags = TRUE, scaleFlags = TRUE)
Train_set <- data.matrix(Train_set)
Test_set  <- data.matrix(Test_set)
# 保存备份
Train_set_bk <- Train_set

## 机器学习算法配置 --------------------------------------------------------
methods <- read.xlsx(file.path(code.path, "methods.xlsx"), sheet = 1)

## 机器学习算法配置 --------------------------------------------------------
methods <- c(
  "Lasso + Stepglm [both]",
  "SVM",
  "glmBoost + SVM",
  "Ridge",
  "Lasso + SVM",
  "glmBoost + Ridge",
  "Enet [alpha=0.1]",
  "glmBoost + Enet [alpha=0.1]",
  "Enet [alpha=0.2]",
  "Enet [alpha=0.3]",
  "glmBoost + Enet [alpha=0.3]",
  "glmBoost + Enet [alpha=0.2]",
  "Enet [alpha=0.4]",
  "glmBoost + Enet [alpha=0.4]",
  "Lasso + glmBoost",
  "Enet [alpha=0.5]",
  "glmBoost",
  "glmBoost + Enet [alpha=0.5]",
  "Enet [alpha=0.6]",
  "glmBoost + Enet [alpha=0.6]",
  "glmBoost + Enet [alpha=0.7]",
  "glmBoost + Enet [alpha=0.8]",
  "Enet [alpha=0.8]",
  "Enet [alpha=0.9]",
  "Lasso",
  "Enet [alpha=0.7]",
  "glmBoost + Enet [alpha=0.9]",
  "glmBoost + Lasso",
  "Lasso + plsRglm",
  "glmBoost + plsRglm",
  "glmBoost + Stepglm [forward]",
  "Lasso + Stepglm [forward]",
  "RF + SVM",
  "Stepglm [forward]",
  "plsRglm",
  "RF + Ridge",
  "RF + Enet [alpha=0.1]",
  "RF + plsRglm",
  "RF + Stepglm [forward]",
  "RF + Enet [alpha=0.2]",
  "RF + Enet [alpha=0.3]",
  "RF + Enet [alpha=0.6]",
  "RF + Lasso",
  "RF + Enet [alpha=0.7]",
  "RF + Enet [alpha=0.5]",
  "RF + glmBoost",
  "RF + Enet [alpha=0.9]",
  "RF + Enet [alpha=0.4]",
  "RF + Enet [alpha=0.8]",
  "RF + Stepglm [both]",
  "RF + Stepglm [backward]",
  "Stepglm [both] + Ridge",
  "Stepglm [backward] + Ridge",
  "Stepglm [both] + plsRglm",
  "Stepglm [backward] + plsRglm",
  "Stepglm [both] + Enet [alpha=0.9]",
  "Stepglm [backward] + Enet [alpha=0.9]",
  "Stepglm [both] + Enet [alpha=0.1]",
  "Stepglm [backward] + Enet [alpha=0.1]",
  "Stepglm [both] + Enet [alpha=0.8]",
  "Stepglm [backward] + Enet [alpha=0.8]",
  "Stepglm [both] + Enet [alpha=0.2]",
  "Stepglm [backward] + Enet [alpha=0.2]",
  "Stepglm [both] + Lasso",
  "Stepglm [backward] + Lasso",
  "Stepglm [both] + Enet [alpha=0.6]",
  "Stepglm [backward] + Enet [alpha=0.6]",
  "glmBoost + GBM",
  "Stepglm [both] + Enet [alpha=0.7]",
  "Stepglm [backward] + Enet [alpha=0.7]",
  "Lasso + Stepglm [backward]",
  "Stepglm [both]",
  "Stepglm [backward]",
  "glmBoost + Stepglm [both]",
  "glmBoost + Stepglm [backward]",
  "Stepglm [both] + Enet [alpha=0.4]",
  "Stepglm [backward] + Enet [alpha=0.4]",
  "Stepglm [both] + Enet [alpha=0.3]",
  "Stepglm [backward] + Enet [alpha=0.3]",
  "Stepglm [both] + glmBoost",
  "Stepglm [backward] + glmBoost",
  "Stepglm [both] + Enet [alpha=0.5]",
  "Stepglm [backward] + Enet [alpha=0.5]",
  "glmBoost + RF",
  "RF",
  "Lasso + GBM",
  "RF + GBM",
  "GBM",
  "Stepglm [both] + SVM",
  "Stepglm [backward] + SVM",
  "Lasso + RF",
  "Stepglm [both] + GBM",
  "Stepglm [backward] + GBM",
  "Stepglm [both] + RF",
  "LDA",
  "glmBoost + LDA",
  "RF+LDA",
  "Stepglm [both]+LDA",
  "Stepglm [backward]+LDA",
  "Lasso + LDA",
  "Stepglm [backward] + RF",
  "XGBoost",
  "Lasso+XGBoost",
  "glmBoost+XGBoost",
  "RF+XGBoost",
  "Stepglm[both]+XGBoost",
  "Stepglm[backward]+XGBoost",
  "NaiveBayes",
  "Lasso+NaiveBayes",
  "glmBoost+NaiveBayes",
  "RF+NaiveBayes",
  "Stepglm[both]+NaiveBayes",
  "Stepglm[backward]+NaiveBayes"
)

methods <- methods[methods != ""]
classVar          <- "outcome"  # 指定需要预测的列名(支持[0,1]二元分类格式)
min.selected.var  <- 2        # 构建模型最少需要的特征数

## 预训练 - 特征筛选 --------------------------------------------------------
Variable         <- colnames(Train_set)
preTrain.method  <- strsplit(methods, "\\+")
preTrain.method  <- lapply(preTrain.method, function(x) rev(x)[-1])  # 去掉最后一个（建模算法）
preTrain.method  <- unique(unlist(preTrain.method))
preTrain.var <- list()

set.seed(777)
for (method in preTrain.method) {
  preTrain.var[[method]] <- RunML(
    method     = method,
    Train_set  = Train_set,
    Train_label = Train_class,
    mode       = "Variable",
    classVar   = classVar
  )
}

## 组合筛选特征并训练模型 ----------------------------------------------------
model   <- list()
fea_list <- list()

for (method in methods) {
  method_parts <- strsplit(method, "\\+")[[1]]
  fea_method   <- method_parts[-length(method_parts)]
  fit_method   <- tail(method_parts, 1)
  
  cat("当前方法：", method, "\n")
  
  feats <- unique(unlist(preTrain.var[fea_method]))
  max.feat <- 20
  if (length(feats) > max.feat) {
    feats <- feats[1:max.feat]   # 或者根据你前面排序的结果取前 N 个
  }
  if (length(feats) < min.selected.var) {
    cat("  → 特征数不足，跳过该方法\n\n")
    next
  }
  
  fea_list[[method]] <- feats
  
  tmpTrain <- Train_set[, feats, drop = FALSE]
  
  tryCatch({
    model[[method]] <- RunML(
      method      = fit_method,
      Train_set   = tmpTrain,
      Train_label = Train_class,
      mode        = "Model",
      classVar    = classVar
    )
  }, error = function(e) {
    cat("  → 构建模型失败:", e$message, "\n")
    model[[method]] <- NULL
  })
  
  cat("\n")
}

cat("模型训练完成。成功构建",
    sum(sapply(model, function(x) !is.null(x))),
    "个模型\n")

Train_set <- Train_set_bk
rm(Train_set_bk)
saveRDS(model, file.path(res.path, "model.rds"))

if (FinalModel == "multiLogistic") {
  # 只保留成功训练的模型
  non_null_models <- Filter(Negate(is.null), model)
  
  logisticmodel <- lapply(non_null_models, function(fit) {
    # 提取当前模型使用的特征
    feats <- ExtractVar(fit)
    
    # 如果提取不到特征，直接返回 NULL
    if (is.null(feats) || length(feats) == 0) {
      return(NULL)
    }
    
    # 只保留在 Train_set 里真正存在的特征，防止“下标出界”
    valid_feats <- intersect(feats, colnames(Train_set))
    
    # 如果一个合法特征都没有，就跳过这个模型
    if (length(valid_feats) == 0) {
      return(NULL)
    }
    
    # 组装一个包含 outcome + 特征的数据框
    dat <- data.frame(
      Train_class[, classVar, drop = FALSE],
      Train_set[, valid_feats, drop = FALSE]
    )
    colnames(dat)[1] <- classVar  # 第一列列名是 classVar，例如 "outcome"
    
    # 构造公式：outcome ~ .
    form <- as.formula(paste(classVar, "~ ."))
    
    tmp <- glm(
      formula = form,
      family  = binomial(),
      data    = dat
    )
    
    # 记录这个 logistic 模型实际用到的特征
    tmp$subFeature <- valid_feats
    return(tmp)
  })
  
  # 去掉 lapply 中返回的 NULL 项
  logisticmodel <- Filter(Negate(is.null), logisticmodel)
  
}

saveRDS(logisticmodel, file.path(res.path, "logisticmodel.rds"))

## 训练集内部近似 10 折交叉验证 --------------------------------------
## 使用最终筛选出的特征集，在 GSE59446 内部评估各算法模型的稳定性

methodsValid <- names(model)

set.seed(2025)
folds <- createFolds(Train_class[[classVar]], k = 10, list = TRUE)

CV_AUC_mat <- matrix(
  NA,
  nrow = length(methodsValid),
  ncol = length(folds),
  dimnames = list(
    methodsValid,
    paste0("Fold", seq_along(folds))
  )
)

for (method in methodsValid) {
  cat("进行训练集 10 折 CV：", method, "\n")
  
  # 1) 该方法对应的模型如果是 NULL，直接跳过
  fit0 <- model[[method]]
  if (is.null(fit0)) {
    cat("  → 模型为 NULL，跳过该方法\n")
    next
  }
  
  # 2) 提取特征，并确保都在 Train_set 里
  feats <- ExtractVar(fit0)
  if (is.null(feats) || length(feats) == 0) {
    cat("  → 未提取到特征，跳过该方法\n")
    next
  }
  feats <- intersect(feats, colnames(Train_set))
  if (length(feats) == 0) {
    cat("  → 提取的特征在 Train_set 中均不存在，跳过该方法\n")
    next
  }
  
  # 3) 解析出用于建模的算法名称（如 Lasso+SVM → SVM）
  method_parts <- strsplit(method, "\\+")[[1]]
  clf_method   <- tail(method_parts, 1)
  
  for (i in seq_along(folds)) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(Train_set)), test_idx)
    
    x_train <- Train_set[train_idx, feats, drop = FALSE]
    y_train <- Train_class[train_idx, , drop = FALSE]
    
    x_test  <- Train_set[test_idx, feats, drop = FALSE]
    y_test  <- Train_class[test_idx, , drop = FALSE]
    
    # 如果这一折里仍然出现 0 个变量（极端情况），也跳过这一折
    if (ncol(x_train) == 0 || ncol(x_test) == 0) {
      cat("    Fold", i, ": 变量数为 0，跳过本折\n")
      next
    }
    
    # 4) 当前折上重训练模型（仅使用建模算法，不再重复特征筛选）
    fit_cv <- RunML(
      method      = clf_method,
      Train_set   = x_train,
      Train_label = y_train,
      mode        = "Model",
      classVar    = classVar
    )
    
    RS_cv <- CalPredictScore(fit = fit_cv, new_data = x_test)
    
    roc_obj <- roc(
      response  = y_test[[classVar]],
      predictor = as.vector(RS_cv),
      quiet     = TRUE
    )
    CV_AUC_mat[method, i] <- as.numeric(auc(roc_obj))
  }
}

# 计算各方法的平均 AUC（忽略 NA）
CV_AUC_mean <- rowMeans(CV_AUC_mat, na.rm = TRUE)

write.table(
  CV_AUC_mat,
  file.path(res.path, "Train_CV_AUC_by_fold.txt"),
  sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE
)

write.table(
  CV_AUC_mean,
  file.path(res.path, "Train_CV_AUC_mean.txt"),
  sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
)


## 评估模型 -----------------------------------------------------

# 读取已保存的模型列表
model <- readRDS(file.path(res.path, "model.rds"))
# model <- readRDS(file.path(res.path, "logisticmodel.rds"))
methodsValid <- names(model)

# 根据各算法模型计算风险得分（训练集 + 验证集）
# 根据各算法模型计算风险得分（训练集 + 验证集）
RS_list <- list()
for (method in methodsValid) {
  # 取该方法对应的特征集合
  feats <- fea_list[[method]]
  
  # 如果没特征或特征不在 Train_set 里，直接跳过
  if (is.null(feats) || length(feats) == 0) {
    cat("  →", method, "无特征，跳过 RS 计算\n")
    next
  }
  valid_feats <- intersect(feats, colnames(Train_set))
  if (length(valid_feats) == 0) {
    cat("  →", method, "特征在 Train_set 中不存在，跳过 RS 计算\n")
    next
  }
  
  # 训练集 + 验证集合并，然后只保留该模型实际使用的特征列
  all_data <- rbind(Train_set, Test_set)
  all_data <- all_data[, valid_feats, drop = FALSE]
  
  RS_list[[method]] <- CalPredictScore(
    fit      = model[[method]],
    new_data = all_data
  )
}

RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat,
            file.path(res.path, "RS_mat.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)


# 根据各算法模型预测分类（训练集 + 验证集）
# 根据各算法模型预测分类（训练集 + 验证集）
Class_list <- list()
for (method in methodsValid) {
  feats <- fea_list[[method]]
  
  if (is.null(feats) || length(feats) == 0) {
    cat("  →", method, "无特征，跳过分类预测\n")
    next
  }
  valid_feats <- intersect(feats, colnames(Train_set))
  if (length(valid_feats) == 0) {
    cat("  →", method, "特征在 Train_set 中不存在，跳过分类预测\n")
    next
  }
  
  all_data <- rbind(Train_set, Test_set)
  all_data <- all_data[, valid_feats, drop = FALSE]
  
  Class_list[[method]] <- PredictClass(
    fit      = model[[method]],
    new_data = all_data
  )
}

Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
write.table(
  Class_mat,
  file.path(res.path, "Class_mat.txt"),
  sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE
)

# 计算 AUC：第一列为训练集 GSE59446，后面为各验证队列
Train_name <- "GSE59446"   # 训练集名称，将作为 AUC 矩阵中的第一列列名

AUC_list <- list()
for (method in methodsValid) {
  AUC_list[[method]] <- RunEval(
    fit        = model[[method]],
    Test_set   = Test_set,      # 验证集表达矩阵（多个 cartilage cohort 合并）
    Test_label = Test_class,    # 验证集临床信息（含 Cohort、outcome）
    Train_set   = Train_set,    # 训练集表达矩阵（GSE59446）
    Train_label = Train_class,  # 训练集临床信息
    Train_name  = "GSE59446",   # 训练集标签，用于区分 cohort
    cohortVar   = "Cohort",
    classVar    = classVar
  )
}
AUC_mat <- do.call(rbind, AUC_list)
# 先按原来方式算验证集 AUC
# 1. 只保留在 AUC_mat 和 CV_AUC_mean 里“共同存在”的方法
common_methods <- intersect(rownames(AUC_mat), names(CV_AUC_mean))

AUC_mat <- AUC_mat[common_methods, , drop = FALSE]

# 2. 用 CV 的平均 AUC 覆盖训练集这一列（GSE59446）
AUC_mat[, "GSE59446"] <- CV_AUC_mean[common_methods]

# 3. 确保没有 NA 或空行名
AUC_mat <- AUC_mat[!is.na(rownames(AUC_mat)) & rownames(AUC_mat) != "", , drop = FALSE]
write.table(
  AUC_mat,
  file.path(res.path, "AUC_mat.txt"),
  sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE
)

# 绘图部分 --------------------------------------------------------------------

AUC_mat <- read.table(
  file.path(res.path, "AUC_mat.txt"),
  sep = "\t", row.names = 1,
  header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE
)

# 保证训练集列（GSE59446）在最左侧，其余为各验证队列
Train_name <- "GSE59446"
train_col  <- Train_name
test_cols  <- setdiff(colnames(AUC_mat), train_col)
AUC_mat    <- AUC_mat[, c(train_col, test_cols), drop = FALSE]
AUC_mat <- AUC_mat[apply(!is.na(AUC_mat), 1, all), , drop = FALSE]
#################
if (ncol(AUC_mat) < 2) {
  stop("AUC_mat 列数 < 2，无法只用前两列计算平均 AUC，请检查。")
}

blood_cols <- colnames(AUC_mat)[1:2]  # 第 1 列：GSE59446；第 2 列：第一个血液验证队列

avg_AUC    <- rowMeans(AUC_mat[, blood_cols, drop = FALSE], na.rm = TRUE)

# 按平均 AUC 从高到低排序

avg_AUC <- sort(avg_AUC, decreasing = TRUE)
AUC_mat <- AUC_mat[names(avg_AUC), , drop = FALSE]

fea_sel <- fea_list[[rownames(AUC_mat)[1]]]   # 保留你原来的 fea_sel 逻辑

# 数字标签（保留 3 位小数）
avg_AUC_num   <- as.numeric(avg_AUC)
avg_AUC_label <- sprintf("%.3f", avg_AUC_num)





# 设置不同队列的颜色
if (ncol(AUC_mat) < 3) {
  CohortCol <- c("red", "blue")
} else {
  CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")
}
colnames(AUC_mat)=c("Training","Validation(Blood)","Validation(Cartilage)")
names(CohortCol) <- colnames(AUC_mat)

cellwidth  <- 1
cellheight <- 0.5
hm <- SimpleHeatmap(
  AUC_mat,     # 热图数据
  avg_AUC_num,     # 右侧条形图（平均 AUC）
  CohortCol,   # 列标签颜色
  "steelblue", # 右侧条形图颜色
  cellwidth  = cellwidth,
  cellheight = cellheight,
  cluster_columns = FALSE,
  cluster_rows    = FALSE
)

pdf(
  file.path(fig.path, "AUC.pdf"),
  width  = cellwidth * ncol(AUC_mat) + 7,
  height = cellheight * nrow(AUC_mat) * 0.45
)
draw(hm)
invisible(dev.off())

print("中文测试")

best_method <- names(avg_AUC)[1]
cat("最佳方法：", best_method, "\n")
cat("对应平均 AUC：", round(avg_AUC[1], 3), "\n")

# 取该方法对应的特征列表
best_genes <- fea_list[["Stepglm [backward]+LDA"]]
best_genes <- intersect(best_genes, colnames(Train_set))  # 双重保险，防止不存在的列
cat("该方法最终使用的基因数：", length(best_genes), "\n")
print(best_genes)

##======================================================
## 依据 AUC_mat 选出最佳方法，并画 ROC 图
##======================================================

# 确保以下对象在环境中：AUC_mat, model, fea_list, Train_set, Test_set,
#                        Train_class, Test_class, classVar, fig.path

library(pROC)

## 1) 选出平均 AUC 最高的方法
avg_AUC <- apply(AUC_mat, 1, mean, na.rm = TRUE)
avg_AUC <- sort(avg_AUC, decreasing = TRUE)

best_method <- names(avg_AUC)[1]
cat("最佳方法：", best_method, "\n")
cat("平均 AUC：", round(avg_AUC[1], 3), "\n")

best_genes <- fea_list[[best_method]]
best_genes <- intersect(best_genes, colnames(Train_set))
cat("该方法使用的基因：", paste(best_genes, collapse = ", "), "\n\n")

## 2) 组装表达矩阵 & 表型信息
# 表达矩阵：样本 × 基因
x_train <- Train_set[, best_genes, drop = FALSE]
x_test  <- Test_set[,  best_genes, drop = FALSE]
all_x   <- rbind(x_train, x_test)

# 表型：outcome + Cohort
train_anno <- data.frame(
  sample  = rownames(Train_set),
  outcome = Train_class[[classVar]],
  Cohort  = "GSE59446",
  stringsAsFactors = FALSE
)

test_anno <- data.frame(
  sample  = rownames(Test_set),
  outcome = Test_class[[classVar]],
  Cohort  = Test_class$Cohort,
  stringsAsFactors = FALSE
)

all_anno <- rbind(train_anno, test_anno)
rownames(all_anno) <- all_anno$sample

# 对齐顺序
all_x <- all_x[all_anno$sample, , drop = FALSE]

## 3) 用最佳 pan-ML 模型对所有样本打分
scores_all <- CalPredictScore(
  fit      = model[[best_method]],
  new_data = all_x
)

all_anno$score <- as.numeric(scores_all[all_anno$sample])

## 4) 逐队列画 ROC（叠加在一张图上）
## 4) 三个队列分别单独画 ROC ---------------------------------------
library(pROC)

cohorts <- unique(all_anno$Cohort)
cohorts
# 比如: "GSE59446" "GSE32127" "E-MEXP-3474"  (如果还有别的，就按实际为准)

for (co in cohorts) {
  idx <- which(all_anno$Cohort == co)
  
  roc_obj <- roc(
    response  = all_anno$outcome[idx],
    predictor = all_anno$score[idx],
    quiet     = TRUE
  )
  
  # 生成当前 cohort 的单独 ROC 图文件
  fn <- paste0("ROC_", co, "_", gsub("[^A-Za-z0-9]+", "_", best_method), ".pdf")
  pdf(file.path(fig.path, fn), width = 5, height = 5)
  
  plot(
    roc_obj,
    col   = "red",
    lwd   = 2,
    main  = paste0(co, "\nAUC = ", round(auc(roc_obj), 3)),
    xlab  = "1 - Specificity",
    ylab  = "Sensitivity"
  )
  abline(0, 1, lty = 2, col = "grey70")
  
  dev.off()
  cat("已输出：", fn, "\n")
}

###########################################

library(pROC)
library(RColorBrewer)

## 1) 选定 5 个基因 ------------------------------------------------
selected_genes <- c("C4B", "AQP1", "HBA2", "ACSL6")

# 只保留在表达矩阵中确实存在的基因（双保险）
selected_genes <- intersect(selected_genes, colnames(Train_set))
selected_genes <- intersect(selected_genes, colnames(Test_set))
cat("本次绘图使用的基因：", paste(selected_genes, collapse = ", "), "\n")

## 2) 整理表达矩阵 & 表型信息 --------------------------------------

# 合并训练集和验证集表达，行=样本，列=基因
all_expr <- rbind(Train_set, Test_set)
all_expr <- all_expr[, selected_genes, drop = FALSE]

# 训练集注释
train_anno <- data.frame(
  sample  = rownames(Train_set),
  outcome = Train_class[[classVar]],
  Cohort  = "GSE59446",
  stringsAsFactors = FALSE
)

# 验证集注释（包含 Cohort）
test_anno <- data.frame(
  sample  = rownames(Test_set),
  outcome = Test_class[[classVar]],
  Cohort  = Test_class$Cohort,
  stringsAsFactors = FALSE
)

# 合并
all_anno <- rbind(train_anno, test_anno)
rownames(all_anno) <- all_anno$sample

# 对齐顺序
all_expr <- all_expr[all_anno$sample, , drop = FALSE]

## 3) 决定要画的 3 个数据集 ---------------------------------------

# 如果你就画“训练集 + 两个软骨”，手动指定（按你真实的 Cohort 名改）：
cohorts_to_plot <- c("GSE59446", "GSE32127", "E-MEXP-3474")
# 如果名字不完全一致，用 unique(all_anno$Cohort) 看一眼，然后改成对应名称

cat("将为以下 Cohort 绘制 5 基因 ROC 图：\n")
print(cohorts_to_plot)

## 4) 为每个 Cohort 单独画一张“5 曲线 ROC 图” -----------------------

# 5 个基因的配色
cols_gene <- brewer.pal(5, "Dark2")
names(cols_gene) <- selected_genes

for (co in cohorts_to_plot) {
  idx <- which(all_anno$Cohort == co)
  
  if (length(idx) == 0) {
    cat("Cohort", co, "在 all_anno 中不存在，跳过。\n")
    next
  }
  
  y_co <- all_anno$outcome[idx]
  X_co <- all_expr[idx, , drop = FALSE]  # 只 5 个基因
  
  # 防止当前队列只有一种结局
  if (length(unique(y_co)) < 2) {
    cat("Cohort", co, "只有单一结局类别，无法计算 ROC，跳过。\n")
    next
  }
  
  # 打开 PDF 文件
  fn <- paste0("ROC_5genes_", co, ".pdf")
  pdf(file.path(fig.path, fn), width = 5, height = 5)
  
  # 空图
  plot(0:1, 0:1, type = "n",
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity",
       ylab = "Sensitivity",
       main = paste0("ROC of 5 genes - ", co))
  abline(0, 1, lty = 2, col = "grey80")
  
  legend_text <- character(0)
  legend_col  <- character(0)
  
  # 逐基因画 ROC
  for (gene in selected_genes) {
    if (!gene %in% colnames(X_co)) next
    
    expr_g <- X_co[, gene]
    # 如果该基因在该队列中几乎没变化，跳过
    if (length(unique(expr_g)) < 2) {
      next
    }
    
    roc_obj <- roc(
      response  = y_co,
      predictor = as.numeric(expr_g),
      quiet     = TRUE
    )
    
    lines(1 - roc_obj$specificities,
          roc_obj$sensitivities,
          col = cols_gene[gene],
          lwd = 2)
    
    legend_text <- c(
      legend_text,
      paste0(gene, " (AUC = ", round(auc(roc_obj), 3), ")")
    )
    legend_col  <- c(legend_col, cols_gene[gene])
  }
  
  if (length(legend_text) > 0) {
    legend("bottomright", legend = legend_text,
           col = legend_col, lwd = 2, cex = 0.7)
  } else {
    mtext("No valid ROC curves for these genes.", side = 3, col = "red")
  }
  
  dev.off()
  cat("已输出：", fn, "\n")
}


##======================================================
## 5 个基因在三套数据集中 case/control 表达箱线图
##======================================================

library(ggplot2)
library(reshape2)
library(ggplot2)
library(reshape2)
library(ggpubr)

# 1) 选定 5 个基因 ------------------------------------------------
selected_genes <- c("C4B", "AQP1", "HBA2", "ACSL6")

# 只保留在表达矩阵中存在的基因（双保险）
selected_genes <- intersect(selected_genes, colnames(Train_set))
selected_genes <- intersect(selected_genes, colnames(Test_set))
cat("本次箱线图使用的基因：", paste(selected_genes, collapse = ", "), "\n")

# 如果有缺失，直接提示
if (length(selected_genes) == 0) {
  stop("选定的基因在表达矩阵中都不存在，请检查基因名是否和列名一致。")
}

## 2) 合并表达与表型信息 ------------------------------------------

# 合并训练集与验证集表达矩阵：行 = 样本，列 = 基因
all_expr <- rbind(Train_set, Test_set)
all_expr <- all_expr[, selected_genes, drop = FALSE]

# 训练集注释
train_anno <- data.frame(
  sample  = rownames(Train_set),
  outcome = Train_class[[classVar]],
  Cohort  = "GSE59446",
  stringsAsFactors = FALSE
)

# 验证集注释（包含 Cohort）
test_anno <- data.frame(
  sample  = rownames(Test_set),
  outcome = Test_class[[classVar]],
  Cohort  = Test_class$Cohort,
  stringsAsFactors = FALSE
)

# 合并
all_anno <- rbind(train_anno, test_anno)
rownames(all_anno) <- all_anno$sample

# 对齐顺序
all_expr <- all_expr[all_anno$sample, , drop = FALSE]

# 整理 outcome 为因子，并给出清晰标签
# 假设 outcome: 0 = Control, 1 = Case（KBD），可按需要调整
all_anno$group <- factor(
  all_anno$outcome,
  levels = c(0, 1),
  labels = c("Control", "Case")
)

## 3) 转为长表格式 -------------------------------------------------

df_long <- cbind(
  sample = rownames(all_expr),
  all_anno[, c("Cohort", "group")],
  as.data.frame(all_expr, check.names = FALSE)
)

df_melt <- melt(
  df_long,
  id.vars = c("sample", "Cohort", "group"),
  variable.name = "Gene",
  value.name = "Expression"
)

# 确保基因顺序固定（按你提供的顺序）
df_melt$Gene <- factor(df_melt$Gene, levels = selected_genes)
write.csv(df_melt,file="df_melt.csv",row.names = F)
getwd()
## 4) 画图：三套数据集分别显示，X 轴基因，按 group 分箱 -----------
####################################################
cohort_mapping <- c(
  "GSE59446" = "Training set",
  "GSE32127" = "Validation set (Blood)",
  "E-MEXP-3474" = "Validation set (Cartilage)"
)

# 应用重命名
df_melt$Cohort_renamed <- cohort_mapping[df_melt$Cohort]

# 定义数据集的显示顺序（按您指定的顺序）
cohort_order <- c(
  "Training set",
  "Validation set (Blood)", 
  "Validation set (Cartilage)"
)

# 将Cohort转换为因子，控制显示顺序
df_melt$Cohort_renamed <- factor(
  df_melt$Cohort_renamed,
  levels = cohort_order
)

# 检查哪些数据集实际存在
cat("数据集信息：\n")
cat("原始名称 -> 新名称\n")
for (orig in names(cohort_mapping)) {
  if (orig %in% df_melt$Cohort) {
    cat(sprintf("%s -> %s\n", orig, cohort_mapping[orig]))
  }
}
cat("\n实际使用的数据集：", paste(unique(df_melt$Cohort_renamed), collapse = ", "), "\n")

## 5) 自定义颜色方案 ------------------------------------------------
group_colors <- c("Control" = "#2E86AB",  # 深蓝色
                  "Case" = "#F24236")     # 珊瑚红

## 6) 统计检验并获取精确p值 ----------------------------------------
stats_df <- df_melt %>%
  group_by(Cohort, Gene) %>%
  summarise(
    p_value = wilcox.test(Expression ~ group)$p.value,
    y_max = max(Expression) * 1.15,
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    p_label = ifelse(p_value < 0.001, 
                     "p < 0.001", 
                     sprintf("p = %.3f", p_value))
  )

# 为统计结果也应用重命名
stats_df$Cohort_renamed <- cohort_mapping[stats_df$Cohort]
stats_df$Cohort_renamed <- factor(stats_df$Cohort_renamed, levels = cohort_order)

## 7) 创建美化的箱线图（使用重命名后的数据集）--------------------------------
p <- ggplot(df_melt, aes(x = Gene, y = Expression, fill = group)) +
  
  # 1. 箱线图美化
  geom_boxplot(
    outlier.size = 1.5,
    outlier.shape = 21,
    outlier.fill = "white",
    outlier.color = "black",
    outlier.alpha = 0.7,
    width = 0.7,
    alpha = 0.8,
    lwd = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  
  # 2. 分面显示（使用重命名后的Cohort）
  facet_wrap(
    ~ Cohort_renamed,  # 使用重命名后的变量
    nrow = 1, 
    scales = "free_y"
  ) +
  
  # 3. 显著性标记（使用重命名后的统计结果）
  geom_text(
    data = stats_df %>% filter(significance != "ns"),
    aes(x = Gene, y = y_max, label = significance),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold",
    vjust = -0.5
  ) +
  
  # 4. 标签和标题
  labs(
    x = "Gene Symbol",
    y = "Normalized Expression Level",
    fill = "Clinical Group",
    title = "Differential Gene Expression Across Cohorts",
    subtitle = "Wilcoxon rank-sum test: *p<0.05, **p<0.01, ***p<0.001"
  ) +
  
  # 5. 颜色和填充
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  
  # 6. 主题美化 - 无网格线版本（修正margin使用方式）
  theme_classic(base_size = 13) +
  theme(
    # 主标题
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5# 只指定下边距
    ),
    
    # 副标题
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5,
      color = "gray40"  # 只指定下边距
    ),
    
    # 坐标轴
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "italic",
      size = 12,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 11, 
      color = "black"
    ),
    axis.title.x = element_text(
      size = 14,
      face = "bold"
      # 移除了错误的margin参数
    ),
    axis.title.y = element_text(
      size = 14,
      face = "bold"
      # 移除了错误的margin参数
    ),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    
    # 分面标签
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    strip.background = element_rect(
      fill = "gray90",
      color = NA,
      linewidth = 0
    ),
    
    # 图例
    legend.position = "top",
    legend.title = element_text(
      size = 12,
      face = "bold"
      # 移除了错误的margin参数
    ),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "cm"),
    legend.box.spacing = unit(0.2, "cm"),
    
    # 完全去掉所有网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # 背景
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    # 分面间距
    panel.spacing = unit(1, "lines"),
    

  )

## 8) 保存数据 ----------------------------------------
# 保存包含原始和重命名信息的数据
df_melt_output <- df_melt %>%
  mutate(
    Cohort_original = Cohort,
    Cohort_display = Cohort_renamed
  ) %>%
  select(sample, Cohort_original, Cohort_display, group, Gene, Expression)

write.table(df_melt_output, file = "expr_4gene_with_renamed_cohorts.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")

# 保存统计结果
stats_output <- stats_df %>%
  select(Cohort_renamed, Gene, p_value, p_label, significance)

write.table(stats_output, file = "statistics_4gene_renamed.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")

## 9) 打印统计摘要 ----------------------------------------
cat("\n统计检验结果摘要：\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
for (cohort_name in cohort_order) {
  if (cohort_name %in% stats_df$Cohort_renamed) {
    cat(sprintf("\n%s:\n", cohort_name))
    cohort_stats <- stats_df %>% filter(Cohort_renamed == cohort_name)
    for (i in 1:nrow(cohort_stats)) {
      cat(sprintf("  %s: %s (p=%s)\n", 
                  cohort_stats$Gene[i],
                  cohort_stats$significance[i],
                  ifelse(cohort_stats$p_value[i] < 0.001, 
                         "<0.001", 
                         sprintf("%.4f", cohort_stats$p_value[i]))))
    }
  }
}

## 10) 保存图形 ----------------------------------------
n_cohorts <- length(unique(df_melt$Cohort_renamed))
plot_width <- max(8, n_cohorts * 3.5)
plot_height <- 6

pdf(file.path(fig.path, "Boxplot_4genes_renamed_cohorts.pdf"),
    width = plot_width, 
    height = plot_height,
    useDingbats = FALSE)
print(p)
dev.off()

# 保存PNG
png(file.path(fig.path, "Boxplot_4genes_renamed_cohorts.png"),
    width = plot_width * 300,
    height = plot_height * 300,
    res = 300)
print(p)
dev.off()

cat("\n✅ 图形已保存！\n")
cat("📊 PDF文件: Boxplot_4genes_renamed_cohorts.pdf\n")
cat("📊 PNG文件: Boxplot_4genes_renamed_cohorts.png\n")
cat("📈 统计结果: statistics_4gene_renamed.txt\n")
cat("📋 数据: expr_4gene_with_renamed_cohorts.txt\n")


##############################################################
cat("箱线图（含显著性星号）已输出：",
    file.path(fig.path, "Boxplot_5genes_by_cohort_case_control_with_signif.pdf"), "\n")

save.image(file="ML_all20251203.RData")

cat("\n" + paste(rep("=", 60), collapse = "") + "\n")
cat("每个基因在每个数据集中的统计检验结果（Wilcoxon秩和检验）:\n")
cat(paste(rep("=", 60), collapse = "") + "\n\n")

# 计算详细的统计信息
detailed_stats <- df_melt %>%
  group_by(Cohort, Gene) %>%
  summarise(
    # 样本量
    n_control = sum(group == "Control"),
    n_case = sum(group == "Case"),
    # 均值 ± 标准差
    mean_control = mean(Expression[group == "Control"]),
    sd_control = sd(Expression[group == "Control"]),
    mean_case = mean(Expression[group == "Case"]),
    sd_case = sd(Expression[group == "Case"]),
    # 中位数 ± 四分位距
    median_control = median(Expression[group == "Control"]),
    iqr_control = IQR(Expression[group == "Control"]),
    median_case = median(Expression[group == "Case"]),
    iqr_case = IQR(Expression[group == "Case"]),
    # 检验统计量
    p_value = tryCatch(
      wilcox.test(Expression ~ group)$p.value,
      error = function(e) NA
    ),
    test_statistic = tryCatch(
      wilcox.test(Expression ~ group)$statistic,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    # 计算log2 fold change (Case/Control)
    log2FC = log2((mean_case + 1e-10) / (mean_control + 1e-10)),
    # 显著性标注
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value >= 0.05 ~ "ns",
      TRUE ~ "NA"
    ),
    # p值格式化
    p_value_formatted = ifelse(
      p_value < 0.001, "<0.001",
      sprintf("%.4f", p_value)
    )
  ) %>%
  arrange(Cohort, Gene)

# 打印详细结果
for (cohort in unique(detailed_stats$Cohort)) {
  cat(sprintf("\n数据集: %s\n", cohort))
  cat(paste(rep("-", 60), collapse = "") + "\n")
  
  cohort_data <- detailed_stats %>% filter(Cohort == cohort)
  
  for (i in 1:nrow(cohort_data)) {
    gene <- cohort_data$Gene[i]
    cat(sprintf("基因: %s\n", gene))
    cat(sprintf("  样本量: Control=%d, Case=%d\n", 
                cohort_data$n_control[i], cohort_data$n_case[i]))
    cat(sprintf("  均值: Control=%.3f±%.3f, Case=%.3f±%.3f\n",
                cohort_data$mean_control[i], cohort_data$sd_control[i],
                cohort_data$mean_case[i], cohort_data$sd_case[i]))
    cat(sprintf("  中位数: Control=%.3f(IQR=%.3f), Case=%.3f(IQR=%.3f)\n",
                cohort_data$median_control[i], cohort_data$iqr_control[i],
                cohort_data$median_case[i], cohort_data$iqr_case[i]))
    cat(sprintf("  Log2FC(Case/Control): %.3f\n", cohort_data$log2FC[i]))
    cat(sprintf("  Wilcoxon W统计量: %.1f\n", cohort_data$test_statistic[i]))
    cat(sprintf("  P值: %s %s\n\n", 
                cohort_data$p_value_formatted[i], cohort_data$significance[i]))
  }
}


## ========= 1. 定义 4 个核心基因 =========
genes_core <- c("C4B", "AQP1", "HBA2", "ACSL6")  # 按实际基因改

# 检查它们是否在训练集中存在
genes_core <- genes_core[genes_core %in% colnames(Train_set)]
if (length(genes_core) < 4) {
  warning("4 个核心基因有不在 Train_set 中的，请检查基因名。当前有效基因：",
          paste(genes_core, collapse = ", "))
}

## ========= 2. 生成组合（单基因/双基因/三基因/四基因） =========
combo_list <- list()

# 单基因
for (i in seq_along(genes_core)) {
  g <- genes_core[i]
  combo_list[[g]] <- g
}

# 两两组合
cmb2 <- combn(genes_core, 2, simplify = FALSE)
for (v in cmb2) {
  nm <- paste(v, collapse = "+")
  combo_list[[nm]] <- v
}

# 三基因组合（按你说的：a+b+c 和 b+c+d）
# 如果你想要所有 C(4,3) 组合，可以直接用 combn(genes_core, 3, ...)
if (length(genes_core) == 4) {
  # a+b+c
  combo_list[[paste(genes_core[1:3], collapse = "+")]] <- genes_core[1:3]
  # b+c+d
  combo_list[[paste(genes_core[2:4], collapse = "+")]] <- genes_core[2:4]
}

# 四基因组合
combo_list[[paste(genes_core, collapse = "+")]] <- genes_core

# 去重复（防止三基因那里与前面重名）
combo_list <- combo_list[!duplicated(names(combo_list))]

library(pROC)

## ========= 3. 构建“队列列表” =========

# 训练集
cohort_list <- list()
cohort_list[[Train_name]] <- list(
  X = Train_set,
  y = as.numeric(Train_class[[classVar]])
)

# 验证集中各个 Cohort
val_cohorts <- unique(Test_class$Cohort)

for (co in val_cohorts) {
  idx <- which(Test_class$Cohort == co)
  if (length(idx) < 3) next  # 太少了就跳过，防止 ROC 不稳定
  
  Xc <- Test_set[idx, , drop = FALSE]
  yc <- as.numeric(Test_class[idx, classVar])
  
  # 至少要有 0 和 1 两类，否则 ROC 无法计算
  if (length(unique(yc)) < 2) next
  
  cohort_list[[co]] <- list(X = Xc, y = yc)
}

cohort_names <- names(cohort_list)

## ========= 4. 计算各组合在各队列的 AUC =========

combo_names <- names(combo_list)

AUC_combo <- matrix(
  NA,
  nrow = length(combo_names),
  ncol = length(cohort_names),
  dimnames = list(combo_names, cohort_names)
)

for (cmb_name in combo_names) {
  genes <- combo_list[[cmb_name]]
  
  for (co in cohort_names) {
    X <- cohort_list[[co]]$X
    y <- cohort_list[[co]]$y
    
    # 该队列里没有所有基因 → 记 NA
    if (!all(genes %in% colnames(X))) {
      AUC_combo[cmb_name, co] <- NA
      next
    }
    
    df <- data.frame(
      outcome = y,
      X[, genes, drop = FALSE]
    )
    
    # 拟合 logistic 回归
    fit <- try(glm(outcome ~ ., data = df, family = binomial()), silent = TRUE)
    if (inherits(fit, "try-error")) {
      AUC_combo[cmb_name, co] <- NA
      next
    }
    
    # 预测风险得分
    rs <- try(predict(fit, type = "link"), silent = TRUE)
    if (inherits(rs, "try-error")) {
      AUC_combo[cmb_name, co] <- NA
      next
    }
    
    # 计算 AUC
    roc_obj <- try(roc(response = y, predictor = as.numeric(rs), quiet = TRUE), silent = TRUE)
    if (inherits(roc_obj, "try-error")) {
      AUC_combo[cmb_name, co] <- NA
      next
    }
    
    AUC_combo[cmb_name, co] <- as.numeric(auc(roc_obj))
  }
}

## ========= 5. 计算每个组合的平均 AUC =========

# 你如果只想在“血液队列”上求平均，可以自己挑列名，例如：
# blood_cols <- c(Train_name, "GSE32127")  # 举例
# Mean_AUC   <- apply(AUC_combo[, blood_cols, drop = FALSE], 1, mean, na.rm = TRUE)

# 如果先不区分，先对所有队列（非 NA）求平均 AUC：
Mean_AUC <- apply(AUC_combo, 1, function(x) mean(as.numeric(x), na.rm = TRUE))

## ========= 6. 组装结果表并导出 =========

AUC_combo_df <- data.frame(
  Combo   = rownames(AUC_combo),
  AUC_combo,
  Mean_AUC = Mean_AUC,
  row.names = NULL,
  check.names = FALSE
)

# 按平均 AUC 从高到低排序
AUC_combo_df <- AUC_combo_df[order(AUC_combo_df$Mean_AUC, decreasing = TRUE), ]

# 写表
write.table(
  AUC_combo_df,
  file = file.path(res.path, "AUC_gene_combinations.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)








