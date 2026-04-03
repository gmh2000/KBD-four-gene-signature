# all ML algorithms
RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){
  # for example: Enet [alpha=0.4]
  method = gsub(" ", "", method) # 
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)  # get name of ML algorithm, e.g., Enet
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method) # get parameter of ML algorithm, e.g., alpha=0.4
  
  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "Stepglm" = list("direction" = method_param),
    NULL
  )
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, ";",
          " using ", ncol(Train_set), " Variables")
  
  args = list("Train_set" = Train_set,
              "Train_label" = Train_label,
              "mode" = mode,
              "classVar" = classVar)
  args = c(args, method_param)
  
  obj <- do.call(what = paste0("Run", method_name),
                 args = args) 
  
  if(mode == "Variable"){
    message(length(obj), " Variables retained;\n")
  }else{message("\n")}
  return(obj)
}

RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
  cv.fit = cv.glmnet(x = Train_set,
                     y = Train_label[[classVar]], 
                     family = "binomial", 
                     type.measure ="auc",
                     alpha=alpha,
                     standardize = FALSE,
                     nfolds = 10)
  fit = glmnet(x = Train_set,
               y = Train_label[[classVar]], 
               family = "binomial", 
               alpha=alpha,
               lambda=cv.fit$lambda.min,
               standardize = FALSE) 
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLasso <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
}

RunRidge <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
}

RunStepglm <- function(Train_set, Train_label, mode, classVar, direction){
  fit <- step(glm(formula = Train_label[[classVar]] ~ .,
                  family  = "binomial",
                  data = as.data.frame(Train_set)), 
              k = 2, direction = direction)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunSVM <- function(Train_set, Train_label, mode, classVar){
  fit <- svm(x = Train_set, y = as.factor(Train_label[[classVar]]),
             kernel = "linear", scale=FALSE, 
             cross = 10, fitted = FALSE)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunglmBoost <- function(Train_set, Train_label, mode, classVar){
  Train_label[[classVar]] <- as.numeric(as.character(Train_label[[classVar]]))
  formula <- as.formula(paste(classVar, "~", ".", sep = ""))
  fit = glmboost(formula,
                 data = cbind(as.data.frame(Train_set),
                              Train_label), center = FALSE)
  cp = cvrisk(fit, folds = cv(model.weights(fit), 
                              type = "kfold", B = 5)) # 
  fit = fit[mstop(cp)] 
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunNB <- function(Train_set, Train_label, mode, classVar){
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  fit <- naiveBayes(x = Train_set, y = Train_label[[classVar]])
  if (mode == "Model") return(fit)
  if (mode == "Variable") {
    # NB doesn't do feature selection; return all variables
    return(colnames(Train_set))
  }
}

RunGBM <- function(Train_set, Train_label, mode, classVar){
  Train_label[[classVar]] <- as.numeric(as.character(Train_label[[classVar]]))
  d <- cbind(Train_label, Train_set)
  fit_gbm <- gbm(
    formula = as.formula(paste(classVar, "~ .")),
    data = d,
    distribution = "bernoulli",
    n.trees = 10000,
    interaction.depth = 3,
    shrinkage = 0.001,
    cv.folds = 5,
    n.cores = 1,
    verbose = FALSE
  )
  best_iter <- gbm.perf(fit_gbm, method = "cv", plot.it = FALSE)
  fit_gbm$best_iter <- best_iter
  
  if (mode == "Model") return(fit_gbm)
  if (mode == "Variable") {
    var_imp <- summary(fit_gbm, plotit = FALSE, n.trees = best_iter)
    return(var_imp$var[var_imp$rel.inf > 0])
  }
}

RunRF <- function(Train_set, Train_label, mode, classVar){
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  fit <- randomForest(x = Train_set, y = Train_label[[classVar]],
                      ntree = 500, importance = TRUE)
  if (mode == "Model") return(fit)
  if (mode == "Variable") {
    var_imp <- importance(fit, type = 2)
    return(rownames(var_imp)[var_imp[,1] > 0])
  }
}

RunXGBoost <- function(Train_set, Train_label, mode, classVar){
  label <- as.numeric(as.character(Train_label[[classVar]]))
  dtrain <- xgb.DMatrix(data = as.matrix(Train_set), label = label)
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 2,
    eta = 1
  )
  nround <- 10
  
  fit <- xgb.train(params = params, data = dtrain,
                   nrounds = nround, verbose = 0)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") {
    imp <- xgb.importance(model = fit)
    return(imp$Feature)
  }
}

RunLDA <- function(Train_set, Train_label, mode, classVar){
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  fit <- lda(x = Train_set, grouping = Train_label[[classVar]])
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(colnames(Train_set))
}

RunplsRglm <- function(Train_set, Train_label, mode, classVar){
  Train_label[[classVar]] <- as.numeric(as.character(Train_label[[classVar]]))
  fit <- plsRglm(
    Y = Train_label[[classVar]],
    X = Train_set,
    nt = 10,
    modele = "pls-glm-family",
    family = binomial()
  )
  if (mode == "Model") return(fit)
  if (mode == "Variable") {
    coefs <- coef(fit)
    vars <- names(coefs)[coefs != 0]
    vars <- vars[vars != "(Intercept)"]
    return(vars)
  }
}

# -------------------------------------------------------------
# 工具函数：提取变量，计算预测值，评估 AUC，绘制热图
# -------------------------------------------------------------
ExtractVar <- function(fit){
  if (inherits(fit, "glmnet")) {
    coefs <- coef(fit)
    vars <- rownames(coefs)[which(coefs != 0)]
    vars <- vars[vars != "(Intercept)"]
  } else if (inherits(fit, "glm")) {
    coefs <- coef(fit)
    vars <- names(coefs)[coefs != 0]
    vars <- vars[vars != "(Intercept)"]
  } else if (inherits(fit, "glmboost")) {
    coefs <- coef(fit, off2int = TRUE)
    vars <- names(coefs)[coefs != 0]
    vars <- vars[vars != "(Intercept)"]
  } else if (inherits(fit, "randomForest")) {
    # 这里前面已经改过：importance 的“行名”才是变量名
    vars <- rownames(fit$importance)
  } else if (inherits(fit, "gbm")) {
    var_imp <- summary(fit, plotit = FALSE, n.trees = fit$best_iter)
    vars <- var_imp$var[var_imp$rel.inf > 0]
  } else if (inherits(fit, "xgb.Booster")) {
    imp <- xgb.importance(model = fit)
    vars <- imp$Feature
  } else if (inherits(fit, "svm")) {
    # ★ 关键修改：用支持向量矩阵的列名作为特征名
    vars <- colnames(fit$SV)
  } else if (inherits(fit, "lda")) {
    vars <- rownames(fit$scaling)
  } else if (inherits(fit, "plsRglm")) {
    coefs <- coef(fit)
    vars <- names(coefs)[coefs != 0]
    vars <- vars[vars != "(Intercept)"]
  } else {
    vars <- NULL
  }
  return(vars)
}


CalPredictScore <- function(fit, new_data){
  # 统一计算风险评分
  if (inherits(fit, "glmnet")) {
    RS <- predict(fit, newx = as.matrix(new_data), type = "response")[,1]
  } else if (inherits(fit, "glm")) {
    RS <- predict(fit, newdata = as.data.frame(new_data), type = "response")
  } else if (inherits(fit, "glmboost")) {
    RS <- predict(fit, newdata = as.data.frame(new_data), type = "response")
  } else if (inherits(fit, "randomForest")) {
    RS <- predict(fit, newdata = as.data.frame(new_data), type = "prob")[,2]
  } else if (inherits(fit, "gbm")) {
    RS <- predict(fit, newdata = as.data.frame(new_data),
                  n.trees = fit$best_iter, type = "response")
  } else if (inherits(fit, "xgb.Booster")) {
    dtest <- xgb.DMatrix(data = as.matrix(new_data))
    RS <- predict(fit, newdata = dtest)
  } else if (inherits(fit, "svm")) {
    RS <- attr(predict(fit, newdata = new_data, decision.values = TRUE), 
               "decision.values")
  } else if (inherits(fit, "lda")) {
    RS <- predict(fit, newdata = new_data)$posterior[,2]
  } else if (inherits(fit, "plsRglm")) {
    RS <- predict(fit, newdata = as.data.frame(new_data), type = "response")
  } else if (inherits(fit, "naiveBayes")) {
    RS <- predict(fit, newdata = new_data, type = "raw")[,2]
  } else {
    stop("Unsupported model type in CalPredictScore.")
  }
  return(RS)
}

PredictClass <- function(fit, new_data, cutoff = 0.5){
  RS <- CalPredictScore(fit, new_data)
  pred <- ifelse(RS >= cutoff, 1, 0)
  return(pred)
}

RunEval <- function(fit, Test_set, Test_label,
                    Train_set = NULL, Train_label = NULL,
                    Train_name = "Training",
                    cohortVar = "Cohort", classVar = "outcome") {
  # 1) 针对 glmnet：使用训练时的“全部变量维度”
  if (inherits(fit, "glmnet")) {
    feats <- rownames(fit$beta)
  } else {
    feats <- ExtractVar(fit)
  }
  
  # 2) 如果模型没有返回特征，直接返回 NA
  if (is.null(feats) || length(feats) == 0) {
    warning("RunEval: no features extracted from model; return NA AUCs.")
    if (!is.null(Train_set) && !is.null(Train_label)) {
      cohorts <- c(Train_name, as.character(unique(Test_label[[cohortVar]])))
    } else {
      cohorts <- as.character(unique(Test_label[[cohortVar]]))
    }
    out <- rep(NA_real_, length(cohorts))
    names(out) <- cohorts
    return(out)
  }
  
  # 3) 确保特征名存在于数据矩阵中
  if (!is.null(Train_set)) {
    feats <- intersect(feats, colnames(Train_set))
  } else {
    feats <- intersect(feats, colnames(Test_set))
  }
  
  if (length(feats) == 0) {
    warning("RunEval: extracted features not found in data; return NA AUCs.")
    if (!is.null(Train_set) && !is.null(Train_label)) {
      cohorts <- c(Train_name, as.character(unique(Test_label[[cohortVar]])))
    } else {
      cohorts <- as.character(unique(Test_label[[cohortVar]]))
    }
    out <- rep(NA_real_, length(cohorts))
    names(out) <- cohorts
    return(out)
  }
  
  # 4) 构建 new_data 和标签矩阵
  if (!is.null(Train_set) && !is.null(Train_label)) {
    Test_label$Cohort  <- Test_label[[cohortVar]]
    Train_label$Cohort <- Train_name
    
    new_data <- rbind.data.frame(
      Train_set[, feats, drop = FALSE],
      Test_set[, feats, drop = FALSE]
    )
    label_all <- rbind(
      Train_label[, c("Cohort", classVar), drop = FALSE],
      Test_label[, c("Cohort", classVar), drop = FALSE]
    )
  } else {
    new_data <- Test_set[, feats, drop = FALSE]
    label_all <- Test_label[, c(cohortVar, classVar), drop = FALSE]
    colnames(label_all)[1] <- "Cohort"
  }
  
  # 5) 统一计算风险分数 —— 这里加 tryCatch，任何预测报错全记 NA
  RS_all <- tryCatch(
    CalPredictScore(fit, new_data),
    error = function(e) {
      warning("RunEval: prediction failed for this model: ",
              conditionMessage(e))
      return(rep(NA_real_, nrow(new_data)))
    }
  )
  
  # 6) 分 cohort 计算各自 AUC
  AUC_vec <- c()
  for (co in unique(label_all$Cohort)) {
    idx <- which(label_all$Cohort == co)
    # 如果这一队列预测值全 NA，就直接 NA
    if (all(is.na(RS_all[idx]))) {
      AUC_vec[co] <- NA_real_
      next
    }
    if (length(unique(label_all[idx, classVar])) == 2) {
      roc_obj <- roc(
        response  = label_all[idx, classVar],
        predictor = RS_all[idx],
        quiet     = TRUE
      )
      AUC_vec[co] <- as.numeric(auc(roc_obj))
    } else {
      AUC_vec[co] <- NA_real_
    }
  }
  return(AUC_vec)
}

SimpleHeatmap <- function(Cindex_mat, avg_Cindex, CohortCol,
                          barCol = "steelblue",
                          cellwidth = 1, cellheight = 0.5,
                          cluster_columns = FALSE, cluster_rows = FALSE){
  # 绘制算法 × 队列的 AUC 热图，并在右侧添加平均 AUC 条形图 + 数字标签
  
  # 平均 AUC 的数字标签
  avg_label <- sprintf("%.3f", avg_Cindex)
  
  row_ha <- rowAnnotation(
    "Average AUC" = anno_barplot(
      avg_Cindex,
      gp = gpar(fill = barCol, col = NA),
      axis_param = list(side = "bottom"),
      width = unit(2, "cm")
    ),
    " " = anno_text(                     # 在条形图右侧加数字
      avg_label,
      just     = "left",
      location = 0.2,
      gp       = gpar(fontsize = 8),
      width    = unit(1, "cm")
    )
  )
  
  col_ha <- HeatmapAnnotation(
    "Cohort" = colnames(Cindex_mat),
    col = list("Cohort" = CohortCol),
    annotation_legend_param = list(
      "Cohort" = list(title = "Cohort")
    ),
    show_annotation_name = FALSE
  )
  
  Heatmap(as.matrix(Cindex_mat), name = "AUC",
          right_annotation = row_ha, 
          top_annotation   = col_ha,
          col = c("#4195C1", "#FFFFFF", "#CB5746"),
          rect_gp = gpar(col = "black", lwd = 1),
          cluster_columns = cluster_columns, cluster_rows = cluster_rows,
          show_column_names = FALSE, 
          show_row_names    = TRUE,
          row_names_side    = "left",
          width  = unit(cellwidth  * ncol(Cindex_mat) + 2, "cm"),
          height = unit(cellheight * nrow(Cindex_mat), "cm"),
          column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
          column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, col) {
            grid.text(
              label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
              x, y, gp = gpar(fontsize = 10)
            )
          }
  )
}
