if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(pROC)
library(dplyr)
# 假设你已经有：
# fit（最终模型）
# Train_set / Test_set / Train_label / Test_label

Test_label$Cohort <- Test_label$Cohort
Train_label$Cohort <- "Training"

new_data <- rbind(
  Train_set[, feats],
  Test_set[, feats]
)

label_all <- rbind(
  Train_label[, c("Cohort", "outcome")],
  Test_label[, c("Cohort", "outcome")]
)

RS_all <- CalPredictScore(fit, new_data)

all_pred <- data.frame(
  cohort = label_all$Cohort,
  outcome = label_all$outcome,
  score = RS_all
)
range(scores_all)
RS_all <- CalPredictScore(fit, new_data)
range(RS_list$`Stepglm [both]+LDA`)
score <- RS_list$`Stepglm [both]+LDA`
class(score)
label_all=rbind(test_anno,train_anno)
length(score)
nrow(label_all)
test_pred <- label_all
test_pred$score <- as.numeric(score[match(test_pred$sample, names(score))])
head(test_pred)
sum(is.na(test_pred$score))
 
library(pROC)
library(dplyr)

eval_binary_model <- function(df,
                              outcome_col = "outcome",
                              score_col = "score",
                              cohort_name = "Cohort",
                              ci_method = c("delong", "bootstrap"),
                              boot_n = 2000) {
  ci_method <- match.arg(ci_method)
  
  dat <- df[, c(outcome_col, score_col)]
  colnames(dat) <- c("outcome", "score")
  dat <- dat[complete.cases(dat), ]
  
  dat$outcome <- as.numeric(as.character(dat$outcome))
  if (!all(dat$outcome %in% c(0, 1))) {
    stop("outcome must be coded as 0/1.")
  }
  
  roc_obj <- roc(
    response = dat$outcome,
    predictor = dat$score,
    quiet = TRUE,
    direction = "<"
  )
  
  auc_val <- as.numeric(auc(roc_obj))
  
  if (ci_method == "delong") {
    auc_ci <- ci.auc(roc_obj, method = "delong")
  } else {
    auc_ci <- ci.auc(roc_obj, method = "bootstrap", boot.n = boot_n)
  }
  
  # 强制转为纯numeric向量
  auc_ci <- as.numeric(auc_ci)
  
  best_coords <- coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  # coords有时返回data.frame/list，这里统一转numeric
  threshold_val   <- as.numeric(best_coords[1, "threshold"])
  sensitivity_val <- as.numeric(best_coords[1, "sensitivity"])
  specificity_val <- as.numeric(best_coords[1, "specificity"])
  
  out <- data.frame(
    Cohort = cohort_name,
    AUC = round(auc_val, 3),
    AUC_lower_95CI = round(auc_ci[1], 3),
    AUC_upper_95CI = round(auc_ci[3], 3),
    Cutoff = round(threshold_val, 4),
    Sensitivity = round(sensitivity_val, 3),
    Specificity = round(specificity_val, 3),
    stringsAsFactors = FALSE
  )
  
  return(list(
    roc = roc_obj,
    summary = out
  ))
}

eval_multiple_cohorts <- function(all_data,
                                  cohort_col = "cohort",
                                  outcome_col = "outcome",
                                  score_col = "score",
                                  ci_method = "delong",
                                  boot_n = 2000) {
  cohort_list <- split(all_data, all_data[[cohort_col]])
  
  res_list <- lapply(names(cohort_list), function(nm) {
    eval_binary_model(
      df = cohort_list[[nm]],
      outcome_col = outcome_col,
      score_col = score_col,
      cohort_name = nm,
      ci_method = ci_method,
      boot_n = boot_n
    )$summary
  })
  
  res_table <- bind_rows(res_list)
  return(res_table)
}


format_auc_ci <- function(df) {
  df %>%
    mutate(
      `AUC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", AUC, AUC_lower_95CI, AUC_upper_95CI),
      Sensitivity = sprintf("%.1f%%", Sensitivity * 100),
      Specificity = sprintf("%.1f%%", Specificity * 100)
    ) %>%
    select(Cohort, `AUC (95% CI)`, Cutoff, Sensitivity, Specificity)
}



results_table <- eval_multiple_cohorts(
  all_data = test_pred,
  cohort_col = "Cohort",
  outcome_col = "outcome",
  score_col = "score",
  ci_method = "delong"
)

final_table <- format_auc_ci(results_table)
print(final_table)
head(test_pred)



table(test_pred$outcome)
boxplot(score ~ outcome, data = test_pred)
score
cor(test_pred$outcome, test_pred$score)
table(test_pred$Cohort, test_pred$outcome)
boxplot(score ~ outcome, data = test_pred)
head(test_pred)





if (!requireNamespace("rmda", quietly = TRUE)) install.packages("rmda")
library(rmda)
train_data <- subset(test_pred, Cohort == "GSE59446")
dca_train <- decision_curve(
  outcome ~ score,
  data = train_data,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = FALSE,
  study.design = "cohort"
)
png("DCA_Training_SCI.png", width = 2200, height = 1700, res = 300)

par(
  bg = "white",
  mar = c(5, 5, 2, 2),
  bty = "l",
  las = 1,
  cex.lab = 1.3,
  cex.axis = 1.1,
  family = "sans"
)

plot_decision_curve(
  dca_train,
  curve.names = "Four-gene model",
  xlab = "Threshold probability",
  ylab = "Net benefit",
  legend.position = "topright",
  col = c("#C00000", "grey55", "black"),
  lwd = 2,
  confidence.intervals = FALSE,
  standardize = FALSE
)

dev.off()

library(rmda)

















val_data <- subset(test_pred, Cohort == "GSE32127")

dca_val <- decision_curve(
  outcome ~ score,
  data = val_data,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = FALSE,
  study.design = "cohort"
)

png("DCA_validation_SCI.png", width = 2200, height = 1700, res = 300)

par(
  bg = "white",
  mar = c(5, 5, 2, 2),
  bty = "l",
  las = 1,
  cex.lab = 1.3,
  cex.axis = 1.1,
  family = "sans"
)

plot_decision_curve(
  dca_val,
  curve.names = "Four-gene model",
  xlab = "Threshold probability",
  ylab = "Net benefit",
  legend.position = "topright",
  col = c("#C00000", "grey55", "black"),
  lwd = 2,
  confidence.intervals = FALSE,
  standardize = FALSE
)

dev.off()



