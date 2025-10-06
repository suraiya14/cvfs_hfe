# ============================================================
# XGBoost + SHAP Workflow (Train/Test, Plots, Grouped SHAP, Metrics)
# ------------------------------------------------------------
# - Loads training/validation CSVs
# - Trains XGBoost (binary:logistic) with watchlist
# - Computes SHAP values and summary plots
# - Aggregates mean |SHAP| per feature family (dist/comp/tran -> CTD)
# - Evaluates AUC, confusion matrix, MCC
# - Draws a confusion matrix figure
#
# NOTE: Functionality, paths, and hyperparameters are unchanged.
#       Parallelization uses (cores - 2) as in the original.
# ============================================================

# ---------------------------
# Libraries
# ---------------------------
library(doParallel)
library(e1071)
library(caret)
library(ROCR)
library(mltools)
library(pROC)
library(parallel)
library(ggplot2)

# ---------------------------
# Load data
# ---------------------------
# Training data; relabel -1 -> 0 (kept)
df<-read.csv("D:\\Research_Work\\Disertation Project 3\\RawData\\FeatureExtraction\\hypergraph\\bin10\\data30/selected_training_merged_file.csv", header = TRUE)
df$Output[df$Output == -1] <- 0

# Validation data; relabel -1 -> 0 (kept)
df_test<-read.csv("D:\\Research_Work\\Disertation Project 3\\RawData\\FeatureExtraction\\hypergraph\\bin10\\data30/selected_validation_merged_file.csv", header = TRUE)
df_test$Output[df_test$Output == -1] <- 0

# ---------------------------
# XGBoost setup
# ---------------------------
set.seed(123)
library(xgboost)

# Matrix construction (exclude label column "Output")
X_train <- data.matrix(df[, !colnames(df) %in% c("Output")])
y_train <- as.numeric(as.character(df$Output))
dtrain  <- xgb.DMatrix(data = X_train, label = y_train)

X_test <- data.matrix(df_test[, !colnames(df_test) %in% c("Output")])
y_test <- df_test$Output
dvalid <- xgb.DMatrix(data = X_test, label = y_test)

# Parallel backend (kept exactly: PSOCK, cores - 2)
cores <- detectCores()
cl <- makePSOCKcluster(cores - 2)
registerDoParallel(cl)

# XGBoost params (unchanged)
params <- list(
  objective = "binary:logistic",
  learning_rate = 0.05,
  subsample = 0.9,
  colsample_bynode = 1,
  reg_lambda = 2,
  max_depth = 5
)

# Training (unchanged watchlist, nrounds, print interval)
fit_xgb <- xgb.train(
  params,
  data = dtrain,
  watchlist = list(valid = dvalid),
  print_every_n = 100,
  nrounds = 10000
)

# Stop cluster (first phase)
stopCluster(cl)

# ---------------------------
# SHAP values
# ---------------------------
library(tidyr)
library(SHAPforxgboost)

# Feature names (exclude label)
x <- colnames(df)
x <- x[x != "Output"]

# Randomized row order for SHAP matrix (kept)
X <- data.matrix(df[sample(nrow(df), nrow(df)), x])

# Parallel backend again for SHAP (unchanged)
cores <- detectCores()
cl <- makePSOCKcluster(cores - 2)
registerDoParallel(cl)

# Compute SHAP (unchanged)
shap <- shap.prep(fit_xgb, X_train = X)

# Stop cluster (second phase)
stopCluster(cl)

# ---------------------------
# SHAP summary plots
# ---------------------------
shap_plot <- shap.plot.summary(shap)
shap_plot <- shap_plot + theme(axis.title.y = element_text(face = "bold"))
print(shap_plot)

# Top-10 features plot (kept)
shap.plot.summary.wrap1(fit_xgb, X, top_n = 10)

# ---------------------------
# NEW: Per-family summed mean |SHAP|
#  - Combine dist/comp/tran into CTD (kept)
# ---------------------------
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

shap_dt <- as.data.table(shap)

# Determine SHAP value column name defensively (unchanged logic)
shap_value_col <- if ("contribution" %in% names(shap_dt)) {
  "contribution"
} else if ("value" %in% names(shap_dt)) {
  "value"
} else if ("phi" %in% names(shap_dt)) {
  "phi"
} else {
  stop("Could not find SHAP value column (expected one of: contribution, value, phi)")
}

# Mean |SHAP| per feature, then map to family and collapse
per_feature <- shap_dt %>%
  as_tibble() %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = mean(abs(.data[[shap_value_col]]), na.rm = TRUE), .groups = "drop") %>%
  mutate(family = str_replace(variable, "_.*$", "")) %>%
  mutate(family = case_when(
    family %in% c("dist", "comp", "tran") ~ "CTD",
    TRUE ~ family
  ))

per_family <- per_feature %>%
  group_by(family) %>%
  summarise(per_family_summed_mean_abs_shap = sum(mean_abs_shap), .groups = "drop") %>%
  arrange(desc(per_family_summed_mean_abs_shap))

print(per_family)

# Bar plot (unchanged styling except the original choices)
if (nrow(per_family) > 0) {
  gg_per_family <- ggplot(per_family,
                          aes(x = reorder(family, per_family_summed_mean_abs_shap),
                              y = per_family_summed_mean_abs_shap)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(x = "Feature family",
         y = "Per-family summed mean |SHAP|") +
    # title was commented in the original; keep that behavior
    theme_minimal(base_size = 13) +
    theme(axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5))
  print(gg_per_family)
} else {
  cat("⚠️ No data in per_family — check SHAP column name (value/contribution/phi).\n")
}

# ---------------------------
# Prediction + Metrics
# ---------------------------
library(caret)

# Model predictions on validation set (unchanged)
df_test$Output <- as.factor(df_test$Output)
dtest <- xgb.DMatrix(data = as.matrix(X_test))
xgbPredictions <- predict(fit_xgb, dtest, type = "prob")
print(xgbPredictions)

# ROC/AUC (unchanged)
roc_curve <- roc(df_test$Output, xgbPredictions)
auc_value <- auc(roc_curve)
print(paste("XGB AUC:", auc_value))

# Class predictions with 0.5 threshold (kept)
xgbpredicted_labels <- ifelse(xgbPredictions > 0.5, 1, 0)
xgbpredicted_labels <- factor(xgbpredicted_labels, levels = levels(df_test$Output))

# Confusion matrix (kept: positive='1', mode='everything')
xgbconfusion <- confusionMatrix(xgbpredicted_labels, df_test$Output, positive = '1', mode = "everything")
print(xgbconfusion)
xgbconfusion$byClass

# MCC (manual calc preserved)
confusion_matrix <- table(Actual = df_test$Output, Predicted = xgbpredicted_labels)
tp <- confusion_matrix["1","1"]
tn <- confusion_matrix["0","0"]
fp <- confusion_matrix["0","1"]
fn <- confusion_matrix["1","0"]
mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
mcc

# ---------------------------
# Confusion-matrix drawing (unchanged)
# ---------------------------
draw_confusion_matrix <- function(cm) {
  layout(matrix(1))
  par(mar = c(2, 2, 2, 2))
  plot(c(123, 345), c(300, 452), type = "n", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  rect(150, 430, 240, 370, col = '#AF97D0')
  text(195, 435, -1, cex = 1.2)
  rect(250, 430, 340, 370, col = '#A7AD50')
  text(295, 435, 1, cex = 1.2)
  text(125, 370, 'Predicted', cex = 1.3, srt = 90, font = 2)
  text(245, 450, 'Actual', cex = 1.3, font = 2)
  rect(150, 305, 240, 365, col = '#A7AD50')
  rect(250, 305, 340, 365, col = '#AF97D0')
  text(140, 400, -1, cex = 1.2, srt = 90)
  text(140, 335, 1, cex = 1.2, srt = 90)
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex = 1.6, font = 2, col = 'black')
  text(195, 335, res[2], cex = 1.6, font = 2, col = 'black')
  text(295, 400, res[3], cex = 1.6, font = 2, col = 'black')
  text(295, 335, res[4], cex = 1.6, font = 2, col = 'black')
}

# Draw using the computed confusion matrix
draw_confusion_matrix(xgbconfusion)
