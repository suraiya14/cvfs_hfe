library(doParallel)
library(e1071)
library(caret)
library(ROCR)
#library(pROC)
#library(klaR)
library(mltools)
library(pROC)
library(parallel)
library(ggplot2)
df<-read.csv("D:\\Research_Work\\Disertation Project 3\\RawData\\FeatureExtraction\\hypergraph\\bin10\\data30/selected_training_merged_file.csv", header = TRUE)
#df<-read.csv("D:\\Research_Work\\Disertation Project 3\\RawData\\FeatureExtraction\\CVFS\\C3E5P0.2/selected_training_merged_file.csv", header = TRUE)
df$Output[df$Output == -1] <- 0

df_test<-read.csv("D:\\Research_Work\\Disertation Project 3\\RawData\\FeatureExtraction\\hypergraph\\bin10\\data30/selected_validation_merged_file.csv", header = TRUE)
#df_test<-read.csv("D:\\Research_Work\\Disertation Project 3\\RawData\\FeatureExtraction\\CVFS\\C3E5P0.2/selected_validation_merged_file.csv", header = TRUE)
df_test$Output[df_test$Output == -1] <- 0
#XGB

set.seed(123)

library(xgboost)

#glimpse(df)

X_train = data.matrix(df[,!colnames(df)%in%c("Output")]) # independent variables for train
X_train
y_train = as.numeric(as.character(df$Output)) #
y_train
dtrain = xgb.DMatrix(data=X_train, label=y_train)
dtrain

X_test = data.matrix(df_test[,!colnames(df_test)%in%c("Output")]) # independent variables for train
X_test
y_test = df_test$Output #

dvalid = xgb.DMatrix(data=X_test, label=y_test)
dvalid

#used XGB
cores <- detectCores() # Detect the number of available cores
cores
#cl <- makeCluster(cores - 1) # Create a cluster with all but one core
cl <- makePSOCKcluster(cores - 2)
cl
registerDoParallel(cl) # Register the cluster for parallel processing

params <- list(
  objective = "binary:logistic",
  learning_rate = 0.05,
  subsample = 0.9,
  colsample_bynode = 1,
  reg_lambda = 2,
  max_depth = 5
)


fit_xgb <- xgb.train(
  params,
  data = dtrain,
  watchlist = list(valid = dvalid),
  #early_stopping_rounds = 20,
  print_every_n = 100,
  nrounds = 10000 # early stopping
)

stopCluster(cl)

library(tidyr)
library(SHAPforxgboost)

x<-colnames(df)
x<-x[x != "Output"]

# Step 1: Select some observations
#X <- data.matrix(df[sample(nrow(df), 5000), x])
X <- data.matrix(df[sample(nrow(df), nrow(df)), x])



# Step 2: Crunch SHAP values

cores <- detectCores() # Detect the number of available cores
cores
#cl <- makeCluster(cores - 1) # Create a cluster with all but one core
cl <- makePSOCKcluster(cores - 2)
registerDoParallel(cl) # 
shap <- shap.prep(fit_xgb, X_train = X)
shap
stopCluster(cl)
shap_back<-shap

shap<-shap_back


shap_plot<-shap.plot.summary(shap)
#top first 10 data displaying
#shap.plot.summary.wrap1(fit_xgb, X, top_n = 10)
shap_plot <- shap_plot + theme(axis.title.y = element_text(face = "bold"))

# Show the modified plot
print(shap_plot)


library(caret)
#prediction
df_test$Output<-as.factor(df_test$Output)
df_test$Output
dtest <- xgb.DMatrix(data = as.matrix( X_test))

xgbPredictions <-predict(fit_xgb, dtest, type = "prob")

print(xgbPredictions)

#write.csv(xgbPredictions, file = "D:\\Research_Work\\Disertation Project 3\\RawData\\FeatureExtraction/CVFS\\probability\\probability_C2E10P0.4.csv", row.names = FALSE)

# 'results' now contains the predicted values


roc_curve <- roc(df_test$Output, xgbPredictions)
auc_value <- auc(roc_curve)
print(paste("SVM AUC:", auc_value))
# Convert predicted values to class labels
xgbpredicted_labels <- ifelse(xgbPredictions > 0.5, 1, 0)
xgbpredicted_labels
df_test$Output<-as.factor(df_test$Output)

summary(df_test$Output)



xgbpredicted_labels <- factor(xgbpredicted_labels, levels = levels(df_test$Output))

# Create the confusion matrix
xgbconfusion <- confusionMatrix(xgbpredicted_labels, df_test$Output, positive ='1', mode = "everything")

# Display the confusion matrix
print(xgbconfusion)
xgbconfusion$byClass

# Assuming 'results' contains your predicted values and 'test' is your test dataset
preds <- xgbpredicted_labels
preds
actuals <- df_test$Output
actuals

# Calculate MCC
confusion_matrix <- table(Actual = actuals, Predicted = preds)
confusion_matrix
tp <- confusion_matrix["1","1"]  # True Positives
tp
tn <- confusion_matrix["0","0"]  # True Negatives
fp <- confusion_matrix["0","1"]  # False Positives
fn <- confusion_matrix["1","0"]  # False Negatives

mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
mcc
# Print MCC

# to rename variables in the shap plot
#y_axis_labels <- attr(shap$variable, "labels")

y_axis_labels <- unique(shap$variable)
y_axis_labels


#shap.plot.summary(shap_updated)
draw_confusion_matrix <- function(cm) {
  layout(matrix(1))
  par(mar=c(2,2,2,2))
  
  # Create the matrix
  rect(200, 350, 300, 400, col='#AF97D0')
  rect(300, 350, 400, 400, col='#A7AD50')
  text(250, 405, "Reference", cex=1.3, font=2)
  text(175, 375, "Prediction", cex=1.3, srt=90, font=2)
  
  # Add in the cm results
  res <- as.numeric(cm$table)
  text(250, 375, "-1", cex=1.3, font=2, col='black')
  text(250, 350, "1", cex=1.3, font=2, col='black')
  text(300, 400, "-1", cex=1.3, font=2, col='black')
  text(300, 375, res[1], cex=1.6, font=2, col='black')
  text(300, 350, res[3], cex=1.6, font=2, col='black')
  text(400, 400, "1", cex=1.3, font=2, col='black')
  text(400, 375, res[2], cex=1.6, font=2, col='black')
  text(400, 350, res[4], cex=1.6, font=2, col='black')
}


draw_confusion_matrix <- function(cm) {
  
  layout(matrix(1))
  par(mar=c(2,2,2,2))
  #plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  plot(c(123, 345), c(300, 452), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  #title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#AF97D0')
  text(195, 435, -1, cex=1.2)
  rect(250, 430, 340, 370, col='#A7AD50')
  text(295, 435, 1, cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#A7AD50')
  rect(250, 305, 340, 365, col='#AF97D0')
  text(140, 400, -1, cex=1.2, srt=90)
  text(140, 335, 1, cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='black')
  text(195, 335, res[2], cex=1.6, font=2, col='black')
  text(295, 400, res[3], cex=1.6, font=2, col='black')
  text(295, 335, res[4], cex=1.6, font=2, col='black')
  
  
  
}  
draw_confusion_matrix(xgbconfusion)
