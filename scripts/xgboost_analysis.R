# Load the required packages
library(tidyverse)
library(caret)
library(randomForest)
library(xgboost)
library(pROC)

set.seed(1234)
# Load the data
data <- read.csv2("~/Documentos/Thaina/gene_expression_dct_FCGR1A_only.csv", sep = ",") # gene_expression_dct.csv
View(data)

# gene_expression_dct_noIRAK3.csv

#data$diag <- as.factor(data$diag)
data$diag <- factor(data$diag, levels = c("NO_TB", "TB"), ordered = TRUE)

# Split the data into training and testing sets
trainIndex <- createDataPartition(data$diag, p = 0.8, list = FALSE)
trainData <- data[trainIndex,]
testData <- data[-trainIndex,]

# Cross-validation ----
fitControl <- trainControl(method = "cv", number = 5, savePredictions = 'final',
                           classProbs = TRUE, summaryFunction=twoClassSummary)


# SVM Model based analysis ----
xgbGrid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(3, 6),
  eta = c(0.01, 0.1),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

set.seed(123)
xgb_model <- train(diag ~ ., data = trainData, method = "xgbTree", trControl = fitControl,
  tuneGrid = xgbGrid, metric = "ROC")

# Make predictions on test set
predictions <- predict(xgb_model, newdata = testData)
prob_predictions <- predict(xgb_model, newdata = testData, type = "prob")

# Calculate confusion matrix
confMatrix <- confusionMatrix(predictions, testData$diag)
print(confMatrix)

# Feature importance
importance <- varImp(xgb_model)
plot(importance)


phase3 <- read.csv2("~/Documentos/Thaina/gene_expression_phase3_FCGR1A_only.csv", sep = ",")

predictions <- predict(xgb_model, newdata = phase3)

phase3$predicted_class <- predictions
View(phase3)

write.csv(phase3, "FCGR1A/phase3_predictions_xgb_FCGR1A.csv")

# FCGR1A, PDCD1LG2, GBP5
data$sig <- data$FCGR1A+data$PDCD1LG2+data$GBP5
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "FCGR1A+PDCD1LG2+GBP5", print.thres = "best", col = "blue", legacy.axes = TRUE)
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

# FCGR1A, PDCD1LG2
data$sig <- data$FCGR1A+data$PDCD1LG2
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "FCGR1A+PDCD1LG2", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

# GBP5, FCGR1A
data$sig <- data$GBP5+data$FCGR1A
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "GBP5+FCGR1A", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

# GBP5, PDCD1LG2
data$sig <- data$GBP5+data$PDCD1LG2
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "GBP5+PDCD1LG2", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")