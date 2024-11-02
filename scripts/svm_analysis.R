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

set.seed(123)
svm_reg <- train(diag ~ ., data = trainData, method = "svmRadial", trControl = fitControl, metric = "ROC")


# Make predictions on test set
predictions <- predict(svm_reg, newdata = testData)
prob_predictions <- predict(svm_reg, newdata = testData, type = "prob")

# Calculate confusion matrix
confMatrix <- confusionMatrix(predictions, testData$diag)
print(confMatrix)

# Feature importance
importance <- varImp(svm_reg)
plot(importance)


phase3 <- read.csv2("~/Documentos/Thaina/gene_expression_phase3_FCGR1A_only.csv", sep = ",")

predictions <- predict(svm_reg, newdata = phase3)

phase3$predicted_class <- predictions
View(phase3)

write.csv(phase3, "FCGR1A/phase3_predictions_svm_FCGR1A.csv")


roc_plot <- plot_roc_curve(
  actual = data$diag,      # True labels (0/1)
  predicted = phase3$predicted_class,    # Predicted probabilities
  title = "Your Model ROC Curve"
)

plot.roc(data$diag, data$sig)

# FCGR1A
data$sig <- data$FCGR1A
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "FCGR1A", print.thres = "best", col = "blue", legacy.axes = TRUE)
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

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