# Load the required packages
library(tidyverse)
library(caret)
library(randomForest)
library(xgboost)
library(pROC)


set.seed(1234)
# Load the data
data <- read.csv2("~/Documentos/Thaina/gene_expression_dct_FCGR1A_only.csv", sep = ",") 
View(data)

# gene_expression_dct_noIRAK3.csv

#data$diag <- as.factor(data$diag)
data$diag <- factor(data$diag, levels = c("NO_TB", "TB"), ordered = TRUE)

# Split the data into training and testing sets
trainIndex <- createDataPartition(data$diag, p = 0.8, list = FALSE)
trainData <- data[trainIndex,]
testData <- data[-trainIndex,]


# Cross-validation ----
fitControl <- trainControl(method = "cv",  
                           number = 5,      
                           savePredictions = 'final',
                           classProbs = TRUE,
                           summaryFunction=twoClassSummary)

# Model Analysis ----
set.seed(123)
#rf <- randomForest(diag~., data = trainData, tuneLength = 10, trControl = fitControl)
rf <- train(diag~., data = trainData, method="rf", tuneLength = 10, metric='ROC', trControl = fitControl)

p1 <- predict(rf, trainData)
confusionMatrix(p1, trainData$diag)

p2 <- predict(rf, testData)
confusionMatrix(p2, testData$diag)

model <- varImp(rf)
plot(model)

# AUC visualization ----
rocobj1 <- plot.roc(data$diag, data$sig_mean,
                    main="Statistical comparison",
                    percent=TRUE,
                    col="#1c61b6")

rocobj2 <- lines.roc(data$diag, data$sig_weighted, 
                     percent=TRUE, 
                     col="#008600")

testobj <- roc.test(rocobj1, rocobj2)

# PDCD1LG2
data$sig <- data$PDCD1LG2
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "PDCD1LG2", print.thres = "best", col = "blue", legacy.axes = TRUE)
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

# GBP5
data$sig <- data$GBP5
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "GBP5", print.thres = "best", col = "blue", legacy.axes = TRUE)
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

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

# IRAK3, PDCD1LG2
data$sig <- data$IRAK3+data$PDCD1LG2
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "IRAK3+PDCD1LG2", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 2), "-", round(ci_value[2], 2)))), bty = "n")

# FCGR1A, PDCD1LG2
data$sig <- data$FCGR1A+data$PDCD1LG2
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "FCGR1A+PDCD1LG2", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

# FCGR1A, PDCD1LG2, GBP5, IRAK3
data$sig <- data$FCGR1A+data$PDCD1LG2+data$GBP5+data$IRAK3
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "FCGR1A+PDCD1LG2+GBP5+IRAK3", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

# GBP5, FCGR1A
data$sig <- data$GBP5+data$FCGR1A
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "GBP5+FCGR1A", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

# GBP5, IRAK3, PDCD1LG2
data$sig <- data$GBP5+data$IRAK3+data$PDCD1LG2
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "GBP5+IRAK3+PDCD1LG2", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

# GBP5, PDCD1LG2
data$sig <- data$GBP5+data$PDCD1LG2
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "GBP5+PDCD1LG2", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

# IRAK3, FCGR1A
data$sig <- data$IRAK3+data$FCGR1A
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "IRAK3+FCGR1A", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

# IRAK3, GBP5 
data$sig <- data$IRAK3+data$GBP5
roc_curve <- roc(diag ~ sig, data)
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "IRAK3+GBP5", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")

phase3 <- read.csv2("~/Documentos/Thaina/gene_expression_phase3_FCGR1A_only.csv", sep = ",")

predictions <- predict(rf, newdata = phase3)

phase3$predicted_class <- predictions
View(phase3)

write.csv(phase3, "FCGR1A/phase3_predictions_rf_FCGR1A.csv")