library(DMwR)
library(smotefamily)

# Read the data
data <- read.csv2("gene_expression_dct_.csv", sep = ",") # gene_expression_dct_noIRAK3.csv

table(data$diag)

# Create signature score
data$sig <- data$FCGR1A + data$GBP5

# Select relevant columns
data <- data %>%
  select(diag, sig)

# Convert diagnosis to factor with specific levels
data$diag <- factor(data$diag, levels = c("TB", "NO_TB"), ordered = TRUE)

# Create training/testing split
set.seed(123)  # Set seed before createDataPartition
trainIndex <- createDataPartition(data$diag, p = 0.8, list = FALSE)
trainData <- data[trainIndex,]
testData <- data[-trainIndex,]

trainData <- SMOTE(diag ~ ., data = trainData, perc.over = 100, perc.under = 200)


# Set up cross-validation
fitControl <- trainControl(method = "cv", 
                           number = 5, 
                           savePredictions = 'final',
                           classProbs = TRUE, 
                           summaryFunction = twoClassSummary)

# Train SVM model
set.seed(123)  # Set seed before training
svm_reg <- train(diag ~ ., 
                 data = trainData, 
                 method = "svmRadial", 
                 trControl = fitControl, 
                 metric = "ROC")


svm_reg <- train(diag ~ ., data = trainData, method = "glm", trControl = fitControl,
                 metric = "ROC")



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
svm_reg <- train(diag ~ ., data = trainData, method = "xgbTree", trControl = fitControl,
                   tuneGrid = xgbGrid, metric = "ROC")

# Make predictions
predictions <- predict(svm_reg, newdata = testData)
prob_predictions <- predict(svm_reg, newdata = testData, type = "prob")

# Calculate ROC
roc_obj <- roc(testData$diag, prob_predictions$TB)
auc_value <- auc(roc_obj)

# Calculate confusion matrix
confMatrix <- confusionMatrix(predictions, testData$diag, mode = "everything")
print(confMatrix)

cal_obj <- calibration(testData$diag ~ prob_predictions$TB, class = "TB")
plot(cal_obj)

# FCGR1A, PDCD1LG2
data$sig <- data$FCGR1A+data$PDCD1LG2
roc_curve <- roc(testData$diag, prob_predictions$TB) 
auc_value <- auc(roc_curve)
ci_value <- ci.auc(roc_curve)
plot(roc_curve, main = "FCGR1A+PDCD1LG2+GPB5", print.thres = "best", col = "blue", legacy.axes = TRUE) 
legend("bottom", legend = c(paste("AUC:", round(auc_value, 2)),paste("95% CI:", paste(round(ci_value[1], 3), "-", round(ci_value[2], 3)))), bty = "n")
