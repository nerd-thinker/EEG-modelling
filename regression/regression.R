library(dplyr)     # Data wrangling (select, mutate, pipe)
library(glmnet)    # LASSO logistic regression — cv.glmnet()
library(caret)     # Train/test splitting, cross-validation framework
library(pROC)      # AUC-ROC calculation and plotting

# Load features_clean.rds -> model_data
features_clean <- readRDS("~/eeg/features_clean.rds")

# Including O/Y predictor in the data, replacing participant ID
model_data <- features_clean %>%
  mutate(`O/Y` = factor(substr(participant_id, 9, 9)), .before = participant_id) %>%
  select(-participant_id)
## first col = O/Y predictor

# How many Y vs O participants there is
table(model_data$`O/Y`)

# Check is there is any Na's left
sum(is.na(model_data))
colSums(is.na(model_data))[colSums(is.na(model_data)) > 0]

# Train data on the same amount of old/young (approx 80% of the data)
# createDataPartition stratifies on O/Y so the split keeps the class balance
set.seed(42)
train_idx <- createDataPartition(model_data$`O/Y`, p = 0.8, list = FALSE)
train_data <- model_data[train_idx, ]
test_data  <- model_data[-train_idx, ]

table(train_data$`O/Y`)
table(test_data$`O/Y`)

# --- Repeated stratified k-fold CV instead of a single train/test split ---
# With only 49 participants, one 80/20 holdout leaves ~9 test rows — far too
# few to trust a single AUC. Repeated CV refits across many different fold
# assignments and averages the result, which is much more stable here, and
# is also less likely to be fooled by the near-perfect separation you get
# almost for free with 10k+ predictors and <50 samples.
x <- as.matrix(select(model_data, -`O/Y`))
y <- model_data$`O/Y`  # factor, levels "O", "Y"

ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,                # 5 folds, stratified on O/Y by default
  repeats = 20,               # repeat the whole 5-fold split 20x
  classProbs = TRUE,
  summaryFunction = twoClassSummary,  # ROC / Sens / Spec per fold
  savePredictions = "final"
)

# alpha = 1 -> LASSO; search a range of lambdas instead of fixing one
grid <- expand.grid(alpha = 1, lambda = 10^seq(-3, 0, length.out = 30))

set.seed(42)
cv_model <- train(
  x = x, y = y,
  method = "glmnet",
  family = "binomial",
  metric = "ROC",
  trControl = ctrl,
  tuneGrid = grid
)

cv_model
cv_model$bestTune

# Mean ROC / sensitivity / specificity across all folds x repeats, at the
# best lambda — this is the number to report, not a single-split AUC
cv_model$results[which.max(cv_model$results$ROC), ]

# Distribution of per-resample ROC (49 fitted folds x repeats) — check how
# much it varies fold to fold, not just the mean
cv_model$resample$ROC
summary(cv_model$resample$ROC)

# Features that survive shrinkage in the final model (refit on all 49 rows
# at the CV-selected lambda) — treat as candidates to investigate further,
# not as confirmed effects
final_coefs <- coef(cv_model$finalModel, s = cv_model$bestTune$lambda)
final_coefs <- final_coefs[final_coefs[, 1] != 0, , drop = FALSE]
final_coefs
