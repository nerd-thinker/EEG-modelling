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

# --- Leave-one-out cross-validation instead of a train/test split ---
# With only 49 participants, an 80/20 holdout leaves ~9 test rows — far too
# few to trust a single AUC, and repeated k-fold still pools multiple
# participants per fold. LOOCV holds out exactly one participant per fold,
# refits on everyone else, and predicts that one held-out participant — so
# every participant is used for both training and testing (never in the
# same fold), getting the most data-efficient estimate possible at this
# sample size.
x <- as.matrix(select(model_data, -`O/Y`))
y <- model_data$`O/Y`  # factor, levels "O", "Y"

ctrl <- trainControl(
  method = "LOOCV",           # n folds, each holding out exactly 1 participant
  classProbs = TRUE,
  summaryFunction = twoClassSummary,  # ROC / Sens / Spec, pooled across folds
  savePredictions = "final"
)

# alpha = 1 -> LASSO; search a range of lambdas instead of fixing one
grid <- expand.grid(alpha = 1, lambda = 10^seq(-3, 0, length.out = 30))
plot(grid$lambda)

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

# ROC / sensitivity / specificity at the best lambda, computed by pooling
# the one-held-out prediction from every LOOCV fold together — this is the
# number to report, not a single-split AUC
cv_model$results[which.max(cv_model$results$ROC), ]

# NOTE: unlike k-fold CV, each LOOCV fold holds out exactly 1 participant,
# so a per-fold ROC/AUC isn't computable (you need both classes present in
# a fold to compute ROC). There's no fold-to-fold "distribution" to inspect
# here the way there was with repeated 5-fold CV — cv_model$resample will
# have per-fold Sens/Spec but not a meaningful ROC spread. The pooled ROC
# above (from savePredictions = "final") is the metric to rely on.
cv_model$pred  # per-participant held-out class + predicted probability

# Features that survive shrinkage in the final model (refit on all
# participants at the CV-selected lambda) — treat as candidates to
# investigate further, not as confirmed effects
final_coefs <- coef(cv_model$finalModel, s = cv_model$bestTune$lambda)
final_coefs <- final_coefs[final_coefs[, 1] != 0, , drop = FALSE]
final_coefs

# Check values. The final coef might be ballooned due to coefficient effect.
# Also only 3 features???
summary(features_clean$theta_T8_none_mean)
summary(features_clean$theta_Fp2_none_mean)
summary(features_clean$alpha_Fp2_log10_mean)
# Because the scale is really small the ceoficient effect might have been enormours.
# But the scaling of the mode, which is default in this case should have taken care of that?

# Let's try only one scaling, z-score and train the model to see what features we are gonna get -------

library(dplyr)
library(stringr)

# All z_score feature columns, keeping O/Y (and cross_ columns, which aren't 
# scaling-variant to begin with)
z_score_data <- model_data %>%
  select(`O/Y`, matches("_z_score_"), starts_with("cross_"))

dim(z_score_data)
