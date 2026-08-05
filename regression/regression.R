library(dplyr)
library(glmnet)

# Including O/Y predictor in the data, replacing participant ID
model_data <- features_clean %>%
  mutate('O/Y' = factor(substr(participant_id, 9, 9)), .before = participant_id, type = "factor") %>%
  select(-participant_id)
## first col = O/Y predictor
