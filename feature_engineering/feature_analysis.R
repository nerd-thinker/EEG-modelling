library(readr)
library(tidyverse)
library(dtplyr)

features_all <- read_csv("~/eeg/features/ALL_PARTICIPANTS_features.csv")
feat001Y <- read_csv("~/eeg/features/001Y_features.csv")
nrow(features_all)
anyNA(features_all)
is.na(features_all)
sum(is.na(feat001Y))

na_summary <- colSums(is.na(features_all))
na_summary[na_summary > 0]

na_summary1 <- colSums(is.na(feat001Y))
na_summary1[na_summary1 > 0]

na_rows <- rowSums(is.na(feat001Y)) > 0
features_all$participant_id[na_rows]

saveRDS(na_summary, file = "na_summary.rds")
write.csv(na_summary, file = "~/eeg/na_summary.csv")

