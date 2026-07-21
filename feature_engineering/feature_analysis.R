library(readr)
library(tidyverse)
library(dtplyr)

features_all <- read_csv("~/eeg/features/ALL_PARTICIPANTS_features.csv")
nrow(features_all)
anyNA(features_all)
is.na(features_all)
sum(is.na(features_all))

na_summary <- colSums(is.na(features_all))
na_summary[na_summary > 0]


na_rows <- rowSums(is.na(features_all)) > 0
features_all$participant_id[na_rows]

saveRDS(na_summary, file = "na_summary.rds")
write.csv(na_summary, file = "~/eeg/na_summary.csv")

# Find the badly named columns - these are the broken ones
bad_cols <- names(na_summary)[grepl("^NA\\.\\.\\.", names(na_summary))]
length(bad_cols)  # how many

# Find which band_channel_scaling combos are fully broken
# (NA for most participants)
high_na <- na_df[na_df$na_count > 10, ]
head(high_na, 20)

# Extract the unique band_channel combos that are broken
broken_combos <- unique(
  gsub("_(none|z_score|01|log10)_.*$", "", high_na$feature)
)
print(broken_combos)