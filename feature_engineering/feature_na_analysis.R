library(readr)
library(tidyverse)
library(dtplyr)
# intial surface level analysis -----------------
features_all <- read_csv("~/eeg/features/ALL_PARTICIPANTS_features.csv")

feat001Y <- read_csv("~/eeg/features/001Y_features.csv")
nrow(features_all)
anyNA(features_all)
is.na(features_all) 
sum(is.na(features_all))

na_summary <- colSums(is.na(features_all))
na_summary[na_summary > 0]

na_summary1 <- colSums(is.na(feat001Y))
na_summary1[na_summary1 > 0]

na_rows <- rowSums(is.na(feat001Y)) > 0
features_all$participant_id[na_rows]

saveRDS(na_summary, file = "na_summary.rds")
write.csv(na_summary, file = "~/eeg/na_summary.csv")

sum(apply(features_all, 2, function(x) any(is.na(x))))
ncol(features_all)
apply(features_all, 1, function(x) sum(any(is.na(x))))

# Na analysis -----------------
x <- read_csv("~/eeg/features/ALL_PARTICIPANTS_features.csv")
# to showcase what the following code found use the original dataset as the above dataset 
# has been manipulated to remove the NA features

# x <- read_csv("~/eeg/features/ALL_PARTICIPANTS_features_original.csv")
(ncol(x) -1) *  nrow(x)
anyNA(x)
sum(is.na(x))

na_summary <- colSums(is.na(x))
na_summary1 <- na_summary[na_summary > 1]
length(na_summary[na_summary1])

## rows that only have 1 Na: -----------------------
only1_na <- na_summary[na_summary == 1] ## only 394 of mean/max_peak_amplitude 
## only1_na output is a numerical vector
only1_na_names <- names(only1_na)

bands <- str_split_i(only1_na_names, "_", 1)
nodes <- str_split_i(only1_na_names, "_", 2)

features_with_1na <- data.frame(
  feature = only1_na_names,
  band = bands,
  node = nodes,
  stringsAsFactors = FALSE
)

table(features_with_1na$band)    # How many 1-NA cases per band?
table(features_with_1na$node)

# How many unique participants are affected by these 394 1-NA features?
affected_participants <- unique(which(is.na(x[, only1_na_names]), arr.ind = TRUE)[, "row"])
cat("Participants affected:", length(affected_participants), "\n")

## conclusion: as the only affected features are mean/max peak_amplitude, and only 394 of them are affected, we can safely 
## remove these features from the dataset without losing much information.

## more than one NA features: -----------------------
table(summary1)
