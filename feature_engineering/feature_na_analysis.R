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

na_summary1 <- colSums(is.na(features_all))
na_summary1[na_summary1 > 1]

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
table(na_summary1) 
## plot to see if some number of feature Na repeat to detect a possible pattern
plot(table(na_summary1)) #, xlab = "feature Na's frequency", ylab = "Amount of the frequency occurrence"

# Check: are all 2+ NA features peak-related?
features_2plus_names <- names(na_summary1)
all_peak <- all(grepl("peak", features_2plus_names))
cat("All 2,304 features are peak-related:", all_peak, "\n")

# Count peak vs non-peak
num_peak <- sum(grepl("peak", features_2plus_names))
num_non_peak <- length(features_2plus_names) - num_peak
cat("Peak-related:", num_peak, "\n")
cat("Non-peak-related:", num_non_peak, "\n")

# Deleting all features that contain Na ----------------
features_with_any_na <- names(na_summary[na_summary > 0])

## keep features without any Na's
features_clean <- x[, !colnames(x) %in% features_with_any_na]

## verify cleaning worked
cat("Features deleted:", length(features_with_any_na), "\n")
cat("Features remaining:", ncol(features_clean) - 1, "\n")
cat("Participants:", nrow(features_clean), "\n")
cat("NAs remaining:", sum(is.na(features_clean)), "\n")

saveRDS(features_clean, file = "features_clean.rds")
