# batch_feature_extraction.R
# ─────────────────────────────────────────────────────────────────────────────
# Parallel batch feature extraction — one participant per core.
# Designed to be run from the terminal with:
#
#   R CMD BATCH --no-save --no-restore batch_feature_extraction.R batch.log &
#
# Progress is written to batch.log. To watch it live:
#   tail -f batch.log
#
# Uses 5 of 6 cores (leaves 1 free for the OS).
# Each participant saves its own CSV immediately — safe to resume after crash:
# set skip_done = TRUE (default) to skip already-completed participants.
# ─────────────────────────────────────────────────────────────────────────────

library(parallel)
library(dplyr)

source("~/eeg/feature_engineering/complete_feature_function.R")
source("~/eeg/filepaths_data_for_each_participant.R")

# ── Config ────────────────────────────────────────────────────────────────────
OUT_DIR   <- "~/eeg/features"
K         <- 70
SKIP_DONE <- TRUE
N_CORES   <- 5   # Ryzen 5 2600 has 6 cores — leave 1 free for the OS

# ── Setup ─────────────────────────────────────────────────────────────────────
out_dir <- path.expand(OUT_DIR)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

participant_ids <- names(all_filepaths)
n               <- length(participant_ids)

# Filter out already-done participants if skip_done = TRUE
if (SKIP_DONE) {
  already_done <- file.exists(
    file.path(out_dir, paste0(participant_ids, "_features.csv"))
  )
  if (any(already_done)) {
    message(sprintf("Skipping %d already-completed participants: %s",
                    sum(already_done),
                    paste(participant_ids[already_done], collapse = ", ")))
  }
  todo_ids <- participant_ids[!already_done]
} else {
  todo_ids <- participant_ids
}

message(sprintf("\n%d participants to process across %d cores\n", length(todo_ids), N_CORES))

if (length(todo_ids) == 0) {
  message("Nothing to do — all participants already complete.")
  quit(save = "no")
}

# ── Worker function ───────────────────────────────────────────────────────────
# Runs on each core — must be self-contained (sources its own dependencies)

run_one_participant <- function(pid, all_filepaths, out_dir, k) {
  
  # Each worker needs its own sourced environment
  source("~/eeg/feature_engineering/complete_feature_function.R")
  
  out_file <- file.path(out_dir, paste0(pid, "_features.csv"))
  t0       <- proc.time()
  
  result <- tryCatch({
    row <- extract_all_features(
      filepaths      = all_filepaths[[pid]],
      participant_id = pid,
      k              = k
    )
    write.csv(row, out_file, row.names = FALSE)
    list(ok = TRUE, elapsed_min = round((proc.time() - t0)["elapsed"] / 60, 1))
    
  }, error = function(e) {
    list(ok = FALSE, error = e$message,
         elapsed_min = round((proc.time() - t0)["elapsed"] / 60, 1))
  })
  
  result$pid <- pid
  result
}

# ── Run in parallel ───────────────────────────────────────────────────────────
batch_start <- proc.time()

cl <- makeCluster(min(N_CORES, length(todo_ids)), type = "FORK")
# FORK copies the current R session to each worker (Linux/Mac only).
# If you are on Windows, change to type = "PSOCK" and the worker function
# will source() everything it needs via the source() call inside run_one_participant().

results_list <- parLapply(cl, todo_ids, function(pid) {
  run_one_participant(pid, all_filepaths, out_dir, K)
})

stopCluster(cl)

# ── Summarise results ─────────────────────────────────────────────────────────
summary_df <- data.frame(
  participant_id = sapply(results_list, `[[`, "pid"),
  status         = ifelse(sapply(results_list, `[[`, "ok"), "done", "failed"),
  elapsed_min    = sapply(results_list, `[[`, "elapsed_min"),
  error          = sapply(results_list, function(x) if (!x$ok) x$error else ""),
  stringsAsFactors = FALSE
)

# Add skipped participants back into the summary
if (SKIP_DONE && any(already_done)) {
  skipped_df <- data.frame(
    participant_id = participant_ids[already_done],
    status         = "skipped",
    elapsed_min    = NA_real_,
    error          = "",
    stringsAsFactors = FALSE
  )
  summary_df <- bind_rows(skipped_df, summary_df)
}

total_min <- (proc.time() - batch_start)["elapsed"] / 60
message(sprintf("\n══ BATCH COMPLETE in %.1f min (%.1f hrs) ══", total_min, total_min / 60))
message("\nPer-participant summary:")
print(summary_df)

failed <- summary_df[summary_df$status == "failed", ]
if (nrow(failed) > 0) {
  message("\nFailed participants:")
  print(failed)
}

# ── Combine all CSVs into master file ─────────────────────────────────────────
message("\nCombining participant CSVs into master file...")

csv_files <- list.files(out_dir, pattern = "_features\\.csv$", full.names = TRUE)
csv_files <- csv_files[!grepl("ALL_PARTICIPANTS", csv_files)]

if (length(csv_files) == 0) {
  warning("No feature CSVs found — nothing to combine.")
} else {
  combined    <- bind_rows(lapply(csv_files, read.csv, check.names = FALSE))
  master_path <- file.path(out_dir, "ALL_PARTICIPANTS_features.csv")
  write.csv(combined, master_path, row.names = FALSE)
  message(sprintf("Master file saved: %s", master_path))
  message(sprintf("  %d participants x %d feature columns",
                  nrow(combined), ncol(combined) - 1))
}

message("\nDone.")