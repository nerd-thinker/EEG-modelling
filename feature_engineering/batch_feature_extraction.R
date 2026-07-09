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
# Safe to interrupt at any time with Ctrl+C — completed participant CSVs
# are already saved to disk. Re-run with SKIP_DONE = TRUE to resume.
# ─────────────────────────────────────────────────────────────────────────────

library(parallel)
library(dplyr)

source("~/eeg/feature_engineering/complete_feature_function.R")
source("~/eeg/raw_data/filepaths_data_for_each_participant.R")

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

# Filter already-done participants
if (SKIP_DONE) {
  already_done <- file.exists(
    file.path(out_dir, paste0(participant_ids, "_features.csv"))
  )
  if (any(already_done)) {
    message(sprintf("Skipping %d already-completed participant(s): %s",
                    sum(already_done),
                    paste(participant_ids[already_done], collapse = ", ")))
  }
  todo_ids <- participant_ids[!already_done]
} else {
  already_done <- rep(FALSE, n)
  todo_ids     <- participant_ids
}

n_todo <- length(todo_ids)

if (n_todo == 0) {
  message("Nothing to do — all participants already complete.")
  quit(save = "no")
}

message(sprintf("\n%d participant(s) to process across %d cores", n_todo, N_CORES))
message(sprintf("Estimated time: ~%.0f hrs (assuming ~115 min per participant)\n",
                ceiling(n_todo / N_CORES) * 115 / 60))

# ── Worker function ───────────────────────────────────────────────────────────
run_one_participant <- function(pid, all_filepaths, out_dir, k) {
  
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
    list(ok          = TRUE,
         elapsed_min = round((proc.time() - t0)["elapsed"] / 60, 1),
         error       = "")
  }, error = function(e) {
    list(ok          = FALSE,
         elapsed_min = round((proc.time() - t0)["elapsed"] / 60, 1),
         error       = e$message)
  })
  
  result$pid <- pid
  result
}

# ── Progress bar helpers ───────────────────────────────────────────────────────
# Prints a text progress bar to the log that looks like:
#   [=====>    ] 6/19 | done: 5 failed: 1 | elapsed: 12.3 min | ETA: 14:42

format_eta <- function(eta_time) {
  format(as.POSIXct(eta_time, origin = "1970-01-01"), "%H:%M")
}

print_progress <- function(completed, n_todo, n_done, n_failed,
                           batch_start, elapsed_per_batch) {
  
  pct        <- completed / n_todo
  bar_width  <- 30
  filled     <- round(pct * bar_width)
  bar        <- paste0(strrep("=", max(0, filled - 1)),
                       if (filled > 0 && filled < bar_width) ">" else "",
                       strrep(" ", bar_width - filled))
  
  elapsed_min  <- (proc.time() - batch_start)["elapsed"] / 60
  batches_left <- ceiling((n_todo - completed) / N_CORES)
  
  if (length(elapsed_per_batch) > 0) {
    avg_batch_min <- mean(elapsed_per_batch)
    eta_min       <- batches_left * avg_batch_min
    eta_clock     <- format_eta(as.numeric(Sys.time()) + eta_min * 60)
    eta_str       <- sprintf("ETA: %s (~%.0f min)", eta_clock, eta_min)
  } else {
    eta_str <- "ETA: calculating..."
  }
  
  message(sprintf("\n[%s] %d/%d | done: %d  failed: %d | elapsed: %.1f min | %s\n",
                  bar, completed, n_todo, n_done, n_failed, elapsed_min, eta_str))
}

# ── Run in parallel batches ───────────────────────────────────────────────────
# Split todo_ids into batches of N_CORES. After each batch completes,
# update the progress bar before starting the next batch.
# This is what makes progress reporting possible with parallel workers.

batches <- split(todo_ids,
                 ceiling(seq_along(todo_ids) / N_CORES))

batch_start      <- proc.time()
all_results      <- list()
elapsed_per_batch <- numeric(0)
n_done           <- 0
n_failed         <- 0
completed        <- 0

for (batch_num in seq_along(batches)) {
  
  batch     <- batches[[batch_num]]
  n_in_batch <- length(batch)
  
  message(sprintf("── Batch %d/%d: processing %s",
                  batch_num, length(batches),
                  paste(batch, collapse = ", ")))
  
  batch_t0 <- proc.time()
  
  cl <- makeCluster(min(N_CORES, n_in_batch), type = "FORK")
  batch_results <- parLapply(cl, batch, function(pid) {
    run_one_participant(pid, all_filepaths, out_dir, K)
  })
  stopCluster(cl)
  
  # Record results
  for (res in batch_results) {
    all_results[[res$pid]] <- res
    if (res$ok) {
      n_done <- n_done + 1
      message(sprintf("  OK      %-15s  %.1f min", res$pid, res$elapsed_min))
    } else {
      n_failed <- n_failed + 1
      message(sprintf("  FAILED  %-15s  %.1f min  -- %s",
                      res$pid, res$elapsed_min, res$error))
    }
  }
  
  completed         <- completed + n_in_batch
  batch_elapsed_min <- (proc.time() - batch_t0)["elapsed"] / 60
  elapsed_per_batch <- c(elapsed_per_batch, batch_elapsed_min)
  
  print_progress(completed, n_todo, n_done, n_failed,
                 batch_start, elapsed_per_batch)
}

# ── Final summary ─────────────────────────────────────────────────────────────
total_min <- (proc.time() - batch_start)["elapsed"] / 60

message(sprintf("══ BATCH COMPLETE in %.1f min (%.1f hrs) ══", total_min, total_min / 60))

summary_df <- data.frame(
  participant_id = c(participant_ids[already_done],
                     sapply(all_results, `[[`, "pid")),
  status         = c(rep("skipped", sum(already_done)),
                     ifelse(sapply(all_results, `[[`, "ok"), "done", "failed")),
  elapsed_min    = c(rep(NA_real_, sum(already_done)),
                     sapply(all_results, `[[`, "elapsed_min")),
  error          = c(rep("", sum(already_done)),
                     sapply(all_results, `[[`, "error")),
  stringsAsFactors = FALSE
)

message("\nPer-participant summary:")
print(summary_df)

failed <- summary_df[summary_df$status == "failed", ]
if (nrow(failed) > 0) {
  message("\nFailed participants:")
  print(failed[, c("participant_id", "error")])
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