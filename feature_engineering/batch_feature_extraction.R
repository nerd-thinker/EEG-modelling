# batch_feature_extraction.R
# ─────────────────────────────────────────────────────────────────────────────
# Can be used two ways:
#
# 1. From the terminal (background, logs to file):
#      R CMD BATCH --no-save --no-restore batch_feature_extraction.R batch.log &
#      tail -f batch.log
#
# 2. From the R console interactively:
#      source("~/eeg/batch_feature_extraction.R")
#      run_batch_parallel(all_filepaths)
#
# Safe to interrupt at any time — completed CSVs are already on disk.
# Resume by re-running with skip_done = TRUE (default).
# ─────────────────────────────────────────────────────────────────────────────

library(parallel)
library(dplyr)

# ── Worker function ───────────────────────────────────────────────────────────
# Self-contained: sources its own dependencies so it works inside FORK workers

.run_one_participant <- function(pid, all_filepaths, out_dir, k) {
  
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

# ── Progress bar ──────────────────────────────────────────────────────────────
.print_progress <- function(completed, n_todo, n_done, n_failed,
                            batch_start, elapsed_per_batch, n_cores) {
  
  pct       <- completed / n_todo
  bar_width <- 30
  filled    <- round(pct * bar_width)
  bar       <- paste0(strrep("=", max(0, filled - 1)),
                      if (filled > 0 && filled < bar_width) ">" else "",
                      strrep(" ", bar_width - filled))
  
  elapsed_min  <- (proc.time() - batch_start)["elapsed"] / 60
  batches_left <- ceiling((n_todo - completed) / n_cores)
  
  if (length(elapsed_per_batch) > 0) {
    avg_batch_min <- mean(elapsed_per_batch)
    eta_min       <- batches_left * avg_batch_min
    eta_clock     <- format(Sys.time() + eta_min * 60, "%H:%M")
    eta_str       <- sprintf("ETA: %s (~%.0f min)", eta_clock, eta_min)
  } else {
    eta_str <- "ETA: calculating..."
  }
  
  message(sprintf("\n[%s] %d/%d | done: %d  failed: %d | elapsed: %.1f min | %s\n",
                  bar, completed, n_todo, n_done, n_failed, elapsed_min, eta_str))
}

# ── Main batch function ───────────────────────────────────────────────────────

run_batch_parallel <- function(all_filepaths,
                               out_dir   = "~/eeg/features",
                               k         = 70,
                               skip_done = TRUE,
                               n_cores   = 5) {
  
  out_dir <- path.expand(out_dir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  participant_ids <- names(all_filepaths)
  n               <- length(participant_ids)
  
  # ── Skip already-done participants ─────────────────────────────────────────
  if (skip_done) {
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
    return(invisible(NULL))
  }
  
  message(sprintf("\n%d participant(s) to process across %d cores", n_todo, n_cores))
  message(sprintf("Estimated total time: ~%.0f hrs (assuming ~115 min per participant)\n",
                  ceiling(n_todo / n_cores) * 115 / 60))
  
  # ── Split into batches of n_cores ──────────────────────────────────────────
  batches           <- split(todo_ids, ceiling(seq_along(todo_ids) / n_cores))
  batch_start       <- proc.time()
  all_results       <- list()
  elapsed_per_batch <- numeric(0)
  n_done            <- 0
  n_failed          <- 0
  completed         <- 0
  
  for (batch_num in seq_along(batches)) {
    
    batch      <- batches[[batch_num]]
    n_in_batch <- length(batch)
    
    message(sprintf("── Batch %d/%d: processing %s",
                    batch_num, length(batches),
                    paste(batch, collapse = ", ")))
    
    batch_t0 <- proc.time()
    
    cl <- makeCluster(min(n_cores, n_in_batch), type = "FORK")
    batch_results <- parLapply(cl, batch, function(pid) {
      .run_one_participant(pid, all_filepaths, out_dir, k)
    })
    stopCluster(cl)
    
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
    elapsed_per_batch <- c(elapsed_per_batch,
                           (proc.time() - batch_t0)["elapsed"] / 60)
    
    .print_progress(completed, n_todo, n_done, n_failed,
                    batch_start, elapsed_per_batch, n_cores)
  }
  
  # ── Final summary ──────────────────────────────────────────────────────────
  total_min <- (proc.time() - batch_start)["elapsed"] / 60
  message(sprintf("══ BATCH COMPLETE in %.1f min (%.1f hrs) ══\n",
                  total_min, total_min / 60))
  
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
  
  message("Per-participant summary:")
  print(summary_df)
  
  failed <- summary_df[summary_df$status == "failed", ]
  if (nrow(failed) > 0) {
    message("\nFailed participants:")
    print(failed[, c("participant_id", "error")])
  }
  
  # ── Combine all CSVs into master file ──────────────────────────────────────
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
  invisible(summary_df)
}


# ── Execution ─────────────────────────────────────────────────────────────────
# This block runs automatically under R CMD BATCH.
# When sourcing interactively, skip this block and call run_batch_parallel()
# manually after sourcing your other scripts.

if (!interactive()) {
  source("~/eeg/feature_engineering/complete_feature_function.R")
  source("~/eeg/filepaths_data_for_each_participant.R")
  
  run_batch_parallel(
    all_filepaths = all_filepaths,
    out_dir       = "~/eeg/features",
    k             = 70,
    skip_done     = TRUE,
    n_cores       = 5
  )
}