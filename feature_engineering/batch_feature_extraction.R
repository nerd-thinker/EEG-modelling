# batch_feature_extraction.R
# ─────────────────────────────────────────────────────────────────────────────
# Runs extract_all_features() for every participant in all_filepaths,
# saves each result to disk immediately, and combines into one master CSV.
#
# REQUIRES:
#   - complete_feature_function.R sourced (provides extract_all_features)
#   - filepaths_data_for_each_participant.R sourced (provides all_filepaths)
#
# USAGE:
#   source("batch_feature_extraction.R")   # defines the function
#   results <- run_batch(all_filepaths, out_dir = "~/eeg/features", k = 70)
# ─────────────────────────────────────────────────────────────────────────────

run_batch <- function(all_filepaths,
                      out_dir    = "~/eeg/features",
                      k          = 70,
                      skip_done  = TRUE) {
  # all_filepaths : named list of filepath vectors, one entry per participant
  # out_dir       : folder where per-participant CSVs and the master CSV are saved
  # k             : GAM basis dimension passed to extract_all_features()
  # skip_done     : if TRUE, skips participants whose CSV already exists in out_dir
  #                 (lets you resume a crashed batch without re-running anyone)
  
  # ── Setup ──────────────────────────────────────────────────────────────────
  out_dir <- path.expand(out_dir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  participant_ids <- names(all_filepaths)
  n               <- length(participant_ids)
  
  # Tracking vectors
  status      <- rep(NA_character_, n)   # "done" | "skipped" | "failed"
  elapsed_min <- rep(NA_real_,      n)
  names(status) <- names(elapsed_min) <- participant_ids
  
  batch_start <- proc.time()
  
  # ── Loop ───────────────────────────────────────────────────────────────────
  for (i in seq_along(participant_ids)) {
    
    pid      <- participant_ids[i]
    out_file <- file.path(out_dir, paste0(pid, "_features.csv"))
    
    # ── Skip if already done ────────────────────────────────────────────────
    if (skip_done && file.exists(out_file)) {
      message(sprintf("[%d/%d] SKIP  %s — file already exists", i, n, pid))
      status[pid] <- "skipped"
      next
    }
    
    message(sprintf("\n[%d/%d] ── START %s ── %s",
                    i, n, pid, format(Sys.time(), "%H:%M:%S")))
    
    t0 <- proc.time()
    
    # ── Run feature extraction with full error catching ──────────────────────
    result <- tryCatch({
      
      row <- extract_all_features(
        filepaths      = all_filepaths[[pid]],
        participant_id = pid,
        k              = k
      )
      
      # Save immediately — if the batch crashes later this participant is safe
      write.csv(row, out_file, row.names = FALSE)
      message(sprintf("  Saved → %s", out_file))
      
      list(ok = TRUE, data = row)
      
    }, error = function(e) {
      message(sprintf("  !! FAILED: %s", e$message))
      list(ok = FALSE, data = NULL)
    })
    
    elapsed         <- (proc.time() - t0)["elapsed"] / 60
    elapsed_min[pid] <- round(elapsed, 1)
    
    if (result$ok) {
      status[pid] <- "done"
      message(sprintf("  Done in %.1f min", elapsed))
    } else {
      status[pid] <- "failed"
      message(sprintf("  Failed after %.1f min", elapsed))
    }
    
    # ── Remaining time estimate ──────────────────────────────────────────────
    done_so_far <- sum(status %in% c("done", "skipped"), na.rm = TRUE)
    if (done_so_far > 0) {
      total_elapsed <- (proc.time() - batch_start)["elapsed"] / 60
      avg_per_p     <- total_elapsed / done_so_far
      remaining     <- avg_per_p * (n - done_so_far)
      message(sprintf("  Batch progress: %d/%d | ~%.0f min remaining",
                      done_so_far, n, remaining))
    }
  }
  
  # ── Combine all saved CSVs into one master file ───────────────────────────
  message("\n── Combining all participant CSVs into master file... ──")
  
  csv_files <- list.files(out_dir,
                          pattern    = "_features\\.csv$",
                          full.names = TRUE)
  
  # Exclude the master file itself if it exists from a previous run
  csv_files <- csv_files[!grepl("ALL_PARTICIPANTS", csv_files)]
  
  if (length(csv_files) == 0) {
    warning("No feature CSVs found in ", out_dir, " — nothing to combine.")
    master <- NULL
  } else {
    master <- tryCatch({
      all_rows <- lapply(csv_files, read.csv, check.names = FALSE)
      combined <- dplyr::bind_rows(all_rows)   # fills missing columns with NA
      
      master_path <- file.path(out_dir, "ALL_PARTICIPANTS_features.csv")
      write.csv(combined, master_path, row.names = FALSE)
      message(sprintf("Master file saved → %s", master_path))
      message(sprintf("  %d participants × %d feature columns",
                      nrow(combined), ncol(combined) - 1))
      combined
    }, error = function(e) {
      message("  !! Could not combine CSVs: ", e$message)
      NULL
    })
  }
  
  # ── Summary ───────────────────────────────────────────────────────────────
  total_min <- (proc.time() - batch_start)["elapsed"] / 60
  message(sprintf("\n══ BATCH COMPLETE in %.1f min (%.1f hrs) ══", total_min, total_min / 60))
  
  summary_df <- data.frame(
    participant_id = participant_ids,
    status         = status,
    elapsed_min    = elapsed_min,
    stringsAsFactors = FALSE
  )
  
  message("\nPer-participant summary:")
  print(summary_df)
  
  # ── Return ────────────────────────────────────────────────────────────────
  list(
    master  = master,       # combined data frame (NULL if combining failed)
    summary = summary_df    # who ran, who was skipped, who failed
  )
}


# ── Run ───────────────────────────────────────────────────────────────────────
# Source your other scripts first, then call run_batch():

source("~/eeg/complete_feature_function.R")
source("~/eeg/filepaths_data_for_each_participant.R")

results <- run_batch(
  all_filepaths = all_filepaths,
  out_dir       = "~/eeg/features",
  k             = 70,
  skip_done     = TRUE     # change to FALSE to re-run everyone from scratch
)

# Access outputs:
results$master     # full participant × feature data frame
results$summary    # per-participant status table

# Check who failed:
results$summary[results$summary$status == "failed", ]

# Check who was slowest:
results$summary[order(-results$summary$elapsed_min), ]