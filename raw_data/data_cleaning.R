# ============================================================================
# EEG DATA VALIDATOR / NORMALIZER
# ----------------------------------------------------------------------------
# Brings inconsistently-shaped EEG CSVs into the form complete_feature_function
# needs: 32 rows (one per node, correctly named) x time-point columns, with a
# strictly positive time vector. Detects orientation automatically rather than
# assuming a fixed layout, since your files come in several different forms.
# ============================================================================

# ---- 0. CONFIG: edit this to match your actual 32-channel montage ----------
# If you leave this NULL, the function falls back to "whichever axis has the
# most non-numeric, non-time-like labels" instead of checking exact names.
expected_nodes <- c(
  "Fp1","Fp2","F7","F3","Fz","F4","F8",
  "FC5","FC1","FC2","FC6",
  "T7","C3","Cz","C4","T8",
  "CP5","CP1","CP2","CP6",
  "P7","P3","Pz","P4","P8",
  "PO9","O1","Oz","O2","PO10",
  "AF3","AF4"
)  # <-- replace with your real 32 node names if different

time_keywords <- c("time", "t", "timestamp", "sec", "secs", "seconds", "ms")

# ---- 1. Helpers --------------------------------------------------------

.is_numeric_like <- function(x) {
  suppressWarnings(!any(is.na(as.numeric(x))))
}

.read_raw <- function(path) {
  read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

# ---- 2. Core: validate + fix a single file ------------------------------

#' @param path Path to one CSV file
#' @param expected_nodes Character vector of valid node names (or NULL to infer)
#' @param sampling_rate Hz, used ONLY if no time info exists at all in the file
#' @param positive_offset_policy "shift" (default) shifts the whole time vector
#'   up so the minimum becomes positive; "error" leaves it and just flags it invalid
check_and_fix_eeg <- function(path,
                              expected_nodes = NULL,
                              sampling_rate = NULL,
                              positive_offset_policy = c("shift", "error")) {
  
  positive_offset_policy <- match.arg(positive_offset_policy)
  file_name <- basename(path)
  notes <- character(0)
  
  out <- tryCatch({
    
    raw <- .read_raw(path)
    header_labels <- names(raw)
    first_col <- raw[[1]]
    
    # Decide if column 1 holds row labels (e.g. "Fp1","Time") rather than data
    first_col_is_labels <- !.is_numeric_like(first_col)
    
    if (first_col_is_labels) {
      rn  <- as.character(first_col)
      mat <- as.matrix(raw[, -1, drop = FALSE])
      rownames(mat) <- rn
      colnames(mat) <- header_labels[-1]
      mat <- apply(mat, 2, as.numeric)
      rownames(mat) <- rn
    } else {
      mat <- as.matrix(raw)
      colnames(mat) <- header_labels
      rownames(mat) <- NULL
    }
    storage.mode(mat) <- "numeric"
    
    # --- Orientation check: which axis actually holds the electrode names? ---
    rn_lc <- tolower(rownames(mat)); cn_lc <- tolower(colnames(mat))
    
    if (!is.null(expected_nodes)) {
      exp_lc <- tolower(expected_nodes)
      rn_match <- sum(rn_lc %in% exp_lc)
      cn_match <- sum(cn_lc %in% exp_lc)
    } else {
      # generic fallback: count non-numeric, non-time-keyword labels on each axis
      rn_match <- sum(!.is_numeric_like(rownames(mat)) & !(rn_lc %in% time_keywords))
      cn_match <- sum(!.is_numeric_like(colnames(mat)) & !(cn_lc %in% time_keywords))
    }
    
    was_transposed <- FALSE
    if (cn_match > rn_match) {
      mat <- t(mat)
      was_transposed <- TRUE
      rn_lc <- tolower(rownames(mat)); cn_lc <- tolower(colnames(mat))
      notes <- c(notes, "transposed: electrodes were in columns")
    }
    
    # --- Pull a "Time" row out of the data block if one snuck in as a node ---
    had_time_row <- any(rn_lc %in% time_keywords)
    time_vec <- NULL
    
    if (had_time_row) {
      time_row_idx <- which(rn_lc %in% time_keywords)[1]
      time_vec <- as.numeric(mat[time_row_idx, ])
      mat <- mat[-time_row_idx, , drop = FALSE]
      rn_lc <- rn_lc[-time_row_idx]
      notes <- c(notes, "extracted a 'Time' row out of the node block")
    } else if (.is_numeric_like(colnames(mat))) {
      # columns are already labelled with real time values
      time_vec <- as.numeric(colnames(mat))
    } else {
      # no time info anywhere - fall back to sampling_rate or sample index
      if (!is.null(sampling_rate)) {
        time_vec <- seq_len(ncol(mat)) / sampling_rate
        notes <- c(notes, sprintf("no time info found; generated from sampling_rate = %g Hz", sampling_rate))
      } else {
        time_vec <- seq_len(ncol(mat))
        notes <- c(notes, "no time info found and no sampling_rate given; used sample index as placeholder")
      }
    }
    
    # --- Sort columns chronologically (in case file wasn't already ordered) ---
    ord <- order(time_vec)
    mat <- mat[, ord, drop = FALSE]
    time_vec <- time_vec[ord]
    colnames(mat) <- as.character(time_vec)
    
    # --- Validate / reorder node rows against expected_nodes -----------------
    found_names <- rownames(mat)
    n_target <- if (!is.null(expected_nodes)) length(expected_nodes) else 32
    
    if (!is.null(expected_nodes)) {
      match_idx <- match(tolower(expected_nodes), tolower(found_names))
      missing_nodes <- expected_nodes[is.na(match_idx)]
      mat <- mat[match_idx[!is.na(match_idx)], , drop = FALSE]
      rownames(mat) <- expected_nodes[!is.na(match_idx)]   # canonical casing
    } else {
      missing_nodes <- character(0)
      if (nrow(mat) > 32) notes <- c(notes, "more than 32 candidate node rows found; left as-is, check manually")
    }
    
    # --- Ensure time vector is strictly positive ------------------------------
    time_was_shifted <- FALSE
    shift_amount <- 0
    if (any(time_vec <= 0)) {
      if (positive_offset_policy == "shift") {
        step <- if (length(time_vec) > 1) median(diff(sort(unique(time_vec)))) else 1
        shift_amount <- abs(min(time_vec)) + step
        time_vec <- time_vec + shift_amount
        colnames(mat) <- as.character(time_vec)
        time_was_shifted <- TRUE
        notes <- c(notes, sprintf("time vector had non-positive values; shifted by +%.4g", shift_amount))
      } else {
        notes <- c(notes, "time vector has non-positive values (policy = 'error')")
      }
    }
    
    is_valid <- (nrow(mat) == n_target) &&
      length(missing_nodes) == 0 &&
      all(time_vec > 0) &&
      !anyNA(mat)
    
    diag <- data.frame(
      file              = file_name,
      nrow              = nrow(mat),
      ncol              = ncol(mat),
      n_nodes_expected  = n_target,
      n_nodes_missing   = length(missing_nodes),
      missing_nodes     = paste(missing_nodes, collapse = "; "),
      was_transposed    = was_transposed,
      had_time_row      = had_time_row,
      time_was_shifted  = time_was_shifted,
      shift_amount      = shift_amount,
      min_time          = min(time_vec),
      max_time          = max(time_vec),
      mean_value        = mean(mat, na.rm = TRUE),
      sd_value          = sd(as.numeric(mat), na.rm = TRUE),
      has_na            = anyNA(mat),
      is_valid          = is_valid,
      notes             = paste(notes, collapse = " | "),
      error             = NA_character_,
      stringsAsFactors  = FALSE
    )
    
    list(data = mat, time = time_vec, diagnostics = diag)
    
  }, error = function(e) {
    diag <- data.frame(
      file = file_name, nrow = NA, ncol = NA, n_nodes_expected = NA,
      n_nodes_missing = NA, missing_nodes = NA, was_transposed = NA,
      had_time_row = NA, time_was_shifted = NA, shift_amount = NA,
      min_time = NA, max_time = NA, mean_value = NA, sd_value = NA,
      has_na = NA, is_valid = FALSE, notes = NA,
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
    list(data = NULL, time = NULL, diagnostics = diag)
  })
  
  out
}

# ---- 3. Wrapper: loop over every CSV in a folder -------------------------

process_eeg_folder <- function(dir,
                               expected_nodes = NULL,
                               pattern = "\\.csv$",
                               sampling_rate = NULL,
                               out_dir = NULL) {
  
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No CSV files found in: ", dir)
  
  results <- vector("list", length(files))
  names(results) <- basename(files)
  diag_rows <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    res <- check_and_fix_eeg(files[i], expected_nodes = expected_nodes,
                             sampling_rate = sampling_rate)
    results[[i]] <- res
    diag_rows[[i]] <- res$diagnostics
    
    if (!is.null(out_dir) && res$diagnostics$is_valid) {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      write.csv(res$data, file.path(out_dir, basename(files[i])), row.names = TRUE)
    }
  }
  
  summary_table <- do.call(rbind, diag_rows)
  rownames(summary_table) <- NULL
  
  list(
    cleaned   = results,        # named list: each$data, each$time, each$diagnostics
    summary   = summary_table   # one row per file, like your summarise() idea
  )
}