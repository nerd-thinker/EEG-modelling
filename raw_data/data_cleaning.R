# =============================================================================
# normalise_eeg.R
# -----------------------------------------------------------------------------
# Brings two known EEG data structures into one standard form:
#
#   TYPE 1 (long):  rows = time points  |  col 1 = Time, remaining = nodes
#                   (may have extra columns like "X" that are dropped)
#   TYPE 2 (wide):  rows = 32 nodes     |  cols  = time points (no Time col)
#                   (may have no header row at all — first row is data)
#
#   TARGET:         rows = time points  |  col 1 = Time (fresh, positive)
#                                       |  cols 2-33 = nodes (named)
#
# Since this is resting EEG, the time vector carries no event-locked meaning
# and is regenerated identically for every file from the sampling rate and
# sample count. This eliminates any cross-participant timing inconsistency.
# =============================================================================

# ---- CONFIG -----------------------------------------------------------------

NODE_NAMES <- c(
  "Fp1","AF3","F7","F3","FC1","FC5","T7","C3","CP1","CP5",
  "P7","P3","Pz","PO3","O1","Oz","O2","PO4","P4","P8",
  "CP6","CP2","C4","T8","FC6","FC2","F4","F8","AF4","Fp2",
  "Fz","Cz"
)

SAMPLING_RATE <- 128   # Hz

# ---- SINGLE FILE ------------------------------------------------------------

normalise_eeg <- function(path) {
  
  file <- basename(path)
  
  # ---- Smart read: detect whether file has a real header row ---------------
  # If every field in row 1 parses as a number, there is no header and
  # read.csv(header=TRUE) would silently eat the first data row.
  first_line <- strsplit(readLines(path, n = 1), ",")[[1]]
  has_header <- any(is.na(suppressWarnings(as.numeric(trimws(first_line)))))
  
  raw <- read.csv(path, header = has_header, check.names = FALSE,
                  stringsAsFactors = FALSE)
  nr  <- nrow(raw)
  nc  <- ncol(raw)
  
  # ---- Detect type by Time column, not exact column count -----------------
  first_col_name <- tolower(trimws(names(raw)[1]))
  is_time_col    <- first_col_name %in% c("time","t","timestamp","sec","secs","ms")
  
  if (is_time_col) {
    
    # ---- TYPE 1 -------------------------------------------------------------
    # Keep only columns whose name matches NODE_NAMES (drops "X" and any other
    # stray columns). Match case-insensitively for safety.
    col_names_upper <- toupper(trimws(names(raw)))
    node_idx        <- which(col_names_upper %in% toupper(NODE_NAMES))
    
    dropped <- names(raw)[!seq_along(names(raw)) %in% c(1, node_idx)]
    ## delete if not needed
    if (length(dropped) > 0) 
      message(sprintf("[%s] Dropped %d unrecognised column(s): %s",
                      file, length(dropped), paste(dropped, collapse = ", ")))
    
    if (length(node_idx) != 32) {
      warning(sprintf("[%s] Type 1 detected but found %d recognised node columns (expected 32). Skipping.",
                      file, length(node_idx)))
      return(NULL)
    }
    
    node_data        <- raw[, node_idx, drop = FALSE]
    # Re-apply canonical casing in the fixed NODE_NAMES order
    canonical_order  <- match(toupper(NODE_NAMES), toupper(names(node_data)))
    node_data        <- node_data[, canonical_order, drop = FALSE]
    names(node_data) <- NODE_NAMES
    
  } else if (nr == 32) {
    
    # ---- TYPE 2 -------------------------------------------------------------
    node_data           <- as.data.frame(t(raw))
    rownames(node_data) <- NULL
    names(node_data)    <- NODE_NAMES
    
  } else {
    warning(sprintf(
      "[%s] Could not identify as Type 1 (Time header column) or Type 2 (32 rows). Skipping. nrow=%d ncol=%d has_header=%s",
      file, nr, nc, has_header))
    return(NULL)
  }
  
  # ---- Build a fresh time vector ------------------------------------------
  n_samples   <- nrow(node_data)
  time_vector <- seq(1/SAMPLING_RATE, n_samples/SAMPLING_RATE, by = 1/SAMPLING_RATE)
  
  # ---- Assemble target data.frame -----------------------------------------
  result           <- cbind(Time = time_vector, node_data)
  rownames(result) <- NULL
  result
}

# ---- FOLDER LOOP ------------------------------------------------------------

normalise_eeg_folder <- function(dir, out_dir = NULL, pattern = "\\.csv$") {
  
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No CSV files found in: ", dir)
  
  if (!is.null(out_dir) && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  cleaned      <- vector("list", length(files))
  names(cleaned) <- basename(files)
  summary_rows <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    
    df <- tryCatch(normalise_eeg(files[i]), error = function(e) {
      warning(sprintf("[%s] Error: %s", basename(files[i]), e$message))
      NULL
    })
    
    cleaned[[i]] <- df
    
    if (!is.null(df)) {
      ok <- ncol(df) == 33 && names(df)[1] == "Time" &&
        all(NODE_NAMES %in% names(df)) && all(df$Time > 0)
      
      fl     <- strsplit(readLines(files[i], n = 1), ",")[[1]]
      t1_hdr <- any(is.na(suppressWarnings(as.numeric(trimws(fl))))) &&
        tolower(trimws(fl[1])) %in% c("time","t","timestamp","sec","secs","ms")
      
      summary_rows[[i]] <- data.frame(
        file          = basename(files[i]),
        type_detected = if (t1_hdr) "Type1" else "Type2",
        n_samples     = nrow(df),
        duration_s    = max(df$Time),
        n_nodes       = ncol(df) - 1,
        is_valid      = ok,
        error         = NA_character_,
        stringsAsFactors = FALSE
      )
      if (!is.null(out_dir))
        write.csv(df, file.path(out_dir, basename(files[i])), row.names = FALSE)
      
    } else {
      summary_rows[[i]] <- data.frame(
        file = basename(files[i]), type_detected = "unknown",
        n_samples = NA, duration_s = NA, n_nodes = NA,
        is_valid = FALSE, error = "see warning",
        stringsAsFactors = FALSE
      )
    }
  }
  
  list(
    data    = cleaned,
    summary = do.call(rbind, summary_rows)
  )
}
# Execution -------------
normalised_data <- normalise_eeg_folder("~/eeg/raw_data/SPUR-EEG-data", out_dir = "~/eeg/raw_data")

