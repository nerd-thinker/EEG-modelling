# libraries -----------
install.packages(c('pak', 'mgcv', 'pracma', 'readxl'))
library(mgcv)
library(readxl)
library(pracma)
library(tidyverse)

# ── Helper functions ──────────────────────────────────────────────────────────

hjorth <- function(x) {
  x   <- na.omit(x)
  dx  <- diff(x)
  ddx <- diff(dx)
  activity   <- var(x)
  mobility   <- sqrt(var(dx) / var(x))
  complexity <- (sqrt(var(ddx) / var(dx))) / mobility
  list(activity = activity, mobility = mobility, complexity = complexity)
}

skewness     <- function(x) mean(((x - mean(x)) / sd(x))^3)
kurtosis_fn  <- function(x) mean(((x - mean(x)) / sd(x))^4)
band_power   <- function(x) sum(x^2) / length(x)

spectral_entropy <- function(x) {
  spec <- abs(fft(x))^2
  spec <- spec / sum(spec)
  spec <- spec[spec > 0]
  -sum(spec * log(spec))
}

# ── Single-channel feature extractor ─────────────────────────────────────────
# Returns a named numeric vector of features for ONE (signal, scaling) combo.
# 'scaling' is one of: "none" | "z_score" | "01" | "log10"
# log10 transform: applied to raw x BEFORE further processing.
#   Because EEG band-power values can be ≤ 0, we shift so min(x) > 0 first:
#   log10(x - min(x) + 1).  This keeps the transform well-defined and
#   preserves the shape while compressing large spikes.

get_features_single <- function(x,
                                time    = NULL,
                                scaling = c("none", "z_score", "01", "log10"),
                                k       = 20,
                                plot    = FALSE) {
  
  if (is.null(time)) time <- seq_along(x)
  scaling <- match.arg(scaling)
  
  # ── 1. Scale / transform ───────────────────────────────────────────────────
  norm_x <- switch(
    scaling,
    
    "none"   = x,
    
    "z_score" = as.numeric(scale(x)),
    
    "01" = if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
      rep(0, length(x))
    } else {
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    },
    
    # log10 of shifted signal (guarantees all values > 0 before log)
    "log10" = {
      x_shifted <- x - min(x, na.rm = TRUE) + 1   # shift so min = 1
      log10(x_shifted)
    }
  )
  
  # ── 2. GAM smooth ─────────────────────────────────────────────────────────
  fit      <- gam(norm_x ~ s(time, k = k))
  smooth_x <- as.numeric(predict(fit, data.frame(time = time)))
  
  # ── 3. Peak detection on smooth signal ────────────────────────────────────
  peaks <- findpeaks(smooth_x,
                     minpeakheight   = mean(smooth_x) + sd(smooth_x),
                     minpeakdistance = 10)
  
  if (!is.null(peaks)) {
    peak_amplitudes <- peaks[, 1]
    peak_indices    <- peaks[, 2]
    peak_times      <- time[peak_indices]
    
    num_peaks <- nrow(peaks)
    mean_amp  <- mean(peak_amplitudes)
    max_amp   <- max(peak_amplitudes)
    sd_amp    <- ifelse(num_peaks > 1, sd(peak_amplitudes), NA_real_)
    
    ipi       <- diff(peak_times)
    mean_ipi  <- ifelse(length(ipi) > 0, mean(ipi), NA_real_)
    sd_ipi    <- ifelse(length(ipi) > 1, sd(ipi),   NA_real_)
    peak_freq <- num_peaks / (max(time) - min(time))
  } else {
    num_peaks <- 0L
    mean_amp  <- NA_real_; max_amp  <- NA_real_; sd_amp   <- NA_real_
    mean_ipi  <- NA_real_; sd_ipi   <- NA_real_; peak_freq <- 0
    peak_indices <- integer(0)
  }
  
  # ── 4. Optional plot ──────────────────────────────────────────────────────
  if (plot) {
    plot(time, norm_x, type = "l", xlab = "Time", ylab = "Signal", col = "grey87")
    lines(time, smooth_x, col = 2, lwd = 2)
    if (length(peak_indices) > 0)
      abline(v = time[peak_indices], lty = 2, col = "grey20")
  }
  
  # ── 5. Extra features ─────────────────────────────────────────────────────
  smooth_ds <- smooth_x[seq(1, length(smooth_x), by = 8)]
  se <- tryCatch(
    pracma::sample_entropy(smooth_ds, edim = 2, r = 0.2 * sd(smooth_ds)),
    error = function(e) NA_real_
  )
  
  spec     <- abs(fft(smooth_x))^2 / length(smooth_x)
  spec     <- spec[1:(length(spec) / 2)]
  peak_psd <- which.max(spec)
  mean_psd <- mean(spec)
  
  h <- hjorth(smooth_x)
  
  # ── 6. Return named vector ─────────────────────────────────────────────────
  c(
    mean                     = mean(smooth_x),
    skewness                 = skewness(smooth_x),
    kurtosis                 = kurtosis_fn(smooth_x),
    q25                      = unname(quantile(smooth_x, .25)),
    q75                      = unname(quantile(smooth_x, .75)),
    n_peaks                  = num_peaks,
    peak_frequency           = peak_freq,
    mean_peak_amplitude      = mean_amp,
    max_peak_amplitude       = max_amp,
    sd_peak_amplitude        = sd_amp,
    mean_inter_peak_interval = mean_ipi,
    sd_inter_peak_interval   = sd_ipi,
    activity                 = h$activity,
    mobility                 = h$mobility,
    complexity               = h$complexity,
    sample_entropy           = se,
    spectral_entropy         = spectral_entropy(smooth_x),
    abs_band_power           = band_power(smooth_x),
    peak_psd_freq            = peak_psd,
    mean_psd                 = mean_psd
  )
}


# ── Master function ───────────────────────────────────────────────────────────
# Returns a single-row data frame:
#   columns = <band>_<channel>_<scaling>_<feature>   (very wide)
#   rows    = one per participant (participant_id passed in)
#
# All four scalings (none, z_score, 01, log10) are always computed.

extract_all_features <- function(filepaths,
                                 participant_id = "P001",
                                 k              = 20) {
  
  scalings   <- c("none", "z_score", "01", "log10")
  start_time <- proc.time()
  
  # ── 1. Read files ──────────────────────────────────────────────────────────
  message("Reading files...")
  bands     <- names(filepaths)
  band_data <- lapply(setNames(filepaths, bands), function(path) {
    df <- if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) read_excel(path)
    else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
      read_csv(path, show_col_types = FALSE)
    }
    else stop("Unsupported file type: ", path)
    
    # Drop phantom columns created by a trailing comma in the header row
    df[ , !grepl("^\\.\\.\\.\\d+$", colnames(df)), drop = FALSE]
  })
  
  time_vec <- band_data[[1]]$Time
  channels <- setdiff(colnames(band_data[[1]]), "Time")
  
  message(sprintf("  %d channels, %d bands, %d scalings → %d combinations",
                  length(channels), length(bands), length(scalings),
                  length(channels) * length(bands) * length(scalings)))
  
  # ── 2. Per-band × per-channel × per-scaling features ──────────────────────
  message("Extracting per-band features...")
  
  feature_list <- list()   # accumulate named vectors
  
  for (b in bands) {
    for (ch in channels) {
      x <- as.numeric(band_data[[b]][[ch]])
      
      for (sc in scalings) {
        col_prefix <- paste(b, ch, sc, sep = "_")   # e.g. "alpha_Fp1_z_score"
        
        feats <- tryCatch(
          get_features_single(x, time = time_vec, scaling = sc, k = k),
          error = function(e) {
            message(sprintf("  FAILED: %s — %s", col_prefix, e$message))
            rep(NA_real_, 20)
          }
        )
        
        # Prepend the prefix to every feature name
        names(feats) <- paste(col_prefix, names(feats), sep = "_")
        feature_list <- c(feature_list, as.list(feats))
      }
    }
  }
  
  # ── 3. Cross-band features (computed on raw signal, once per channel) ──────
  message("Computing cross-band features...")
  
  abs_powers <- sapply(setNames(bands, bands), function(b) {
    sapply(channels, function(ch) band_power(as.numeric(band_data[[b]][[ch]])))
  })  # matrix: channels × bands
  
  total_power <- rowSums(abs_powers)
  
  # Relative powers per channel-band
  rel_power_feats <- as.list(
    setNames(
      as.vector(abs_powers / total_power),
      paste0("cross_", rep(channels, length(bands)),
             "_rel_power_", rep(bands, each = length(channels)))
    )
  )
  
  # Band-ratio features per channel
  ratio_feats <- lapply(setNames(channels, channels), function(ch) {
    i <- which(channels == ch)
    alpha_p <- abs_powers[i, "alpha"]
    beta_p  <- abs_powers[i, "beta"]
    theta_p <- abs_powers[i, "theta"]
    delta_p <- abs_powers[i, "delta"]
    setNames(
      c(theta_p / alpha_p,
        alpha_p / beta_p,
        delta_p / beta_p,
        theta_p / beta_p,
        (alpha_p + theta_p) / beta_p),
      paste0("cross_", ch, "_",
             c("ratio_theta_alpha", "ratio_alpha_beta",
               "ratio_delta_beta",  "ratio_theta_beta",
               "ratio_stress"))
    )
  }) %>% unlist()
  
  cross_list <- c(rel_power_feats, as.list(ratio_feats))
  
  # ── 4. Assemble single wide row ────────────────────────────────────────────
  all_feats <- c(list(participant_id = participant_id),
                 feature_list,
                 cross_list)
  
  result_df <- as.data.frame(all_feats, check.names = FALSE)
  
  elapsed <- (proc.time() - start_time)["elapsed"]
  message(sprintf("\nDone! %.1f seconds (%.1f min) | %d feature columns",
                  elapsed, elapsed / 60, ncol(result_df) - 1))
  
  result_df
}


# ── Run ───────────────────────────────────────────────────────────────────────

filepaths <- c(
  alpha = "~/eeg/raw_data/CTEEG022Y_Alpha.xlsx",
  beta  = "~/eeg/raw_data/CTEEG022Y_Beta.xlsx",
  gamma = "~/eeg/raw_data/CTEEG022Y_Gamma1.xlsx",
  delta = "~/eeg/raw_data/CTEEG022Y_Delta.xlsx",
  theta = "~/eeg/raw_data/CTEEG022Y_Theta.xlsx"
)
# "raw_data/SPUR EEG data/001Y_Alpha.csv"
filepaths <- c(
  alpha = "~/eeg/raw_data/SPUR EEG data/001Y_Alpha.csv",
  beta  = "~/eeg/raw_data/SPUR EEG data/001Y_Beta.csv",
  gamma = "~/eeg/raw_data/SPUR EEG data/001Y_Delta.csv",
  delta = "~/eeg/raw_data/SPUR EEG data/001Y_Gamma1.csv",
  theta = "~/eeg/raw_data/SPUR EEG data/001Y_Theta.csv"
)

participant_row <- extract_all_features(
  filepaths      = filepaths,
  participant_id = "001Y",
  k              = 70
)

# Result: 1 row, very many columns
dim(participant_row)          # 1 × (1 + n_features)
names(participant_row)[1:10]  # preview column names

# Save to CSV for later stacking with other participants
write.csv(participant_row,
          file      = "CTEEG022Y_features.csv",
          row.names = FALSE)


# ── Stacking participants later ───────────────────────────────────────────────
# When you have multiple participants, just run extract_all_features() for each
# and bind the rows:
#
# all_participants <- bind_rows(
#   extract_all_features(filepaths_P001, participant_id = "P001"),
#   extract_all_features(filepaths_P002, participant_id = "P002"),
#   ...
# )
#
# Or in a loop / lapply over a list of filepath-sets:
#
# participant_file_list <- list(
#   CTEEG022Y = c(alpha = "...", beta = "...", ...),
#   CTEEG023Y = c(alpha = "...", beta = "...", ...)
# )
#
# all_participants <- lapply(names(participant_file_list), function(pid) {
#   extract_all_features(participant_file_list[[pid]], participant_id = pid)
# }) %>% bind_rows()

## NA test --------------
lapply(participant_row, function(x) {
  y <- apply(x, function(x) which(any(is.na(x))))
  return(y)
}) %>% unlist

# NA troubleshooting ------------------
which(any(is.na(participant_row)), arr.ind = T)
any(is.na(participant_row))

# Find which columns are NA
na_cols <- names(participant_row)[is.na(participant_row)]
data.frame(na_columns = na_cols) ## all of the inter_peak_interals

# How many NAs total?
sum(is.na(participant_row))

# Check how many peaks these channels have
get_features_single(
  x       = as.numeric(band_data$beta$F4),
  time    = time_vec,
  scaling = "z_score",
  k       = 20,
  plot    = TRUE   # visualise
)

