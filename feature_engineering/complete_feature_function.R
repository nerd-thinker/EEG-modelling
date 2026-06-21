# libraries ----
library(readxl)
library(mgcv)
library(pracma)
library(tidyverse)

# ── Helper functions ──────────────────────────────────────────────────────────

hjorth <- function(x) {
  x   <- na.omit(x)
  dx  <- diff(x)
  ddx <- diff(dx)
  list(
    activity   = var(x),
    mobility   = sqrt(var(dx) / var(x)),
    complexity = (sqrt(var(ddx) / var(dx))) / sqrt(var(dx) / var(x))
  )
}

skewness         <- function(x) mean(((x - mean(x)) / sd(x))^3)
kurtosis_fn      <- function(x) mean(((x - mean(x)) / sd(x))^4)
band_power       <- function(x) sum(x^2) / length(x)
spectral_entropy <- function(x) {
  spec <- abs(fft(x))^2
  spec <- spec / sum(spec)
  spec <- spec[spec > 0]
  -sum(spec * log(spec))
}

# ── Single-channel feature extraction (one scaling at a time) ─────────────────
# Called inside extract_all_features() for every band × channel × scaling combo.
# scaling: "none" | "z_score" | "01" | "log10"

get_features_single <- function(x,
                                time    = NULL,
                                scaling = c("none", "z_score", "01", "log10"),
                                k       = 70,
                                plot    = FALSE) {
  if (is.null(time)) time <- seq_along(x)
  scaling <- match.arg(scaling)
  
  # Scale / transform
  norm_x <- switch(
    scaling,
    "none"    = x,
    "z_score" = as.numeric(scale(x)),
    "01"      = if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
      rep(0, length(x))
    } else {
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    },
    "log10"   = log10(x - min(x, na.rm = TRUE) + 1)
  )
  
  # GAM smooth
  fit      <- gam(norm_x ~ s(time, k = k))
  smooth_x <- as.numeric(predict(fit, data.frame(time = time)))
  
  # Peak detection on smoothed signal
  peaks <- findpeaks(smooth_x,
                     minpeakheight   = mean(smooth_x) + sd(smooth_x),
                     minpeakdistance = 10)
  
  if (!is.null(peaks)) {
    peak_amplitudes <- peaks[, 1]
    peak_indices    <- peaks[, 2]
    peak_times      <- time[peak_indices]
    num_peaks       <- nrow(peaks)
    mean_amp        <- mean(peak_amplitudes)
    max_amp         <- max(peak_amplitudes)
    sd_amp          <- ifelse(num_peaks > 1, sd(peak_amplitudes), NA_real_)
    ipi             <- diff(peak_times)
    mean_ipi        <- ifelse(length(ipi) > 0, mean(ipi), NA_real_)
    sd_ipi          <- ifelse(length(ipi) > 1, sd(ipi),   NA_real_)
    peak_freq       <- num_peaks / (max(time) - min(time))
  } else {
    num_peaks <- 0L
    mean_amp  <- NA_real_; max_amp  <- NA_real_; sd_amp   <- NA_real_
    mean_ipi  <- NA_real_; sd_ipi   <- NA_real_; peak_freq <- 0
    peak_indices <- integer(0)
  }
  
  if (plot) {
    plot(time, norm_x, type = "l", xlab = "Time", ylab = "Signal", col = "grey87")
    lines(time, smooth_x, col = 2, lwd = 2)
    if (length(peak_indices) > 0) abline(v = time[peak_indices], lty = 2, col = "grey20")
  }
  
  # Sample entropy on downsampled smooth signal
  smooth_ds <- smooth_x[seq(1, length(smooth_x), by = 8)]
  se <- tryCatch(
    pracma::sample_entropy(smooth_ds, edim = 2, r = 0.2 * sd(smooth_ds)),
    error = function(e) NA_real_
  )
  
  # Power spectral density
  spec     <- abs(fft(smooth_x))^2 / length(smooth_x)
  spec     <- spec[1:(length(spec) / 2)]
  peak_psd <- which.max(spec)
  mean_psd <- mean(spec)
  
  h <- hjorth(smooth_x)
  
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
# Returns one wide row per participant.
# Columns: <band>_<channel>_<scaling>_<feature>  +  cross_<channel>_<ratio>
# All four scalings (none, z_score, 01, log10) are always computed.

extract_all_features <- function(filepaths,
                                 participant_id = "unknown",
                                 k              = 70) {
  
  scalings    <- c("none", "z_score", "01", "log10")
  start_total <- proc.time()
  bands       <- names(filepaths)
  
  # ── 1. Read files ────────────────────────────────────────────────────────
  message("Reading files...")
  band_data <- lapply(setNames(filepaths, bands), function(path) {
    df <- if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
      read_excel(path)
    } else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
      suppressWarnings(read_csv(path, show_col_types = FALSE))
    } else {
      stop("Unsupported file type: ", path)
    }
    # Drop phantom columns from trailing comma in CSV header
    df[ , !grepl("^\\.\\.\\.\\d+$", colnames(df)), drop = FALSE]
  })
  
  time_vec <- band_data[[1]]$Time
  channels <- setdiff(colnames(band_data[[1]]), "Time")
  
  message(sprintf("  %d channels, %d bands, %d scalings -> %d combinations",
                  length(channels), length(bands), length(scalings),
                  length(channels) * length(bands) * length(scalings)))
  
  # ── 2. Per-band x per-channel x per-scaling features ─────────────────────
  message("Extracting per-band features...")
  
  feature_list <- list()
  
  for (b in bands) {
    for (ch in channels) {
      x          <- as.numeric(band_data[[b]][[ch]])
      col_prefix <- paste(b, ch, sep = "_")
      
      for (sc in scalings) {
        feat_prefix <- paste(col_prefix, sc, sep = "_")
        
        feats <- tryCatch(
          get_features_single(x, time = time_vec, scaling = sc, k = k),
          error = function(e) {
            message(sprintf("  FAILED: %s - %s", feat_prefix, e$message))
            rep(NA_real_, 20)
          }
        )
        
        names(feats) <- paste(feat_prefix, names(feats), sep = "_")
        feature_list <- c(feature_list, as.list(feats))
      }
    }
  }
  
  # ── 3. Cross-band features ────────────────────────────────────────────────
  message("Computing cross-band features...")
  
  abs_powers <- sapply(setNames(bands, bands), function(b) {
    sapply(channels, function(ch) band_power(as.numeric(band_data[[b]][[ch]])))
  })
  
  total_power <- rowSums(abs_powers)
  
  rel_power_feats <- as.list(
    setNames(
      as.vector(abs_powers / total_power),
      paste0("cross_", rep(channels, length(bands)),
             "_rel_power_", rep(bands, each = length(channels)))
    )
  )
  
  # Safe accessor — returns NA vector if band is missing for this participant
  get_power <- function(band) {
    if (band %in% colnames(abs_powers)) abs_powers[, band]
    else rep(NA_real_, length(channels))
  }
  
  ratio_feats <- data.frame(
    ratio_theta_alpha = get_power("theta") / get_power("alpha"),
    ratio_alpha_beta  = get_power("alpha") / get_power("beta"),
    ratio_delta_beta  = get_power("delta") / get_power("beta"),
    ratio_theta_beta  = get_power("theta") / get_power("beta"),
    ratio_stress      = (get_power("alpha") + get_power("theta")) / get_power("beta")
  )
  
  ratio_list <- as.list(
    setNames(
      unlist(ratio_feats),
      paste0("cross_", rep(channels, ncol(ratio_feats)),
             "_", rep(names(ratio_feats), each = length(channels)))
    )
  )
  
  cross_list <- c(rel_power_feats, ratio_list)
  
  # ── 4. Assemble single wide row ───────────────────────────────────────────
  all_feats <- c(list(participant_id = participant_id),
                 feature_list,
                 cross_list)
  
  result_df <- as.data.frame(all_feats, check.names = FALSE)
  
  elapsed <- (proc.time() - start_total)["elapsed"]
  message(sprintf("\nDone! %.1f seconds (%.1f min) | %d feature columns",
                  elapsed, elapsed / 60, ncol(result_df) - 1))
  
  result_df
}