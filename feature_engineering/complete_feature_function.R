# libraries -----------
library(readxl)
library(mgcv)
library(pracma)
library(tidyverse)
# library(dplyr)

# Helper functions -------------
hjorth <- function(x) {
  x <- na.omit(x)
  
  dx  <- diff(x)        # first derivative
  ddx <- diff(dx)       # second derivative
  
  activity   <- var(x)
  mobility   <- sqrt(var(dx) / var(x)) ## estimate of the mean frequency
  complexity <- (sqrt(var(ddx) / var(dx))) / mobility ## estimate of the bandwidth of the signal, which indicates the similarity of the shape of the signal to a pure sine wave
  
  list(activity = activity, mobility = mobility, complexity = complexity)
}

skewness <- function(x) mean(((x - mean(x)) / sd(x))^3)
kurtosis <- function(x) mean(((x - mean(x)) / sd(x))^4)
band_power <- function(x) sum(x^2) / length(x)

spectral_entropy <- function(x) {
  spec <- abs(fft(x))^2
  spec <- spec / sum(spec)
  spec <- spec[spec > 0]
  -sum(spec * log(spec))
}

# Single channel feature extraction -----------

get_features <- function(x,
                         time = NULL,
                         scaling = c("none","z_score","01"),
                         k = 20,
                         plot = FALSE) {
  if(is.null(time)){ 
    time <- seq_along(x)
  } ## time needs to be provided as real time in seconds node$Time != sample indices (1:length(x))
  scaling <- match.arg(scaling)
  norm_x <- switch(scaling,
                   "none" = x,
                   "z_score" = as.numeric(scale(x)),
                   "01" = if (max(x) == min(x)) {   ## catching division by 0
                     rep(0, length(x))
                   } else {
                     (x - min(x)) / (max(x) - min(x))
                   })
  
  # GAM smooth 
  fit <- gam(norm_x ~ s(time, k = k))
  smooth_x    <- as.numeric(predict(fit, data.frame(time = time)))  
  
  # Detect peaks on smooth signal
  peaks <- findpeaks(smooth_x,
                     minpeakheight   = mean(smooth_x) + sd(smooth_x),
                     minpeakdistance = 10)
  
  if (!is.null(peaks)) {
    peak_amplitudes <- peaks[, 1]          # peak heights
    peak_indices    <- peaks[, 2]          # positions smooth_x
    peak_times      <- time[peak_indices]  # map back to original time scale
    
    num_peaks    <- nrow(peaks)            
    mean_amp     <- mean(peak_amplitudes)
    max_amp      <- max(peak_amplitudes)
    sd_amp       <- ifelse(num_peaks > 1, sd(peak_amplitudes), NA)
    
    ipi          <- diff(peak_times)
    mean_ipi     <- ifelse(length(ipi) > 0, mean(ipi), NA)
    sd_ipi       <- ifelse(length(ipi) > 1, sd(ipi),   NA)
    peak_freq    <- num_peaks / (max(time) - min(time))
    
  } else {
    message("No peaks detected.")
    num_peaks <- 0
    mean_amp  <- NA;  max_amp   <- NA;  sd_amp    <- NA
    mean_ipi  <- NA;  sd_ipi    <- NA;  peak_freq <- 0
    peak_amplitudes <- numeric(0)
    peak_indices    <- integer(0)
  }
  
  # Plot
  if(plot == TRUE){
    plot(time, norm_x, type = "l",
         xlab = "Time", ylab = "Signal",
         col = "grey87")                              # raw signal in grey
    lines(time, smooth_x, col = 2, lwd = 2)   # GAM smooth in red
    
    if (length(peak_indices) > 0)
      abline(v = time[peak_indices], lty = 2, col = "grey20")
  }
  # Sample entropy on downsapled signal
  smooth_ds <- smooth_x[seq(1,length(smooth_x), by = 8)]
  se        <- tryCatch(
    pracma::sample_entropy(smooth_ds, edim = 2, r = .2 * sd(smooth_ds)),
    error = function(e) NA
  )
  # Power Spectral Density (PSD) features
  spec      <- abs(fft(smooth_x))^2 / length(smooth_x)
  spec      <- spec[1:(length(spec) / 2)]
  peak_psd  <- which.max(spec)
  mean_psd  <- mean(spec)
  
  h <- hjorth(smooth_x)
  
  c(
    # General smooth features
    "mean"                    = mean(smooth_x),
    "skewness"               = skewness(smooth_x),
    "kurtosis"                = kurtosis(smooth_x),
    "q25"                     = unname(quantile(smooth_x, .25)),
    "q75"                     = unname(quantile(smooth_x, .75)),
    
    # Peak features
    "n_peaks"                 = num_peaks,
    "peak_frequency"          = peak_freq, ##have if time for different participants is nto equal
    "mean_peak_amplitude"     = mean_amp,
    "max_peak_amplitude"      = max_amp,
    "sd_peak_amplitude"       = sd_amp,
    "mean_inter_peak_interval"= mean_ipi,
    "sd_inter_peak_interval"  = sd_ipi,
    
    # Hjorth
    "activity"                = h$activity, ##same as variance
    "mobility"                = h$mobility,
    "complexity"              = h$complexity,
    
    # entropy
    "sample_entropy"          = se,
    "spectral_entropy"        = spectral_entropy(smooth_x),
    
    # band power
    "abs_band_power"          = band_power(smooth_x),
    "peak_psd_freq"           = peak_psd,
    "mean_psd"                = mean_psd
  )
}





# MASTER FUNCTION -------------------
## Data assumptions: 
## data is wide type, data has columns: Time, Fp1, other nodes, Time starts at 0
extract_all_features <- function(filepaths,
                                 scaling = "z_score",
                                 k       = 20) {
  start_total <- proc.time()
  bands <- names(filepaths)  # band names come from the named vector
  
  # Step 1: Read all files
  message("Reading Excel files...")
  band_data <- lapply(setNames(filepaths, bands), function(fp) {
    read_excel(fp)
  })
  
  # Extract time and channels from first file
  time_vec <- band_data[[1]]$Time
  channels <- setdiff(colnames(band_data[[1]]), "Time")
  
  message(sprintf("Found %d channels across %d bands", length(channels), length(bands)))
  
  # Step 2: Per-channel per-band features
  message("Extracting features per channel per band...")
  
  per_band <- lapply(setNames(bands, bands), function(b) {
    
    band_df <- band_data[[b]]
    
    result <- sapply(channels, function(ch) {
      x <- as.numeric(band_df[[ch]])
      tryCatch(
        get_features(x,
                     time    = time_vec,
                     scaling = scaling,
                     k       = k,
                     plot    = FALSE),
        error = function(e) {
          message(sprintf("  Failed: band = %s, channel = %s — %s", b, ch, e$message))
          rep(NA, 20)
        }
      )
    })
    
    as.data.frame(t(result)) %>%
      mutate(channel = channels, band = b) %>%
      relocate(band, channel)
  })
  
  # Step 3: Cross-band features
  message("Computing cross-band features...")
  
  abs_powers <- sapply(setNames(bands, bands), function(b) {
    sapply(channels, function(ch) {
      band_power(as.numeric(band_data[[b]][[ch]]))
    })
  })
  
  total_power <- rowSums(abs_powers)
  rel_powers  <- abs_powers / total_power
  colnames(rel_powers) <- paste0("rel_power_", bands)
  
  power_ratios <- data.frame(
    ratio_theta_alpha = abs_powers[, "theta"] / abs_powers[, "alpha"],
    ratio_alpha_beta  = abs_powers[, "alpha"] / abs_powers[, "beta"],
    ratio_delta_beta  = abs_powers[, "delta"] / abs_powers[, "beta"],
    ratio_theta_beta  = abs_powers[, "theta"] / abs_powers[, "beta"],
    ratio_stress      = (abs_powers[, "alpha"] + abs_powers[, "theta"]) / abs_powers[, "beta"]
  )
  
  peak_psd_by_band <- sapply(setNames(bands, bands), function(b) {
    sapply(channels, function(ch) {
      x    <- as.numeric(band_data[[b]][[ch]])
      # Gam smoothing
      time <- seq_along(x)
      fit  <- gam(x ~ s(time, k=k))
      sx   <- as.numeric(predict(fit))
      
      # PSD - skip bin 1(DC component)
      spec <- abs(fft(sx))^2 / length(sx)
      spec <- spec[2:(length(spec) / 2)] ##start from bin 2  
      which.max(spec) + 1                ##index correction
    })
  })
  colnames(peak_psd_by_band) <- paste0("peak_psd_", bands)
  
  cross_band <- as.data.frame(cbind(rel_powers, power_ratios, peak_psd_by_band)) %>%
    mutate(channel = channels) %>%
    relocate(channel)
  
  total_time <- proc.time() -start_total
  message(sprintf("\nDone! Total time: %.1f seconds (%.1f minutes)", total_time["elapsed"], total_time["elapsed"]/60))
  
  list(per_band   = per_band, cross_band = cross_band)
}

# Run it
# Step 1: name your filepaths by band
filepaths <- c(
  alpha = "/home/ks/EGG-modeling/raw_data/CTEEG022Y_Alpha.xlsx",
  beta  = "/home/ks/EGG-modeling/raw_data/CTEEG022Y_Beta.xlsx",
  gamma = "/home/ks/EGG-modeling/raw_data/CTEEG022Y_Gamma1.xlsx", 
  delta = "/home/ks/EGG-modeling/raw_data/CTEEG022Y_Delta.xlsx",
  theta = "/home/ks/EGG-modeling/raw_data/CTEEG022Y_Theta.xlsx"
)
results <- extract_all_features(
  filepaths = filepaths,
  scaling   = "z_score",
  k         = 20
)

# Results -----------------
results$per_band$alpha    # alpha band — channels x features data frame
results$cross_band         # cross-band features — channels x features data frame

# Combine all bands into one long data frame
all_features <- bind_rows(results$per_band)

# Possible problems --------------

# Next steps -----------
## a bit too many features so will need to run dimentionality reduction algorithms (eg. PCA - principal component analysis)



plot(my_data$)

