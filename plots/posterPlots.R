library(tidyverse)
library(readr)
library(mgcv)
library(pracma)

alpha01O <- read_csv("~eeg/raw_data/SPUR-EEG-data/CTEEG001O_Alpha.csv")

alpha01O <- normalise_eeg("~/eeg/raw_data/SPUR-EEG-data/CTEEG001O_Alpha.csv")
gamma01O <- normalise_eeg("~/eeg/raw_data/SPUR-EEG-data/CTEEG001O_Gamma1.csv")
gamma033Y <- normalise_eeg("~/eeg/raw_data/SPUR-EEG-data/CTEEG033Y_Gamma1.csv")

plot(gamma033Y$Time[200:400], gamma033Y$F3[200:400], type = "l", col = "blue", 
     xlab = "Time (s)", ylab = "EEG Signal (uV)", 
     main = "EEG Signal for CTEEG033Y - Gamma Band")

norm_x <- as.numeric(scale(gamma033Y$F3))
# GAM smooth
fit      <- gam(norm_x ~ s(gamma033Y$Time, k = 70))
smooth_x <- as.numeric(predict(fit, data.frame(time = gamma033Y$Time)))

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