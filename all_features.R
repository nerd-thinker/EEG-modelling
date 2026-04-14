# libraries
library(pracma)
library(dplyr)
library(tidyr)
library(mgcv)
library(moments)

# combined function for all
get_features <- function(x,
                         time = NULL,
                         normalisation = c("none","z_score","01"),
                         k = 20) {
  if(is.null(time)){ 
    time <- seq_along(x)
  } ## time needs to be provided as real time in seconds node$Time != sample indices (1:length(x))
  normalisation <- match.arg(normalisation)
  norm_x <- switch(normalisation,
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
  plot(time, norm_x, type = "l",
       xlab = "Time", ylab = "Signal",
       col = "grey87")                              # raw signal in grey
  lines(time, smooth_x, col = 2, lwd = 2)   # GAM smooth in red
  
  if (length(peak_indices) > 0)
    abline(v = time[peak_indices], lty = 2, col = "grey20")
  
  return(c(
    # General smooth features
    "mean"                    = mean(smooth_x),
    #"var"                     = var(smooth_x),
    "q25"                     = unname(quantile(smooth_x, .25)),
    "q75"                     = unname(quantile(smooth_x, .75)),
    "skewdness"               = moments::skewness(smooth_x),
    "kurtosis"                = moments::kurtosis(smooth_x),
    
    # Peak features
    "n_peaks"                 = num_peaks,
    "peak_frequency"          = peak_freq, ##have if time for different participants is nto equal
    "mean_peak_amplitude"     = mean_amp,
    "max_peak_amplitude"      = max_amp,
    "sd_peak_amplitude"       = sd_amp,
    "mean_inter_peak_interval"= mean_ipi,
    "sd_inter_peak_interval"  = sd_ipi,
    
    # Hjorth
    "activity"                = hjorth(smooth_x)$activity, ##same as variance
    "mobility"                = hjorth(smooth_x)$mobility,
    "complexity"              = hjorth(smooth_x)$complexity
  ))
}
hjorth <- function(x) {
  x <- na.omit(x)
  
  dx  <- diff(x)        # first derivative
  ddx <- diff(dx)       # second derivative
  
  activity   <- var(x)
  mobility   <- sqrt(var(dx) / var(x)) ## estimate of the mean frequency
  complexity <- (sqrt(var(ddx) / var(dx))) / mobility ## estimate of the bandwidth of the signal, which indicates the similarity of the shape of the signal to a pure sine wave
  
  list(activity = activity, mobility = mobility, complexity = complexity)
}

# Test the functions
node <- bands_norm$alpha$FC5
time <- bands_norm$alpha$Time
get_features(node,time = time, normalisation = "z_score")
get_features(node,time = time, normalisation = "01")
get_features(node, normalisation = "none")

smooth_ds <- node[seq(1,length(node), by = 8)]
pracma::sample_entropy(smooth_ds, edim = 2, r = 0.2 * sd(smooth_ds))
