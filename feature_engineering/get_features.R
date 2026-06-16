library(mgcv)
library(tidyverse)

get_features <- function(x,
                         normalisation = c("none","z_score","01"),
                         k = 20) {
  time <- 1:length(x)
  normalisation <- match.arg(normalisation)
  norm_x <- switch(normalisation,
                   "none" = x,
                   "z_score" = as.numeric(scale(x)),
                   "01" = (x - min(x))/(max(x) - min(x)))
  fit <- gam(norm_x ~ s(time, k = k))
  smooth_x <- predict(fit, data.frame(time = seq(1, max(time), length = 1000)))
  peaks <- detect_peaks(smooth_x)
  
  plot(1:1000, smooth_x, col = 2, ty = "l")
  abline(v = peaks, lty = 2)
  
  return(c("mean" = mean(smooth_x),
           "var" = var(smooth_x),
           "q25" = quantile(smooth_x, .25),
           "q75" = quantile(smooth_x, .75),
           "n_peaks" = length(peaks),
           "mean_peak_height" = mean(smooth_x[peaks]),
           "var_peak_height" = var(smooth_x[peaks]),
           "q25_peak_height" = quantile(smooth_x[peaks], .25),
           "q75_peak_height" = quantile(smooth_x[peaks], .75)))
}

detect_peaks <- function(x, min_prominence = 0, allow_plateau = TRUE) {
  n <- length(x)
  peaks <- rep(FALSE, n)
  
  for (i in 2:(n - 1)) {
    
    # Strict peak
    if (x[i] > x[i - 1] && x[i] > x[i + 1]) {
      peaks[i] <- TRUE
    }
    
    # Plateau peak
    if (allow_plateau && x[i] >= x[i - 1] && x[i] >= x[i + 1]) {
      
      # expand plateau
      left <- i
      right <- i
      
      while (left > 1 && x[left] == x[i]) left <- left - 1
      while (right < n && x[right] == x[i]) right <- right + 1
      
      # check edges of plateau
      if (x[left] < x[i] && x[right] < x[i]) {
        peaks[i] <- TRUE
      }
    }
  }
  
  # Optional: filter by prominence
  if (min_prominence > 0) {
    for (i in which(peaks)) {
      left_min <- min(x[1:i])
      right_min <- min(x[i:n])
      prominence <- x[i] - max(left_min, right_min)
      
      if (prominence < min_prominence) {
        peaks[i] <- FALSE
      }
    }
  }
  
  return(which(peaks))
}

my_data <- tibble(delta_fp1 = sin(seq(0,6*pi, length=500)) + rnorm(500, 0, .1),
                  delta_fp2 = cos(seq(0,6*pi, length=500)) + rnorm(500, 0, .1))

get_features(my_data$delta_fp1)
get_features(my_data$delta_fp2)

apply(my_data, 2, get_features) ## this gives us a matrix of features for 1 single participant = every cell is a feature

# we will have
## 32 nodes * 5 waves = 160 node:wave combinations
## 3 normalisations
## X pre-features
## total: 160 * 3 * X per participant


