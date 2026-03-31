library(dplyr)
library(tidyr)
library(mgcv)

## my experimentation with epochs: 
x <- bands_clean$alpha$Fp1

plot(x, type = "l", col = "blue")
model <- gam(x~ s(Time, k= 60), data = bands_clean$alpha)
prediction <- predict(model)
lines(prediction, col = "red", lwd = 2)

epochs <- seq(from =0, to = length(x), by = 200) ## 200 seems good enough to determine if the gam function is increasing or decreasing
abline(v=epochs)

# feature function ------------

#epochs no can change if you want more defined information during smaller interval

time_domain_features <- function(band_df, epochs = 200, k = 60) {
  
  # Convert to long format: Time | channel | signal
  long_df <- band_df %>%
    pivot_longer(
      cols = -Time,
      names_to = "channel",
      values_to = "signal"
    )
  
  # Create epochs
  long_df <- long_df %>%
    mutate(epoch = cut(Time, breaks = epochs, labels = FALSE))
  
  # GAM smoothing per channel
  long_df <- long_df %>%
    group_by(channel) %>%
    group_modify(~ {
      fit <- gam(signal ~ s(Time, k = k), data = .x)
      .x$smooth <- predict(fit, newdata = .x)
      .x
    }) %>%
    ungroup()
  
  # Raw features + Hjorth
  features_raw <- long_df %>%
    group_by(channel, epoch) %>%
    summarise(
      mean_power     = mean(signal, na.rm = TRUE),
      var_power      = var(signal,  na.rm = TRUE),
      sd_power       = sd(signal,   na.rm = TRUE),
      max_power      = max(signal,  na.rm = TRUE),
      min_power      = min(signal,  na.rm = TRUE),
      range_power    = max_power - min_power,
      activity       = hjorth(signal)$activity,
      mobility       = hjorth(signal)$mobility,
      complexity     = hjorth(signal)$complexity,
      .groups = "drop"
    )
  
  # Smoothed features + Hjorth
  features_smooth <- long_df %>%
    group_by(channel, epoch) %>%
    summarise(
      mean_power     = mean(smooth, na.rm = TRUE),
      var_power      = var(smooth,  na.rm = TRUE),
      sd_power       = sd(smooth,   na.rm = TRUE),
      max_power      = max(smooth,  na.rm = TRUE),
      min_power      = min(smooth,  na.rm = TRUE),
      range_power    = max_power - min_power,
      activity       = hjorth(smooth)$activity,
      mobility       = hjorth(smooth)$mobility,
      complexity     = hjorth(smooth)$complexity,
      .groups = "drop"
    )
  
  return(list(raw = features_raw, smooth = features_smooth))
}
# Define Hjorth parameters from scratch
hjorth <- function(x) {
  x <- na.omit(x)
  
  dx  <- diff(x)        # first derivative
  ddx <- diff(dx)       # second derivative
  
  activity   <- var(x)
  mobility   <- sqrt(var(dx) / var(x)) ## estimate of the mean frequency
  complexity <- (sqrt(var(ddx) / var(dx))) / mobility ## estimate of the bandwidth of the signal, which indicates the similarity of the shape of the signal to a pure sine wave
  
  list(activity = activity, mobility = mobility, complexity = complexity)
}

# Test it
test <- c(1.2, 1.5, 1.3, 1.8, 1.1, 1.6)
hjorth(test)

## Compute the basic features for separate wavelengths
alpha_basic_features <- time_domain_features(bands_clean$alpha)
beta_basic_features <- time_domain_features(bands_clean$beta)
delta_basic_fetures <- time_domain_features(bands_celan$delta)
theta_basic_features <- time_domain_features(bands_clean$theta)
gamma_basic_features <- time_domain_features(bands_clean$gamma)
#for all:
wave_features <- lapply(bands_clean, time_domain_features)

## compare results:
alpha_basic_features$raw %>%
  
alpha_basic_features$smooth


### hjorth parameters interpretation
hjorth(bands_clean$alpha$Fp1)
hjorth(bands_norm$alpha$Fp1) ##be careful when using normalized as activity = 1 (var of whole time-series)

##gamma should have higher mobility than alpha:
hjorth(bands_clean$gamma$Fp1) ##mobility: 0.151 > 0.043 as expected


