library(dplyr)
library(tidyr)
library(mgcv)

## my experimentation got deleted unfortunetly, how sad
x <- bands_clean$alpha$Fp1

plot(x, type = "l", col = "blue")
model <- gam(x~ s(Time, k= 60), data = bands_clean$alpha)
prediction <- predict(model)
lines(prediction, col = "red", lwd = 2)

epochs <- seq(from =0, to = length(x), by = 200) ## 200 seems good enough to determine if the gam function is increasing or decreasing
abline(v=epochs)

# feature function ------------

#epochs no can change if you want more defined information during smaller interval
library(pracma) #library for hjorth parameter estimates

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
