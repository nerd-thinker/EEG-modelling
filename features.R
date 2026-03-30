## my experimentation got deleted unfortunetly, how sad

# feature function -------------
library(dplyr)
library(tidyr)
library(mgcv)

time_domain_features <- function(band_df, epochs = 200) { ## maybe add k factor for gams
  
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
      fit <- gam(signal ~ s(Time, k = 60), data = .x)
      .x$smooth <- predict(fit, newdata = .x)
      .x
    }) %>%
    ungroup()
  
  #  Raw features
  features_raw <- long_df %>%
    group_by(channel, epoch) %>%
    summarise(
      mean_power  = mean(signal, na.rm = TRUE),
      var_power   = var(signal, na.rm = TRUE),
      sd_power    = sd(signal, na.rm = TRUE),
      max_power   = max(signal, na.rm = TRUE),
      min_power   = min(signal, na.rm = TRUE),
      range_power = max_power - min_power,
      .groups = "drop"
    )
  
  #  Smoothed features
  features_smooth <- long_df %>%
    group_by(channel, epoch) %>%
    summarise(
      mean_power  = mean(smooth, na.rm = TRUE),
      var_power   = var(smooth, na.rm = TRUE),
      sd_power    = sd(smooth, na.rm = TRUE),
      max_power   = max(smooth, na.rm = TRUE),
      min_power   = min(smooth, na.rm = TRUE),
      range_power = max_power - min_power,
      .groups = "drop"
    )
  
  return(list(raw = features_raw, smooth = features_smooth))
}