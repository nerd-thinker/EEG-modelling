library(dtw)
library(ggplot2)
# Lookup table for each wave ----
## why we need restructured data frames: instead of channels repeating across time, lookup collapses it to one row per channel
## with the entire time series packed into a list 
## also: DTW needs the whole time series for a channel as a single vector to compare against another channel's vector
make_smooth_lookup <- function(df) {
  df %>%
    group_by(Channel) %>%
    summarise(
      Time           = list(Time),
      smoothed_power = list(as.numeric(smoothed_power)),
      .groups = "drop"
    )
}

theta_lookup <- make_smooth_lookup(theta_results)
alpha_lookup <- make_smooth_lookup(alpha_results)
beta_lookup  <- make_smooth_lookup(beta_results)
delta_lookup <- make_smooth_lookup(delta_results)
gamma_lookup <- make_smooth_lookup(gamma_results)

# Downsampling ----
downsample_lookup <- function(lookup, every_n = 8) {
  lookup %>%
    mutate(
      Time           = map(Time,           ~ .x[seq(1, length(.x), by = every_n)]),
      smoothed_power = map(smoothed_power, ~ .x[seq(1, length(.x), by = every_n)])
    )
}

theta_lookup_ds  <- downsample_lookup(theta_lookup)
alpha_lookup_ds  <- downsample_lookup(alpha_lookup)
beta_lookup_ds   <- downsample_lookup(beta_lookup)
delta_lookup_ds  <- downsample_lookup(delta_lookup)
gamma_lookup_ds  <- downsample_lookup(gamma_lookup)

## check how much it changed: 
raw_gam <- alpha_results %>% filter(Channel == "Fp1") %>% pull(smoothed_power)
raw_time_gam <- alpha_results %>% filter(Channel == "Fp1") %>% pull(Time)

ds_signal <- alpha_lookup_ds$smoothed_power[[which(alpha_lookup_ds$Channel == "Fp1")]]
ds_time   <- alpha_lookup_ds$Time[[which(alpha_lookup_ds$Channel == "Fp1")]]

# Plot
plot(x = raw_time_gam, y = raw_gam, type = "l", lwd = 1,
     xlab = "Time", ylab = "Smoothed Power",
     main = "Alpha Fp1 - Raw vs Downsampled")

lines(x = ds_time, y = ds_signal, col = 2, lwd = 2)

legend("topright", legend = c("Raw", "Downsampled (every 8th)"),
       col = c(1, 2), lty = 1, lwd = c(1, 2))

# 32 x 32 channel matrices (5 waves) ----
compute_channel_dtw <- function(lookup) {
  channels <- lookup$Channel
  n <- length(channels)
  mat <- matrix(0, n, n, dimnames = list(channels, channels))
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d <- dtw(
        lookup$smoothed_power[[i]],
        lookup$smoothed_power[[j]],
        distance.only = TRUE
      )$normalizedDistance
      mat[i, j] <- d
      mat[j, i] <- d
    }
  }
  mat
}

theta_channel_mat <- compute_channel_dtw(theta_lookup_ds)
cat("theta done\n")
alpha_channel_mat <- compute_channel_dtw(alpha_lookup_ds)
cat("alpha done\n")
beta_channel_mat  <- compute_channel_dtw(beta_lookup_ds)
cat("beta done\n")
delta_channel_mat <- compute_channel_dtw(delta_lookup_ds)
cat("delta done\n")
gamma_channel_mat <- compute_channel_dtw(gamma_lookup_ds)
cat("gamma done\n")

# boxplot TIME -----
boxplot(theta_channel_mat)
