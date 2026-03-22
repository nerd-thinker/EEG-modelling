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

# boxplot time --------
#### COMMENT: THE SEQUENCE OF NODE IS NOT CORRECT NEED TO FIX IT
boxplot(alpha_channel_mat, main = "Alpha gam's dtw distances")
boxplot(beta_channel_mat, main = "Betta gam's dtw distances")
boxplot(delta_channel_mat, main = "Delta gam's dtw distances")
boxplot(gamma_channel_mat, main = "Gamma gam's dtw distances")
boxplot(theta_channel_mat, main = "Theta gam's dtw distances")


# 5x5 matrices ---
wave_lookups <- list(
  theta = theta_lookup_ds,
  alpha = alpha_lookup_ds,
  beta  = beta_lookup_ds,
  delta = delta_lookup_ds,
  gamma = gamma_lookup_ds
)

compute_wave_dtw <- function(channel_name, wave_lookups) {
  wave_names <- names(wave_lookups)
  n <- length(wave_names)
  mat <- matrix(0, n, n, dimnames = list(wave_names, wave_names))
  
  # Extract this channel's time series from each wave
  ts_list <- lapply(wave_lookups, function(lookup) {
    idx <- which(lookup$Channel == channel_name)
    lookup$smoothed_power[[idx]]
  })
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d <- dtw(ts_list[[i]], ts_list[[j]],
               distance.only = TRUE)$normalizedDistance
      mat[i, j] <- d
      mat[j, i] <- d
    }
  }
  mat
}

# Run for all 32 channels
all_channels <- theta_lookup_ds$Channel

wave_dtw_matrices <- lapply(
  setNames(all_channels, all_channels),
  compute_wave_dtw,
  wave_lookups = wave_lookups
)

cat("All 5x5 matrices done\n")

# Should be a 5x5 symmetric matrix with 0 diagonal
x <-wave_dtw_matrices[["Fp1"]]
boxplot(x)
isSymmetric(wave_dtw_matrices[["Fp1"]])  # TRUE

# Combine all channels into one long data frame
wave_long <- lapply(all_channels, function(ch) {
  wave_dtw_matrices[[ch]] %>%
    as.data.frame() %>%
    rownames_to_column("Wave_1") %>%
    pivot_longer(-Wave_1, names_to = "Wave_2", values_to = "DTW_distance") %>%
    filter(Wave_1 != Wave_2) %>%
    mutate(Channel = ch)
}) %>% bind_rows()

# Boxplot: distribution of wave-to-wave DTW per channel
ggplot(wave_long, aes(x = Channel, y = DTW_distance)) +
  geom_boxplot(fill = "coral", alpha = 0.7, outlier.size = 0.8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Wave DTW Distances per Channel (5x5)",
       x = "Channel", y = "Normalised DTW Distance")
