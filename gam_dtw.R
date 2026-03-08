## run all waves_results to get smoothed_power col of GAM models

## reshape data to wide format (cahnels x Time)

library(tidyverse)
library(dtw)
#library(dtwclust)  # For faster DTW

# Reshape from long to wide format
theta_wide <- theta_results %>%
  select(Channel, Time, smoothed_power) %>%
  pivot_wider(
    id_cols = Time,
    names_from = Channel,
    values_from = smoothed_power
  ) %>%
  arrange(Time)  # Make sure time is ordered

# Check the result
dim(theta_wide)  # Should be [timepoints × 33] (Time + 32 channels)
theta_wide[1:5, 1:5]  # View first few rows and columns

# Extract just the channel data as a matrix
# Rows = channels, Columns = time points
channel_matrix <- theta_wide %>%
  select(-Time) %>%  # Remove Time column
  as.matrix() %>%
  t()  # Transpose so channels are rows, time points are columns

# Get channel names
channels <- colnames(theta_wide)[-1]  # Exclude Time column

# Verify dimensions
dim(channel_matrix)  # Should be 32 × timepoints
rownames(channel_matrix) <- channels

# Basic DTW computation (slower but clear)
compute_dtw_matrix_base <- function(channel_matrix) {
  n_channels <- nrow(channel_matrix)
  channels <- rownames(channel_matrix)
  
  # Initialize 32×32 matrix
  dtw_matrix <- matrix(0, nrow = n_channels, ncol = n_channels,
                       dimnames = list(channels, channels))
  
  # Compute pairwise DTW distances
  for(i in 1:(n_channels-1)) {
    for(j in (i+1):n_channels) {
      # Get time series for channels i and j
      ts_i <- channel_matrix[i, ]
      ts_j <- channel_matrix[j, ]
      
      # Compute DTW distance
      dtw_dist <- dtw(ts_i, ts_j, 
                      #window.type = "sakoeChiba",  # Optional: limit warping
                      window.size = 10)$distance   # Adjust window size
      
      # Fill symmetric matrix
      dtw_matrix[i, j] <- dtw_dist
      dtw_matrix[j, i] <- dtw_dist
    }
  }
  
  return(dtw_matrix)
}

# Compute (this may take a few minutes)
theta_dtw_base <- compute_dtw_matrix_base(channel_matrix)

