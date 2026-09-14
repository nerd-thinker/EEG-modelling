library(tidyverse)
library(readr)
library(mgcv)
library(pracma)
library(patchwork)

# ---------------------------------------------------------------------------
# Poster figure: how peaks are detected in the EEG signal.
# Panel 1 - raw signal. Panel 2 - raw signal + GAM smooth (k = 70) + peaks.
#
# NOTE on data choice: gamma033Y / channel F3 (the original draft here) only
# produces 1 detected peak across the whole ~197s recording, so it can't
# illustrate the method. alpha01O / channel P7 has 11 well-separated peaks
# and is used instead for a clear demo figure.
# ---------------------------------------------------------------------------

df <- read_csv("~/eeg/raw_data/SPUR-EEG-data/CTEEG001O_Alpha.csv", show_col_types = FALSE)
df <- df[, !grepl("^\\.\\.\\.\\d+$", colnames(df)), drop = FALSE]  # drop trailing-comma phantom col

t      <- df$Time
norm_x <- as.numeric(scale(df$P7))

# GAM smooth
fit      <- gam(norm_x ~ s(t, k = 70))
smooth_x <- as.numeric(predict(fit, newdata = data.frame(t = t)))

# Peak detection on smoothed signal
mean_smooth <- mean(smooth_x)
sd_smooth   <- sd(smooth_x)
threshold   <- mean_smooth + sd_smooth

peaks        <- findpeaks(smooth_x,
                           minpeakheight   = threshold,
                           minpeakdistance = 10)
peak_indices <- peaks[, 2]

# Zoom window for the poster panels (4 peaks, clearly separated)
win <- which(t >= 55 & t <= 95)

plot_df <- tibble(
  time   = t[win],
  raw    = norm_x[win],
  smooth = smooth_x[win]
)
peak_df <- tibble(
  time  = t[peak_indices[peak_indices %in% win]],
  value = smooth_x[peak_indices[peak_indices %in% win]]
)

# Threshold lines (drawn instead of per-peak vertical guides): the mean of
# the smoothed signal, and the mean + 1 SD cutoff actually used by findpeaks()
# above to decide what counts as a peak.
threshold_df <- tibble(
  type  = factor(c("Mean", "Mean + 1 SD (peak threshold)"),
                 levels = c("Mean", "Mean + 1 SD (peak threshold)")),
  value = c(mean_smooth, threshold)
)

y_range <- range(c(plot_df$raw, plot_df$smooth, threshold_df$value))

# Poster palette: poppy, high-contrast red/blue for the two data series;
# black for the peak markers so they read clearly against both.
col_raw    <- "#0033ff"
col_smooth <- "#e8112d"
col_peak   <- "#111111"
col_mean   <- "grey30"
col_sd     <- "grey65"

poster_theme <- theme_minimal(base_size = 20) +
  theme(
    plot.title       = element_text(face = "bold", size = 22, margin = margin(b = 4)),
    axis.title       = element_text(color = "grey25", size = 17),
    axis.text        = element_text(color = "grey40", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 15),
    legend.title     = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(10, 16, 10, 10)
  )

# Panel 1: raw signal only
p1 <- ggplot(plot_df, aes(time, raw)) +
  geom_line(color = col_raw, linewidth = 0.75) +
  coord_cartesian(ylim = y_range) +
  labs(title = "1. Raw EEG Signal", x = "Time (s)", y = "Amplitude (z-scored)") +
  poster_theme

# Panel 2: raw + GAM smooth + mean / +1 SD threshold + detected peaks
p2 <- ggplot(plot_df, aes(time, raw)) +
  geom_hline(data = subset(threshold_df, type == "Mean"),
             aes(yintercept = value, color = type), linewidth = 0.8) +
  geom_hline(data = subset(threshold_df, type == "Mean + 1 SD (peak threshold)"),
             aes(yintercept = value, color = type), linewidth = 0.6,
             linetype = "dashed", alpha = 0.8) +
  geom_line(aes(color = "Raw signal"), linewidth = 0.8, alpha = 0.9) +
  geom_line(aes(time, smooth, color = "GAM smooth (k = 70)"), linewidth = 1.3) +
  geom_point(data = peak_df, aes(time, value, color = "Detected peak"), size = 4.5) +
  coord_cartesian(ylim = y_range) +
  scale_color_manual(values = c(
    "Raw signal"                    = col_raw,
    "GAM smooth (k = 70)"           = col_smooth,
    "Mean"                          = col_mean,
    "Mean + 1 SD (peak threshold)"  = col_sd,
    "Detected peak"                 = col_peak
  ), breaks = c("Raw signal", "GAM smooth (k = 70)", "Mean",
                "Mean + 1 SD (peak threshold)", "Detected peak")) +
  labs(title = "2. GAM Smoothing & Peak Detection", x = "Time (s)", y = "Amplitude (z-scored)") +
  guides(color = guide_legend(override.aes = list(
    linetype  = c("solid", "solid", "solid", "dashed", "blank"),
    shape     = c(NA, NA, NA, NA, 16),
    linewidth = c(1, 1.3, 0.8, 0.6, NA)
  ))) +
  poster_theme

combined <- p1 + p2 +
  plot_annotation(
    title = "Peak Detection in EEG Signal",
    theme = theme(plot.title = element_text(face = "bold", size = 26, hjust = 0.5,
                                             margin = margin(b = 6)))
  )

ggsave("plots/poster_peak_detection.png", combined, width = 16, height = 6.5, dpi = 300, bg = "white")
