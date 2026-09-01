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
peaks        <- findpeaks(smooth_x,
                           minpeakheight   = mean(smooth_x) + sd(smooth_x),
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

y_range <- range(c(plot_df$raw, plot_df$smooth))

# Palette (dataviz skill default: categorical slot 1 blue, slot 2 orange,
# status/critical red used as a highlight marker, not a series)
col_raw    <- "#2a78d6"
col_smooth <- "#eb6834"
col_peak   <- "#d03b3b"

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
  geom_line(color = col_raw, linewidth = 0.6) +
  coord_cartesian(ylim = y_range) +
  labs(title = "1. Raw EEG Signal", x = "Time (s)", y = "Amplitude (z-scored)") +
  poster_theme

# Panel 2: raw + GAM smooth + detected peaks
p2 <- ggplot(plot_df, aes(time, raw)) +
  geom_line(aes(color = "Raw signal"), linewidth = 0.5, alpha = 0.55) +
  geom_line(aes(time, smooth, color = "GAM smooth (k = 70)"), linewidth = 1.3) +
  geom_vline(data = peak_df, aes(xintercept = time), linetype = "dashed",
             color = "grey55", linewidth = 0.5) +
  geom_point(data = peak_df, aes(time, value, color = "Detected peak"), size = 4.5) +
  coord_cartesian(ylim = y_range) +
  scale_color_manual(values = c(
    "Raw signal"          = col_raw,
    "GAM smooth (k = 70)" = col_smooth,
    "Detected peak"       = col_peak
  ), breaks = c("Raw signal", "GAM smooth (k = 70)", "Detected peak")) +
  labs(title = "2. GAM Smoothing & Peak Detection", x = "Time (s)", y = "Amplitude (z-scored)") +
  guides(color = guide_legend(override.aes = list(
    linetype  = c("solid", "solid", "blank"),
    shape     = c(NA, NA, 16),
    linewidth = c(1, 1.3, NA)
  ))) +
  poster_theme

combined <- p1 + p2 +
  plot_annotation(
    title = "Peak Detection in EEG Signal",
    theme = theme(plot.title = element_text(face = "bold", size = 26, hjust = 0.5,
                                             margin = margin(b = 6)))
  )

ggsave("plots/poster_peak_detection.png", combined, width = 16, height = 6.5, dpi = 300, bg = "white")
