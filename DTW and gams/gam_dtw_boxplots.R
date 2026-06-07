# requirements to run this:
## run whole script of gam_dtw 
## !!! run all theta_channel_mat <- compute_channel_dtw(theta_lookup_ds) (from gam_dtw)

# 32x32 boxplots across nodes -----
# Define the correct channel order and grouping ----
channel_order <- c(
  "Fp1", "Fp2",
  "AF3", "AF4",
  "F7", "F3", "Fz", "F4", "F8",
  "FC5", "FC1", "FC2", "FC6",
  "T7", "C3", "Cz", "C4", "T8",
  "CP5", "CP1", "CP2", "CP6",
  "P7", "P3", "Pz", "P4", "P8",
  "PO3", "PO4",
  "O1", "Oz", "O2"
)

# Define region groups and assign a colour to each ----
channel_regions <- list(
  Prefrontal  = c("Fp1", "Fp2"),
  AntFrontal  = c("AF3", "AF4"),
  Frontal     = c("F7", "F3", "Fz", "F4", "F8"),
  FrontoCent  = c("FC5", "FC1", "FC2", "FC6"),
  Temporal    = c("T7", "C3", "Cz", "C4", "T8"),
  CentPar     = c("CP5", "CP1", "CP2", "CP6"),
  Parietal    = c("P7", "P3", "Pz", "P4", "P8"),
  ParOccip    = c("PO3", "PO4"),
  Occipital   = c("O1", "Oz", "O2")
)

region_colours <- c(
  Prefrontal  = "#E41A1C",
  AntFrontal  = "#FF7F00",
  Frontal     = "#F0E442",
  FrontoCent  = "#4DAF4A",
  Temporal    = "#00BFC4",
  CentPar     = "#377EB8",
  Parietal    = "#984EA3",
  ParOccip    = "#A65628",
  Occipital   = "#F781BF"
)

# Map each channel to its region colour ----
channel_colours <- setNames(
  sapply(channel_order, function(ch) {
    region_colours[names(which(sapply(channel_regions, function(r) ch %in% r)))]
  }),
  channel_order
)

# Helper: reorder columns of a DTW matrix to match channel_order ----
reorder_mat <- function(mat) {
  ord <- intersect(channel_order, colnames(mat))   # keeps only channels present
  mat[ord, ord]
}

# Base-R boxplots (32x32 matrices) ----
plot_dtw_boxplot <- function(mat, title) {
  mat_ord <- reorder_mat(mat)
  boxplot(mat_ord,
          main   = title,
          col    = channel_colours[colnames(mat_ord)],
          las    = 2,          # rotate x-axis labels
          cex.axis = 0.7)
  
  # Add a legend for regions
  #legend("topright",
         #legend = names(region_colours),
         #fill   = region_colours,
         #cex    = 0.65,
         #bty    = "n")
}

plot_dtw_boxplot(alpha_channel_mat, "Alpha DTW distances")
plot_dtw_boxplot(beta_channel_mat,  "Beta DTW distances")
plot_dtw_boxplot(delta_channel_mat, "Delta DTW distances")
plot_dtw_boxplot(gamma_channel_mat, "Gamma DTW distances")
plot_dtw_boxplot(theta_channel_mat, "Theta DTW distances")

# ggplot boxplot for all 5 waves, colour-coded by region ----

# Stack all 5 matrices into one long data frame
wave_names <- c("alpha", "beta", "delta", "gamma", "theta")

wave_mats <- list(
  alpha = alpha_channel_mat,
  beta  = beta_channel_mat,
  delta = delta_channel_mat,
  gamma = gamma_channel_mat,
  theta = theta_channel_mat
)

all_waves_long <- lapply(names(wave_mats), function(wv) {
  mat <- wave_mats[[wv]]
  mat_ord <- reorder_mat(mat)           # reorder columns to channel_order
  as.data.frame(mat_ord) %>%
    pivot_longer(everything(),
                 names_to  = "Channel",
                 values_to = "DTW_distance") %>%
    filter(DTW_distance > 0) %>%        # drop self-comparisons (diagonal zeros)
    mutate(Wave = wv)
}) %>%
  bind_rows() %>%
  mutate(
    Channel = factor(Channel,
                     levels = intersect(channel_order, unique(Channel))),
    Region  = channel_to_region[as.character(Channel)],
    Region  = factor(Region, levels = names(channel_regions)),
    Wave    = factor(Wave, levels = wave_names)
  )

# Plot
ggplot(all_waves_long, aes(x = Channel, y = DTW_distance, fill = Region)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
  scale_fill_manual(values = region_colours) +
  facet_wrap(~ Wave, ncol = 1, strip.position = "right") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    strip.text       = element_text(size = 9, face = "bold"),
    legend.position  = "top",
    legend.key.size  = unit(0.5, "cm"),
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(0.3, "lines")
  ) +
  labs(
    title = "Channel DTW distances by wave and region",
    x     = NULL,
    y     = "Normalised DTW distance",
    fill  = "Region"
  )
# Separate boxplot per wave ----
plot_wave_boxplot <- function(wave_name) {
  df <- all_waves_long %>% filter(Wave == wave_name)
  
  ggplot(df, aes(x = Channel, y = DTW_distance, fill = Region)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
    scale_fill_manual(values = region_colours) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "top",
      legend.key.size = unit(0.5, "cm"),
      legend.text     = element_text(size = 8),
      legend.title    = element_text(size = 9),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste(tools::toTitleCase(wave_name), "wave — channel DTW distances"),
      x     = NULL,
      y     = "Normalised DTW distance",
      fill  = "Region"
    )
}

# Print all 5 individually
lapply(levels(all_waves_long$Wave), plot_wave_boxplot)
