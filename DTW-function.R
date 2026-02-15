
# Preparing data ----------------------------------------------------------

#libraries
library(dtw)
library(readxl)
library(proxy)
library(dplyr)
library(tidyr)

# reading data in
Delta <- read_excel("eeg/raw_data/Delta.xlsx")
Alpha <- read_excel("eeg/raw_data/Alpha.xlsx")
Beta <- read_excel("eeg/raw_data/Beta.xlsx")
Gamma <- read_excel("eeg/raw_data/Gamma1.xlsx")
Theta <- read_excel("eeg/raw_data/Theta.xlsx")

# bands 
bands <- list(
  alpha = Alpha,
  delta = Delta,
  theta = Theta,
  beta  = Beta,
  gamma = Gamma
)

#clean time
bands_clean <- lapply(bands, function(df) {
  df[df$Time >= 0, ]
})

# resetting row names
bands_clean <- lapply(bands_clean, function(df) {
  rownames(df) <- NULL
  df
})

## normalize function for all bands
normalize_band <- function(df) {
  df_norm <- df
  eeg_cols <- setdiff(names(df), "Time")
  
  df_norm[eeg_cols] <- lapply(df[eeg_cols], function(x) {
    as.numeric(scale(x))
  })
  
  df_norm
}

## apply normalization to all bands
bands_norm <- lapply(bands_clean, normalize_band)

## downsample function
downsample <- function(x, k) {
  x[seq(1, length(x), by = k)]
}

## dtw-distance function
dtw_dist <- function(x, y, window_frac = 0.1) {
  w <- floor(window_frac * min(length(x), length(y)))
  dtw(x, y, window.type = "sakoechiba", window.size = w)$distance
}



# DTW for 5x5 each electrode matrix ---------------------------------------------------------------------

# DTW matrix for one electrode:
compute_dtw_matrix <- function(bands_norm, electrode, k = 5) {
  
  ## extract same electrode across bands
  signals <- lapply(bands_norm, function(df) df[[electrode]])
  
  ## optional sanity check
  # print(sapply(signals, sd))
  
  ## downsample
  signals_ds <- lapply(signals, downsample, k = k)
  
  band_names <- names(signals_ds)
  n <- length(signals_ds)
  
  ## initialize matrix
  dtw_matrix <- matrix(
    NA_real_, n, n,
    dimnames = list(band_names, band_names)
  )
  
  ## DTW computation
  for (i in seq_len(n)) {
    dtw_matrix[i, i] <- 0
    for (j in i:n) {
      if (i != j) {
        d <- dtw_dist(signals_ds[[i]], signals_ds[[j]])
        dtw_matrix[i, j] <- d
        dtw_matrix[j, i] <- d
      }
    }
  }
  
  return(dtw_matrix)
}

# run dtw for all electrodes automatically
electrodes <- setdiff(colnames(Alpha), "Time")

dtw_results <- vector("list", length(electrodes))
names(dtw_results) <- electrodes

for (e in electrodes) {
  message("Processing ", e)
  dtw_results[[e]] <- compute_dtw_matrix(
    bands_norm = bands_norm,
    electrode = e,
    k = 5
  )
  gc(FALSE)
}

dtw_results$Fp2

#FP1 Fp2 | AF3 AF4 | F7 F3 Fz F4 F8 | FC5 FC1 FC2 FC6 | T7 C3 Cz C4 T8 | CP5 CP1 CP2 CP6 | P7 P3 Pz P4 P8 | PO3 PO4 | O1 Oz O2 

mah <- dtw_results$O2 ##edit for specific electrode

# Basic boxplot
boxplot(mah,
        main = "Boxplot of O2 Electrode Data",
        ylab = "Value (μV or similar unit)",
        col = "lightblue",
        border = "darkblue")

# Add mean point
#points(mean(mah), col = "red", pch = 19, cex = 1.2)



# DTW 32x32 matrix each wavelength ----------------------------------------
# frequency specific pattern
compute_band_electrode_dtw <- function(df, k = 5) {
  electrodes <- setdiff(colnames(df), "Time")
  signals <- lapply(df[electrodes], downsample, k = k)
  
  n <- length(signals)
  mat <- matrix(NA_real_, n, n,
                dimnames = list(electrodes, electrodes))
  
  for (i in seq_len(n)) {
    mat[i, i] <- 0
    for (j in i:n) {
      if (i != j) {
        d <- dtw_dist(signals[[i]], signals[[j]])
        mat[i, j] <- d
        mat[j, i] <- d
      }
    }
  }
  mat
}

dtw_by_band <- lapply(bands_norm, compute_band_electrode_dtw)
dtw_by_band$alpha

## Reordering electrodes in order:: 
electrode_order <- c(
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

#reorder a single matrix
reorder_dtw_matrix <- function(mat, electrode_order) {
  mat[electrode_order, electrode_order]
}
#reorder all matrixes
dtw_by_band_ordered <- lapply(
  dtw_by_band,
  reorder_dtw_matrix,
  electrode_order = electrode_order
)

#reordering band_norm:
alpha_norm_alphb<-bands_norm$alpha%>% select("Time", "Fp1", "Fp2",
                                             "AF3", "AF4",
                                             "F7", "F3", "Fz", "F4", "F8",
                                             "FC5", "FC1", "FC2", "FC6",
                                             "T7", "C3", "Cz", "C4", "T8",
                                             "CP5", "CP1", "CP2", "CP6",
                                             "P7", "P3", "Pz", "P4", "P8",
                                             "PO3", "PO4",
                                             "O1", "Oz", "O2")
beta_norm_alphb<-bands_norm$beta%>% select("Time", "Fp1", "Fp2",
                                           "AF3", "AF4",
                                           "F7", "F3", "Fz", "F4", "F8",
                                           "FC5", "FC1", "FC2", "FC6",
                                           "T7", "C3", "Cz", "C4", "T8",
                                           "CP5", "CP1", "CP2", "CP6",
                                           "P7", "P3", "Pz", "P4", "P8",
                                           "PO3", "PO4",
                                           "O1", "Oz", "O2")
theta_norm_alphb<-bands_norm$theta%>% select("Time", "Fp1", "Fp2",
                                             "AF3", "AF4",
                                             "F7", "F3", "Fz", "F4", "F8",
                                             "FC5", "FC1", "FC2", "FC6",
                                             "T7", "C3", "Cz", "C4", "T8",
                                             "CP5", "CP1", "CP2", "CP6",
                                             "P7", "P3", "Pz", "P4", "P8",
                                             "PO3", "PO4",
                                             "O1", "Oz", "O2")
delta_norm_alphb<-bands_norm$delta%>% select("Time", "Fp1", "Fp2",
                                             "AF3", "AF4",
                                             "F7", "F3", "Fz", "F4", "F8",
                                             "FC5", "FC1", "FC2", "FC6",
                                             "T7", "C3", "Cz", "C4", "T8",
                                             "CP5", "CP1", "CP2", "CP6",
                                             "P7", "P3", "Pz", "P4", "P8",
                                             "PO3", "PO4",
                                             "O1", "Oz", "O2")
gamma_norm_alphb<-bands_norm$gamma%>% select("Time", "Fp1", "Fp2",
                                             "AF3", "AF4",
                                             "F7", "F3", "Fz", "F4", "F8",
                                             "FC5", "FC1", "FC2", "FC6",
                                             "T7", "C3", "Cz", "C4", "T8",
                                             "CP5", "CP1", "CP2", "CP6",
                                             "P7", "P3", "Pz", "P4", "P8",
                                             "PO3", "PO4",
                                             "O1", "Oz", "O2")
norm_alphb <- list(alpha_norm_alphb, beta_norm_alphb, delta_norm_alphb, theta_norm_alphb, gamma_norm_alphb)
names(norm_alphb) <- c("alpha", "beta", "delta", "theta", "gamma")

dtw_by_band_ordered
# Basic boxplot
mah <- dtw_by_band_ordered$theta ##edit for specific wavelength
boxplot(mah,
        main = "Boxplot of Theta Wavelength Data",
        ylab = "Value (μV or similar unit)",
        col = "lightblue",
        border = "darkblue")

library(openxlsx)

wb <- createWorkbook()

for (band in names(dtw_by_band_ordered)) {
  addWorksheet(wb, band)
  writeData(
    wb,
    sheet = band,
    x = dtw_by_band_ordered[[band]],
    rowNames = TRUE
  )
}

saveWorkbook(wb, "DTW_32x32_all_bands_ordered.xlsx", overwrite = TRUE)

#Rafael modeling messing

plot.ts(norm_alphb$alpha$AF3, ty='p')

y<-norm_alphb$alpha$AF3
y <- bands_norm$alpha$Fp1

x <- 1:length(y)

fit_linear <- lm(y ~ x)
fit_linear
abline(fit_linear, lwd = 2, col = 2)

library(splines)

fit_splines <- lm(y ~ bs(x, 100))
y_hat <- predict(fit_splines)

lines(x, y_hat, col = 4, lwd = 3)

library(mgcv)
fit_mgcv <- gam(y ~ s(x, k = 100))
summary(fit_mgcv)

library(gratia)
draw(fit_mgcv)

# perhaps pre-process with a moving average (need to figure out optimal window)

y_ma100 <- stats::filter(y, rep(1/100, 100), sides = 1) #https://statisticsglobe.com/moving-average-maximum-median-sum-of-time-series-in-r
plot.ts(y)
lines(y_ma100, col = 2)
plot(y_ma100)

fit_mgcv <- gam(y_ma100 ~ s(x, k = 100))
summary(fit_mgcv)

lines(predict(fit_mgcv), col = 4, lwd = 2)

library(gratia)
draw(fit_mgcv)

