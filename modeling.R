library('forecast')  # moving average
library(splines)     # regression spline function
library(mgcv)        # Mixed GAM Computation Vehicle with Automatic Smoothness Estimation
library(gratia)      # different GAM library 
library(zoo)         # Infrastructure for Regular and Irregular Time Series
library(ggplot2)     
library(tidyverse)


# Running across only Fp1 node of alpha wavelength of normalized E --------

# Moving average Model

ts_data <- bands_norm$alpha$Fp1

sma_result <- TTR::SMA(ts_data, n = 150) ##how to find perfect window?

plot(ts_data, col = 3, main = "Time series with moving average", t = "lines") 
lines(sma_result, col = 2)

## exponential makes no sense
ema_result <- stats::filter(ts_data, filter = .2, method = "recursive")
lines(ema_result, col = 4)

## moving average different approach (ai)
df_alpha <- bands_norm$alpha

df_alpha <- df_alpha %>% 
  mutate(Fp1_smoothed = rollmean(Fp1, k = 200, fil = NA))

ggplot(df_alpha, aes(x = Time)) +
  geom_line(aes(y = Fp1), alpha = 0.3) + # Raw data in background
  geom_line(aes(y = Fp1_smoothed), color = "red") + # Smoothed line
  theme_minimal()

# Generalized Additive model
# k = 20 is the "basis dimension" - higher means more "wiggly"
gam_model <- gam(Fp1 ~ s(Time, k = 20), data = df_alpha)

summary(gam_model)

df_alpha$gam_preds <- predict(gam_model, df_alpha)

## plot
ggplot(df_alpha, aes(x = Time, y = Fp1)) +
  geom_point(alpha = 0.1) +
  geom_line(aes(y = gam_preds), color = "blue", size = 1) +
  labs(title = "GAM Fit for Alpha Band (Fp1)")


# Running GAM on smoothed eeg across all nodes and wavelength -------------
library(tidyverse)
library(mgcv)

# Replace 'bands_norm$alpha' with your actual data source
## pulling all channels together into 3 columns: time, Channel (32 pooled together) 
## and power (normalized) acrooss each channell?? across wavelength
alpha_long <- bands_norm$alpha %>% 
  pivot_longer(
    cols = -Time, 
    names_to = "Channel", 
    values_to = "Power"
  )
##not normalized bands -- raw
alpha_long_raw <- bands_clean$alpha %>%
  pivot_longer(
    cols = -Time,
    names_to = "Channel",
    values_to = "Power"
  )

# GAM computation for all channels across wavelength
alpha_results <- alpha_long_raw %>%
  group_by(Channel) %>%
  nest() %>%
  mutate(
    model = map(data, ~ gam(Power ~ s(Time, k = 20), data = .x)),
    # Extract the smoothed fit
    smoothed_power = map2(model, data, ~ predict(.x, .y))
  ) %>%
  unnest(cols = c(data, smoothed_power))

# Visualizing Frontal vs Occipital channels
alpha_results %>%
  filter(Channel %in% c("Fp1", "Fp2", "O1", "O2")) %>%
  ggplot(aes(x = Time)) +
  geom_line(aes(y = Power), alpha = 0.2) + # Raw noise
  geom_line(aes(y = smoothed_power), color = "blue", size = 1) + # GAM Trend
  facet_wrap(~Channel) +
  theme_minimal() +
  labs(title = "Alpha Band: Raw Data vs. GAM Smoothing")

## beta GAM computation
beta_long <- bands_norm$beta %>% 
  pivot_longer(
    cols = -Time, 
    names_to = "Channel", 
    values_to = "Power"
  )

beta_results <- beta_long %>%
  group_by(Channel) %>%
  nest() %>%
  mutate(
    model = map(data, ~ gam(Power ~ s(Time, k = 20), data = .x)),
    # Extract the smoothed fit
    smoothed_power = map2(model, data, ~ predict(.x, .y))
  ) %>%
  unnest(cols = c(data, smoothed_power))

## delta GAM computation
delta_long <- bands_norm$delta %>% 
  pivot_longer(
    cols = -Time, 
    names_to = "Channel", 
    values_to = "Power"
  )

delta_results <- delta_long %>%
  group_by(Channel) %>%
  nest() %>%
  mutate(
    model = map(data, ~ gam(Power ~ s(Time, k = 20), data = .x)),
    # Extract the smoothed fit
    smoothed_power = map2(model, data, ~ predict(.x, .y))
  ) %>%
  unnest(cols = c(data, smoothed_power))

## gamma GAM computation
gamma_long <- bands_norm$gamma %>% 
  pivot_longer(
    cols = -Time, 
    names_to = "Channel", 
    values_to = "Power"
  )

gamma_results <- gamma_long %>%
  group_by(Channel) %>%
  nest() %>%
  mutate(
    model = map(data, ~ gam(Power ~ s(Time, k = 20), data = .x)),
    # Extract the smoothed fit
    smoothed_power = map2(model, data, ~ predict(.x, .y))
  ) %>%
  unnest(cols = c(data, smoothed_power))

## theta GAM computation
theta_long <- bands_norm$theta %>% 
  pivot_longer(
    cols = -Time, 
    names_to = "Channel", 
    values_to = "Power"
  )

theta_results <- theta_long %>%
  group_by(Channel) %>%
  nest() %>%
  mutate(
    model = map(data, ~ gam(Power ~ s(Time, k = 20), data = .x)),
    # Extract the smoothed fit
    smoothed_power = map2(model, data, ~ predict(.x, .y))
  ) %>%
  unnest(cols = c(data, smoothed_power))
