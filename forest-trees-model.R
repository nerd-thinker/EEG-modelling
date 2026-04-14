y <- ts_data[1:1000]
plot(y)

library(randomForest)

fit <- randomForest(x = poly(1:length(y), 10), y = y, ntree = 10)
plot(fit)

plot(y, cex = .5)
lines(predict(fit), col = 2, lwd = 2)

library(ggplot2)
library(tidyverse)

ggplot(bands_clean$alpha, aes(x = Time, y = Fp1)) + geom_line() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Alpha Frequency Fp1 electrode Time Series", x = "Time in seconds", y = "Amplitude")

ggplot(bands_clean$gamma, aes(x = Time, y = Fp1 )) + geom_line() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Gamma Frequency Fp1 electrode Time Series", x = "Time in seconds", y = "Amplitude")

# For alpha band - 5 specific electrodes on one plot
bands_clean$alpha %>%
  pivot_longer(
    cols = -Time,
    names_to = "Channel",
    values_to = "Power"
  ) %>%
  filter(Channel %in% c("Fp1", "Fp2", "O1", "O2", "Cz")) %>%
  ggplot(aes(x = Time, y = Power, color = Channel)) +
  geom_line(size = 1) +
  labs(
    title = "Alpha Band Raw Activity Across 5 Electrodes",
    x = "Time (seconds)",
    y = "Power"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14)
  )
bands_clean$alpha %>%
  pivot_longer(
    cols = -Time,
    names_to = "Channel",
    values_to = "Power"
  ) %>%
  filter(Channel %in% c("Fp1", "Fp2", "O1", "O2", "Cz")) %>%
  ggplot(aes(x = Time, y = Power, color = Channel)) +
  geom_line(size = 0.8) +
  facet_wrap(~Channel, ncol = 1) +
  labs(
    title = "Alpha Band Raw Activity Across 5 Electrodes",
    x = "Time (seconds)",
    y = "Power"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold")
  )
