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
