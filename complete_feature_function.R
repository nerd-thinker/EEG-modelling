# Helper functions -------------
hjorth <- function(x) {
  x <- na.omit(x)
  
  dx  <- diff(x)        # first derivative
  ddx <- diff(dx)       # second derivative
  
  activity   <- var(x)
  mobility   <- sqrt(var(dx) / var(x)) ## estimate of the mean frequency
  complexity <- (sqrt(var(ddx) / var(dx))) / mobility ## estimate of the bandwidth of the signal, which indicates the similarity of the shape of the signal to a pure sine wave
  
  list(activity = activity, mobility = mobility, complexity = complexity)
}

skewness <- function(x) mean(((x - mean(x)) / sd(x))^3)
kurtosis <- function(x) mean(((x - mean(x)) / sd(x))^4)
band_power <- function(x) sum(x^2) / length(x)

spectral_entropy <- function(x) {
  spec <- abs(fft(x))^2
  spec <- spec / sum(spec)
  spec <- spec[spec > 0]
  -sum(spec * log(spec))
}

