lapply(results$per_band, function(x) {
  y <- apply(x, 2, function(x) which(any(is.na(x))))
  return(y)
  }) %>% unlist

results$per_band$alpha

plot.ts(alpha$P7)

norm_x <- log10(alpha$P7)
time <- 1:length(norm_x)
k <- 40
fit <- gam(norm_x ~ s(time, k = k))
smooth_x    <- as.numeric(predict(fit, data.frame(time = time)))  

plot(time, norm_x, cex = .2)
lines(time, smooth_x, ty = 'l', col = 2, lwd = 4)
abline(h = mean(smooth_x), lwd = 4, col = 4)
abline(h = mean(smooth_x)+sd(smooth_x), lwd = 4, col = 4, lty = 2)

# Detect peaks on smooth signal
peaks <- findpeaks(smooth_x,
                   minpeakheight   = mean(smooth_x) + sd(smooth_x),
                   minpeakdistance = 10)

abline(v = peaks[,2])
abline(v = peaks[,4],col=2)
