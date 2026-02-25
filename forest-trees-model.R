y <- ts_data[1:1000]
plot(y)

library(randomForest)

fit <- randomForest(x = poly(1:length(y), 10), y = y, ntree = 10)
plot(fit)

plot(y, cex = .5)
lines(predict(fit), col = 2, lwd = 2)
