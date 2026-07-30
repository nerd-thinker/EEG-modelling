library(tidyverse)

variances <- apply(features_clean[, -1], 2, var)
features_for_pca <- features_clean[, c(TRUE, variances > 0)]

cat("  Features with zero variance:", sum(variances == 0), "\n")
cat("  Features to keep:", sum(variances > 0), "\n")
cat("  Data dimensions:", nrow(features_for_pca), "participants ×", 
    ncol(features_for_pca) - 1, "features\n")
##PCA
pca_results <- prcomp(features_for_pca[, -1], 
                  center = TRUE, 
                  scale. = TRUE)

pca_all <- readRDS("pca_results.rds")  # Load PCA results 
# Variance explained
var_explained <- cumsum(pca_all$sdev^2 / sum(pca_all$sdev^2))

# How many PCs for 80%, 90%, 95%?
n_pc_80 <- which(var_explained >= 0.80)[1]
n_pc_90 <- which(var_explained >= 0.90)[1]
n_pc_95 <- which(var_explained >= 0.95)[1]

# Scree plot (how variance drops off)
plot(pca_all, type = "l")

# Cumulative variance
plot(cumsum(pca_all$sdev^2 / sum(pca_all$sdev^2)),
     xlab = "PC", ylab = "Cumulative Variance Explained",
     type = "b")
abline(h = 0.8, lty = 2, col = "red")  # Mark 80%
abline(v = 25, lty = 2, col = "blue")  # Mark PC25


