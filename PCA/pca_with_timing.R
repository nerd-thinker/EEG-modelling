library(tidyverse)

variances <- apply(features_clean[, -1], 2, var)
features_for_pca <- features_clean[, c(TRUE, variances > 0)]

cat("  Features with zero variance:", sum(variances == 0), "\n")
cat("  Features to keep:", sum(variances > 0), "\n")
cat("  Data dimensions:", nrow(features_for_pca), "participants ×", 
    ncol(features_for_pca) - 1, "features\n")

pca_all <- prcomp(features_for_pca[, -1], 
                  center = TRUE, 
                  scale. = TRUE)

# Variance explained
var_explained <- cumsum(pca_all$sdev^2 / sum(pca_all$sdev^2))

# How many PCs for 80%, 90%, 95%?
n_pc_80 <- which(var_explained >= 0.80)[1]
n_pc_90 <- which(var_explained >= 0.90)[1]
n_pc_95 <- which(var_explained >= 0.95)[1]

plot(pca_all, type = "l")

library(nFactors) ##doesn't work with high dimension of PCA's
library(paran)
paran_result <- paran(features_for_pca[, -1],
                      iterations = 500,
                      quietly = FALSE,
                      seed = 123)
