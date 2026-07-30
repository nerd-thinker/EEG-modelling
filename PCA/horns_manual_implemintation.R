## Instead of taking specific number of PCs (25/34) implement Horn's parallel analysis to see which are 
## actually important. Horn's parallel analysis compares the eigenvalues of the actual data with those 
## from random data to determine the number of significant components.

# Step 1: extract actual eigenvalues
actual_eigenvalues <- pca_all$sdev^2
cat("Num of eigenvalues:", length(actual_eigenvalues), "\n")

# Step 2: Generate random data and compute eigenvalues
## Source pca_with_timing.R first

n_iterations <- 1000 
n_participants <- nrow(features_for_pca) -1 # Exclude participant ID
n_features <- ncol(features_for_pca) - 1  # Exclude participant ID

random_eigenvalues_all <- matrix(0, nrow = n_iterations, ncol = n_participants) #storage for random eigenvalues

set.seed(123)
for( i in 1:n_iterations){
  random_data <- matrix(rnorm(n_participants * n_features), 
                        nrow = n_participants, 
                        ncol = n_features)
  random_pca <- prcomp(random_data, center = TRUE, scale. = TRUE)
  random_eigenvalues_all[i, ] <- random_pca$sdev^2
}

# Step 3: Compute mean and 95th percentile of random eigenvalues
mean_random_eigenvalues <- colMeans(random_eigenvalues_all)
percentile_95 <- apply(random_eigenvalues_all, 2, quantile, probs = 0.95)

# Step 4: Compare actual eigenvalues with random eigenvalues
n_significant_components <- sum(actual_eigenvalues > percentile_95)

# Results
results_table <- data.frame(
  PC = 1:min(length(actual_eigenvalues), length(percentile_95)),
  Actual = actual_eigenvalues[1:min(length(actual_eigenvalues), length(percentile_95))],
  Mean_Random = mean_random_eigenvalues[1:min(length(actual_eigenvalues), length(percentile_95))],
  Percentile_95 = percentile_95[1:min(length(actual_eigenvalues), length(percentile_95))],
  Significant = actual_eigenvalues[1:min(length(actual_eigenvalues), length(percentile_95))] > 
    percentile_95[1:min(length(actual_eigenvalues), length(percentile_95))]
)
print(results_table) # PC 48 TRUE but artifact due to rounding error (false positive)

# Step 5: Compare to variance-based threshold
var_explained <- cumsum(actual_eigenvalues / sum(actual_eigenvalues))
n_pc_80 <- which(var_explained >= 0.80)[1]
n_pc_90 <- which(var_explained >= 0.90)[1]

# Step 6: extract and save PC scores
pc_scores <- pca_all$x[, 1:15]  # Extract first 15 PCs that are significant
pc_data <- cbind(
  participant_id = features_for_pca$participant_id, 
  as.data.frame(pc_scores)
)

