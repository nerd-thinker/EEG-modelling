library(readxl)



filepaths <- c(
  alpha = "/home/ks/EGG-modeling/raw_data/CTEEG022Y_Alpha.xlsx",
  beta  = "/home/ks/EGG-modeling/raw_data/CTEEG022Y_Beta.xlsx",
  gamma = "/home/ks/EGG-modeling/raw_data/CTEEG022Y_Gamma1.xlsx", 
  delta = "/home/ks/EGG-modeling/raw_data/CTEEG022Y_Delta.xlsx",
  theta = "/home/ks/EGG-modeling/raw_data/CTEEG022Y_Theta.xlsx"
)

bands <- names(filepaths)  # band names come from the named vector

# Step 1: Read all files
message("Reading Excel files...")
eeg_data1 <- lapply(setNames(filepaths, bands), function(fp) {
  read_excel(fp)
})

# Extract time and channels from first file
time_vec <- eeg_data1[[1]]$Time
channels <- setdiff(colnames(eeg_data1[[1]]), "Time")
saveRDS(eeg_data1, file = "eeg_data", compress = F)
