base_dir <- "~/eeg/raw_data/cleaned-data/SPUR-EEG-data"

# Participant: 001Y  (no Alpha file)
# time vector present, nodes present
filepaths_001Y <- c(
  beta  = file.path(base_dir, "001Y_Beta.csv"),
  delta = file.path(base_dir, "001Y_Delta.csv"),
  gamma = file.path(base_dir, "001Y_Gamma1.csv"),
  theta = file.path(base_dir, "001Y_Theta.csv")
)

# Participant: 002Y  (no Beta file)
# no time vector no nodes, 32 rows
filepaths_002Y <- c(
  alpha = file.path(base_dir, "002Y_Alpha.csv"),
  delta = file.path(base_dir, "002Y_Delta.csv"),
  gamma = file.path(base_dir, "002Y_Gamma.csv"),
  theta = file.path(base_dir, "002Y_Theta.csv")
)

# Participant: 004O  (no Delta file)
# no time vector no nodes, 32 rows
filepaths_004O <- c(
  alpha = file.path(base_dir, "004O_Alpha.csv"),
  beta  = file.path(base_dir, "004O_Beta.csv"),
  gamma = file.path(base_dir, "004O_Gamma.csv"),
  theta = file.path(base_dir, "004O_Theta.csv")
)

# Participant: 008O  (no Gamma file)
# no time vector no nodes, 32 rows
filepaths_008O <- c(
  alpha = file.path(base_dir, "008O_Alpha.csv"),
  beta  = file.path(base_dir, "008O_Beta.csv"),
  delta = file.path(base_dir, "008O_Delta.csv"),
  theta = file.path(base_dir, "008O_Theta.csv")
)

# Participant: 008Y (Theta missing)
# no time vector no nodes, 32 rows
filepaths_008Y <- c(
  alpha = file.path(base_dir, "008Y_Alpha.csv"),
  beta  = file.path(base_dir, "008Y_Beta.csv"),
  delta = file.path(base_dir, "008Y_Delta.csv"),
  gamma = file.path(base_dir, "008Y_Gamma.csv")
)


# Participant: 012O 
# no time vector no nodes, 32 rows
filepaths_012O <- c(
  alpha = file.path(base_dir, "012O_Alpha.csv"),
  beta  = file.path(base_dir, "012O_Beta.csv"),
  delta = file.path(base_dir, "012O_Delta.csv"),
  gamma = file.path(base_dir, "012O_Gamma.csv"),
  theta = file.path(base_dir, "012O_Theta.csv")
)

# Participant: 013Y  (no Alpha file)
# no time vector no nodes, 32 rows
filepaths_013Y <- c(
  beta  = file.path(base_dir, "013Y_Beta.csv"),
  delta = file.path(base_dir, "013Y_Delta.csv"),
  gamma = file.path(base_dir, "013Y_Gamma.csv"),
  theta = file.path(base_dir, "013Y_Theta.csv")
)

# Participant: 016O  (no Beta file)
# no time vector no nodes, 32 rows
filepaths_016O <- c(
  alpha = file.path(base_dir, "016O_Alpha.csv"),
  delta = file.path(base_dir, "016O_Delta.csv"),
  gamma = file.path(base_dir, "016O_Gamma.csv"),
  theta = file.path(base_dir, "016O_Theta.csv")
)

# Participant: 017Y  (no Delta file)
# no time vector no nodes, 32 rows
filepaths_017Y <- c(
  alpha = file.path(base_dir, "017Y_Alpha.csv"),
  beta  = file.path(base_dir, "017Y_Beta.csv"),
  gamma = file.path(base_dir, "017Y_Gamma.csv"),
  theta = file.path(base_dir, "017Y_Theta.csv")
)

# Participant: 024Y  (no Gamma file)
# no time vector no nodes, 32 rows
filepaths_024Y <- c(
  alpha = file.path(base_dir, "024Y_Alpha.csv"),
  beta  = file.path(base_dir, "024Y_Beta.csv"),
  delta = file.path(base_dir, "024Y_Delta.csv"),
  theta = file.path(base_dir, "024Y_Theta.csv")
)

# Participant: 029Y (no Theta file)
# no time vector no nodes, 32 rows
filepaths_029Y <- c(
  alpha = file.path(base_dir, "029Y_Alpha.csv"),
  beta  = file.path(base_dir, "029Y_Beta.csv"),
  delta = file.path(base_dir, "029Y_Delta.csv"),
  gamma = file.path(base_dir, "029Y_Gamma.csv")
)

# Participant: 033Y 
# no time vector no nodes, 32 rows
filepaths_033Y <- c(
  alpha = file.path(base_dir, "033Y_Alpha.csv"),
  beta  = file.path(base_dir, "033Y_Beta.csv"),
  delta = file.path(base_dir, "033Y_Delta.csv"),
  gamma = file.path(base_dir, "033Y_Gamma.csv"),
  theta = file.path(base_dir, "033Y_Theta.csv")
)

# Participant: 038Y  (no Alpha file)
# no time vector no nodes, 32 rows
filepaths_038Y <- c(
  beta  = file.path(base_dir, "038Y_Beta.csv"),
  delta = file.path(base_dir, "038Y_Delta.csv"),
  gamma = file.path(base_dir, "038Y_Gamma.csv"),
  theta = file.path(base_dir, "038Y_Theta.csv")
)

# Participant: CTEEG003O
# time vector present, nodes present
filepaths_CTEEG003O <- c(
  alpha = file.path(base_dir, "CTEEG003O_Alpha.csv"),
  beta  = file.path(base_dir, "CTEEG003O_Beta.csv"),
  delta = file.path(base_dir, "CTEEG003O_Delta.csv"),
  gamma = file.path(base_dir, "CTEEG003O_Gamma1.csv"),
  theta = file.path(base_dir, "CTEEG003O_Theta.csv")
)

# Participant: CTEEG006Y
# time vector present, nodes present
filepaths_CTEEG006Y <- c(
  alpha = file.path(base_dir, "CTEEG006Y_Alpha.csv"),
  beta  = file.path(base_dir, "CTEEG006Y_Beta.csv"),
  delta = file.path(base_dir, "CTEEG006Y_Delta.csv"),
  gamma = file.path(base_dir, "CTEEG006Y_Gamma1.csv"),
  theta = file.path(base_dir, "CTEEG006Y_Theta.csv")
)

# Participant: CTEEG012Y
# time vector present, nodes present
filepaths_CTEEG012Y <- c(
  alpha = file.path(base_dir, "CTEEG012Y_Alpha.csv"),
  beta  = file.path(base_dir, "CTEEG012Y_Beta.csv"),
  delta = file.path(base_dir, "CTEEG012Y_Delta.csv"),
  gamma = file.path(base_dir, "CTEEG012Y_Gamma1.csv"),
  theta = file.path(base_dir, "CTEEG012Y_Theta.csv")
)

# Participant: CTEEG014Y
# time vector present, nodes present
filepaths_CTEEG014Y <- c(
  alpha = file.path(base_dir, "CTEEG014Y_Alpha.csv"),
  beta  = file.path(base_dir, "CTEEG014Y_Beta.csv"),
  delta = file.path(base_dir, "CTEEG014Y_Delta.csv"),
  gamma = file.path(base_dir, "CTEEG014Y_Gamma1.csv"),
  theta = file.path(base_dir, "CTEEG014Y_Theta.csv")
)

# Participant: CTEEG022Y (ALREADY HAVE DATA ON THAT)
# time vector present, nodes present
filepaths_CTEEG022Y <- c(
  alpha = file.path(base_dir, "CTEEG022Y_Alpha.csv"),
  beta  = file.path(base_dir, "CTEEG022Y_Beta.csv"),
  delta = file.path(base_dir, "CTEEG022Y_Delta.csv"),
  gamma = file.path(base_dir, "CTEEG022Y_Gamma1.csv"),
  theta = file.path(base_dir, "CTEEG022Y_Theta.csv")
)

# Participant: CTEEG037Y
# time vector present, nodes present
filepaths_CTEEG037Y <- c(
  alpha = file.path(base_dir, "CTEEG037Y_Alpha.csv"),
  beta  = file.path(base_dir, "CTEEG037Y_Beta.csv"),
  delta = file.path(base_dir, "CTEEG037Y_Delta.csv"),
  gamma = file.path(base_dir, "CTEEG037Y_Gamma1.csv"),
  theta = file.path(base_dir, "CTEEG037Y_Theta.csv")
)

# All participants as a named list — use to loop over all with extract_all_features()
all_filepaths <- list(
  "001Y"      = filepaths_001Y,
  "002Y"      = filepaths_002Y,
  "004O"      = filepaths_004O,
  "008O"      = filepaths_008O,
  "008Y"      = filepaths_008Y,
  "012O"      = filepaths_012O,
  "013Y"      = filepaths_013Y,
  "016O"      = filepaths_016O,
  "017Y"      = filepaths_017Y,
  "024Y"      = filepaths_024Y,
  "029Y"      = filepaths_029Y,
  "033Y"      = filepaths_033Y,
  "038Y"      = filepaths_038Y,
  "CTEEG003O" = filepaths_CTEEG003O,
  "CTEEG006Y" = filepaths_CTEEG006Y,
  "CTEEG012Y" = filepaths_CTEEG012Y,
  "CTEEG014Y" = filepaths_CTEEG014Y,
  "CTEEG022Y" = filepaths_CTEEG022Y,
  "CTEEG037Y" = filepaths_CTEEG037Y
)

# Example: run extract_all_features() across all participants
# results_all <- lapply(all_filepaths, function(fp) {
#   extract_all_features(filepaths = fp, scaling = "z_score", k = 20)
# })

# run for the participant ---------------
participant_row <- extract_all_features(
  filepaths      = filepaths_CTEEG003O,
  participant_id = "CTEEG003O",
  k              = 70
)
