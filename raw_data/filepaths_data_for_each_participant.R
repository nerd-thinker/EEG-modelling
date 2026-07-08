# filepaths_data_for_each_participant.R
# -----------------------------------------------------------------------------
# Defines file paths for all 49 participants in the SPUR-EEG dataset.
# All files live in a single flat folder — no subdirectories.
# Generated programmatically from participant-data.txt; do not edit by hand.
# -----------------------------------------------------------------------------

base_dir <- "~/eeg/raw_data/SPUR-EEG-data"

# ---- Participant IDs --------------------------------------------------------
# 49 participants. Note: CTEEG023Y and CTEEG026Y are absent from the dataset.

participant_ids <- c(
  "CTEEG001O", "CTEEG001Y",
  "CTEEG002O", "CTEEG002Y",
  "CTEEG003O",
  "CTEEG004O", "CTEEG004Y",
  "CTEEG005O", "CTEEG005Y",
  "CTEEG006O", "CTEEG006Y",
  "CTEEG007O",
  "CTEEG008O", "CTEEG008Y",
  "CTEEG009O", "CTEEG009Y",
  "CTEEG010O", "CTEEG010Y",
  "CTEEG011O",
  "CTEEG012O", "CTEEG012Y",
  "CTEEG013O", "CTEEG013Y",
  "CTEEG014O", "CTEEG014Y",
  "CTEEG015O",
  "CTEEG016O",
  "CTEEG017O", "CTEEG017Y",
  "CTEEG018O", "CTEEG018Y",
  "CTEEG019O",
  "CTEEG020O",
  "CTEEG021Y", "CTEEG022Y",
  "CTEEG024Y", "CTEEG025Y",
  "CTEEG027Y", "CTEEG028Y", "CTEEG029Y", "CTEEG030Y",
  "CTEEG031Y", "CTEEG032Y", "CTEEG033Y", "CTEEG034Y",
  "CTEEG035Y", "CTEEG036Y", "CTEEG037Y", "CTEEG038Y"
)

# ---- Build named filepath vectors for every participant --------------------
# Each entry: named vector with bands as names and full file paths as values.
# Band name "gamma" uses the "Gamma1" filename suffix (as in the raw data).

make_filepaths <- function(pid) {
  c(
    alpha = file.path(base_dir, paste0(pid, "_Alpha.csv")),
    beta  = file.path(base_dir, paste0(pid, "_Beta.csv")),
    delta = file.path(base_dir, paste0(pid, "_Delta.csv")),
    gamma = file.path(base_dir, paste0(pid, "_Gamma1.csv")),
    theta = file.path(base_dir, paste0(pid, "_Theta.csv"))
  )
}

all_filepaths <- setNames(
  lapply(participant_ids, make_filepaths),
  participant_ids
)

# ---- Quick sanity check (runs on source) -----------------------------------
cat(sprintf("Loaded filepaths for %d participants, %d bands each (%d files total)\n",
            length(all_filepaths),
            length(all_filepaths[[1]]),
            length(all_filepaths) * length(all_filepaths[[1]])))
