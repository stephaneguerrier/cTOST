# Define path
folder <- "ctost/proportion"
path <- "ctost/proportion/data_temp"

# List files
all_files <- list.files(path = path, pattern = "\\.rda$")

if (length(all_files) == 0) {
  stop("No .rda files found in ", path)
}

# Load first file to get dimensions
load(paste0(path, "/", all_files[1]))
ncol_file <- ncol(res)

# Create matrix to save all results
df_all_results <- matrix(NA, ncol = ncol_file)

# For all files, load and bind
for (file_index in seq_along(all_files)) {
  file_i <- all_files[file_index]
  file_name <- paste0(path, "/", file_i)
  load(file_name)
  df_all_results <- rbind(df_all_results, res)

  # Print progress every 1000 files
  if (file_index %% 1000 == 0) {
    cat(sprintf("Processed %d / %d files\n", file_index, length(all_files)))
  }
}

# Remove first row (NA placeholder)
df_all_results <- df_all_results[-1, ]

# Add column names
colnames(df_all_results) <- c(
  "decision_dp", "ci_lower_dp", "ci_upper_dp",
  "decision_classical", "ci_lower_classical", "ci_upper_classical",
  "p1_hat", "p2_hat", "p1_obs", "p2_obs",
  "p1_true", "p2_true", "n", "epsilon", "scenario", "id_slurm", "seed"
)

# Print summary
cat(sprintf("\nTotal simulations: %d\n", nrow(df_all_results)))
cat(sprintf("Total files processed: %d\n", length(all_files)))

# Save matrix of results with timestamp
time <- Sys.time()
time_2 <- gsub(" ", "_", time)
time_3 <- gsub(":", "-", time_2)
file_name_to_save <- paste0(folder, "/simu_dp_proportion_", time_3, ".rda")
cat(sprintf("Saving to: %s\n", file_name_to_save))
save(df_all_results, file = file_name_to_save)
