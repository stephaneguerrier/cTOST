# Build pkgdown site with README figures
# Run this script to build your pkgdown website with proper asciicast output

# Clean build to avoid memory issues
if (dir.exists("docs")) {
  message("Removing old docs directory...")
  unlink("docs", recursive = TRUE)
}

library(pkgdown)

# Build the site with clean cache
message("Building pkgdown site...")
build_site(override = list(destination = "docs"), devel = FALSE, preview = FALSE)

# Copy README figures to docs
if (dir.exists("README_files")) {
  message("Copying README_files to docs/...")

  # Remove old README_files if it exists
  if (dir.exists("docs/README_files")) {
    unlink("docs/README_files", recursive = TRUE)
  }

  # Copy new README_files
  file.copy("README_files", "docs/", recursive = TRUE)

  message("✓ README figures copied successfully!")
} else {
  warning("README_files directory not found. Make sure to knit README.Rmd first.")
}

message("✓ Site build complete!")
