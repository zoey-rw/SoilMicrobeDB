# Test automatic download when file is missing locally
# This simulates what happens when a new user runs a script

# Load configuration
source("scripts/create_database/config.R")

cat("=== Testing Automatic Download from GitHub ===\n\n")

# Temporarily remove the GTDB mapping file to test auto-download
backup_file <- file.path(GTDB_MAPPING_LOCAL_DIR, "GTDB_NCBI_key.rds.backup")
local_file <- file.path(GTDB_MAPPING_LOCAL_DIR, "GTDB_NCBI_key.rds")

if (file.exists(backup_file) && !file.exists(local_file)) {
  cat("Test: GTDB mapping file is missing locally\n")
  cat("Expected: Should download from GitHub automatically\n\n")
  
  # The config should try to download it
  cat("Checking GTDB_NCBI_MAPPING_FILE path...\n")
  cat("  Path:", GTDB_NCBI_MAPPING_FILE, "\n")
  cat("  Exists:", file.exists(GTDB_NCBI_MAPPING_FILE), "\n\n")
  
  # Manually trigger download to test
  if (!file.exists(local_file)) {
    github_url <- paste0(GITHUB_REPO_BASE, "/gtdb_mapping/GTDB_NCBI_key.rds")
    cat("Attempting to download from GitHub...\n")
    cat("  URL:", github_url, "\n")
    
    tryCatch({
      result <- download_from_github(github_url, local_file, "GTDB-NCBI mapping file")
      if (file.exists(local_file) && file.size(local_file) > 0) {
        cat("✓ Successfully downloaded GTDB mapping file (", file.size(local_file), " bytes)\n")
        cat("  File is now available at:", local_file, "\n\n")
      }
    }, error = function(e) {
      cat("✗ Download failed:", e$message, "\n\n")
    })
  }
} else {
  cat("Note: GTDB mapping file already exists locally\n")
  cat("  To test download, remove:", local_file, "\n")
  cat("  Backup exists at:", backup_file, "\n\n")
}

cat("=== Test Complete ===\n")

