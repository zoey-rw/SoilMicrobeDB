# Test script to verify GitHub download functionality
# This script tests that reference data files can be downloaded from GitHub

# Load configuration
source("scripts/create_database/config.R")

cat("=== Testing GitHub Download Functionality ===\n\n")

# Test 1: Check if download_from_github function exists
cat("Test 1: Checking if download_from_github function exists...\n")
if (exists("download_from_github")) {
  cat("✓ download_from_github function found\n\n")
} else {
  cat("✗ download_from_github function not found\n\n")
  stop("Helper functions not loaded correctly")
}

# Test 2: Test downloading a small file (RefSoil)
cat("Test 2: Testing download of RefSoil data file from GitHub...\n")
test_dir <- file.path(GENOME_DB_DIR, "test_downloads")
dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)
test_file <- file.path(test_dir, "RefSoil_s1_test.xlsx")

# Remove test file if it exists
if (file.exists(test_file)) {
  unlink(test_file)
}

github_url <- paste0(GITHUB_REPO_BASE, "/RefSoil/ref_data/RefSoil_s1.xlsx")
cat("  GitHub URL:", github_url, "\n")
cat("  Local path:", test_file, "\n")

tryCatch({
  result <- download_from_github(github_url, test_file, "RefSoil test file")
  if (file.exists(test_file) && file.size(test_file) > 0) {
    cat("✓ Successfully downloaded RefSoil file (", file.size(test_file), " bytes)\n\n")
    # Clean up
    unlink(test_file)
  } else {
    cat("✗ Download failed - file not found or empty\n\n")
  }
}, error = function(e) {
  cat("✗ Download failed with error:", e$message, "\n\n")
})

# Test 3: Test that config paths are set correctly
cat("Test 3: Checking configuration paths...\n")
cat("  GITHUB_REPO_BASE:", GITHUB_REPO_BASE, "\n")
cat("  GTDB_NCBI_MAPPING_FILE:", GTDB_NCBI_MAPPING_FILE, "\n")
cat("  REFSOIL_CONFIG$data_file:", REFSOIL_CONFIG$data_file, "\n")
cat("  SPIRE_CONFIG$taxonomy_file:", SPIRE_CONFIG$taxonomy_file, "\n")
cat("  GEM_CONFIG$supplementary_file:", GEM_CONFIG$supplementary_file, "\n\n")

# Test 4: Verify files exist locally or can be accessed
cat("Test 4: Checking if reference files exist locally...\n")
files_to_check <- list(
  "GTDB mapping" = GTDB_NCBI_MAPPING_FILE,
  "RefSoil data" = REFSOIL_CONFIG$data_file,
  "SPIRE taxonomy" = SPIRE_CONFIG$taxonomy_file,
  "SPIRE microntology" = SPIRE_CONFIG$microntology_file,
  "GEM supplementary" = GEM_CONFIG$supplementary_file
)

for (name in names(files_to_check)) {
  file_path <- files_to_check[[name]]
  if (file.exists(file_path)) {
    cat("  ✓", name, "exists locally (", file.size(file_path), "bytes)\n")
  } else {
    cat("  ✗", name, "not found locally:", file_path, "\n")
    cat("    (Will be downloaded from GitHub when needed)\n")
  }
}

cat("\n=== Test Complete ===\n")

