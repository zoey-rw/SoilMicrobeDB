# Create GTDB-NCBI mapping keys for all GTDB releases (214, 95, and 207)
# This script processes Excel files with GTDB vs NCBI mappings and merges with NCBI taxdump files
# to create comprehensive mapping keys for taxonomy conversion.
#
# Generates:
#   GTDB_NCBI_key.rds - Unified mapping file containing all GTDB versions (214, 95, and 207)
#                       Used by both SPIRE (filters to GTDB207) and SMAG (uses all versions)
#
# Usage:
#   Rscript scripts/create_database/process_genomes/GTDB_mapping/create_all_gtdb_mappings.r
#
# Environment variables:
#   NCBI_TAX_DIR: Directory containing NCBI taxdump files (default: data/genome_database/ncbi_taxonomy)
#   GTDB_MAPPING_OUTPUT: Output file for unified mapping (default: data/genome_database/gtdb_mapping/GTDB_NCBI_key.rds)

# Load configuration and helper functions
# Find the create_database directory
script_dir <- tryCatch({
  # Try to get script path from commandArgs
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg)
    dirname(normalizePath(script_path))
  } else {
    # Fallback: try sys.frame
    dirname(normalizePath(sys.frame(1)$ofile))
  }
}, error = function(e) {
  getwd()
})
if (length(script_dir) == 0 || script_dir == "." || is.na(script_dir)) {
  script_dir <- getwd()
}
# Go up to create_database directory (from process_genomes/GTDB_mapping)
config_path <- file.path(script_dir, "../..", "config.R")
if (!file.exists(config_path)) {
  # Try alternative path from project root
  config_path <- file.path(script_dir, "..", "..", "..", "..", "scripts", "create_database", "config.R")
}
if (!file.exists(config_path)) {
  stop("Cannot find config.R. Please ensure you're running from the project directory.")
}
source(config_path)
source(HELPER_FUNCTIONS)

# Load required libraries
if (!require(readxl, quietly = TRUE)) {
  stop("Package 'readxl' is required. Install with: install.packages('readxl')")
}
if (!require(CHNOSZ, quietly = TRUE)) {
  stop("Package 'CHNOSZ' is required. Install with: install.packages('CHNOSZ')")
}

source_name <- "GTDB_mapping"
log_message("Starting GTDB-NCBI mapping key creation for all releases (214, 95, 207)", source_name = source_name)

# ============================================================================
# Constants
# ============================================================================

# GTDB server base URLs
GTDB_BASE_URL <- "https://data.gtdb.ecogenomic.org/releases"
GTDB_214_URL <- file.path(GTDB_BASE_URL, "release214/214.0/auxillary_files")
GTDB_95_URL <- file.path(GTDB_BASE_URL, "release95/95.0/auxillary_files")
GTDB_207_URL <- file.path(GTDB_BASE_URL, "release207/207.0/auxillary_files")

TAXONOMIC_RANKS <- c("phylum", "class", "order", "family", "genus", "species")
EXCEL_SHEET_NAMES <- 1:6  # Sheets correspond to ranks in TAXONOMIC_RANKS
EXCEL_COL_NAMES <- c("GTDB_taxon", "n_genomes", "n_ncbi_taxa", "NCBI_taxa")
NCBI_TAXON_COLS <- c("NCBI_taxon1", "NCBI_taxon2", "NCBI_taxon3", "NCBI_taxon4", "NCBI_taxon5")

# Filter thresholds for minimum number of genomes
MIN_GENOMES_PHYLUM <- 2
MIN_GENOMES_CLASS <- 3
MIN_GENOMES_ORDER <- 2
MIN_GENOMES_FAMILY <- 1

# Excluded taxa (problematic mappings)
EXCLUDED_TAXA <- list(
  bacteria = c("o__Thermoplasmatales"),
  archaea = c("o__Desulfofervidales", "f__Desulfofervidaceae", "g__Seramator")
)

# Excluded taxa for GTDB 207 (same as general, but defined separately for clarity)
EXCLUDED_TAXA_BACTERIA <- c("o__Thermoplasmatales")
EXCLUDED_TAXA_ARCHAEA <- c("o__Desulfofervidales", "f__Desulfofervidaceae", "g__Seramator")

# NCBI taxdump column names
NCBI_TAXDUMP_COLS <- c("tax_id", "tax_name", "species", "genus", "family",
                       "order", "class", "phylum", "kingdom", "superkingdom", "NA")

# ============================================================================
# Helper Functions
# ============================================================================

# Read Excel file (local, HTTP URL, or remote SSH)
read_excel_file <- function(filepath, source_name = "GTDB_mapping") {
  # Check if local file exists
  if (file.exists(filepath)) {
    return(filepath)
  }
  
  # Check if it's an HTTP/HTTPS URL
  if (grepl("^https?://", filepath)) {
    log_message(paste("Downloading Excel file from URL:", filepath), source_name = source_name)
    temp_file <- tempfile(fileext = ".xlsx")
    
    # Try downloading with download.file (with increased timeout)
    result <- tryCatch({
      # Set timeout to 5 minutes (300 seconds) for large files
      old_timeout <- getOption("timeout")
      options(timeout = 300)
      download.file(filepath, temp_file, mode = "wb", quiet = TRUE)
      options(timeout = old_timeout)
    }, error = function(e) {
      log_message(paste("Download failed, trying with curl:", e$message), source_name = source_name)
      # Fallback to curl with longer timeout (10 minutes)
      system(paste("curl -L --max-time 600 -o", shQuote(temp_file), shQuote(filepath)))
    })
    
    if (file.exists(temp_file) && file.size(temp_file) > 0) {
      return(temp_file)
    } else {
      stop(paste("Failed to download Excel file from URL:", filepath))
    }
  }
  
  # Check if it's a remote SSH path
  if (grepl("^/projectnb", filepath)) {
    log_message(paste("Reading Excel file from remote server via SSH:", filepath), source_name = source_name)
    temp_file <- tempfile(fileext = ".xlsx")
    ssh_cmd <- paste0("ssh zrwerbin@scc2.bu.edu 'cat ", shQuote(filepath), "'")
    result <- system(paste(ssh_cmd, ">", temp_file))
    if (result == 0 && file.exists(temp_file) && file.size(temp_file) > 0) {
      return(temp_file)
    } else {
      stop(paste("Failed to read remote Excel file:", filepath))
    }
  }
  
  stop(paste("Excel file not found:", filepath))
}

# Read NCBI taxdump file (local or remote)
read_ncbi_taxdump <- function(filepath, source_name = "GTDB_mapping") {
  if (file.exists(filepath)) {
    return(data.table::fread(filepath,
                             col.names = NCBI_TAXDUMP_COLS,
                             sep = "|") %>%
           mutate_all(~str_replace_all(., "\t", "")) %>%
           mutate_all(~na_if(., "")) %>%
           select(-`NA`))
  } else if (grepl("^/projectnb", filepath)) {
    log_message(paste("Reading NCBI taxdump from remote:", filepath), source_name = source_name)
    temp_file <- tempfile(fileext = ".dmp")
    ssh_cmd <- paste0("ssh zrwerbin@scc2.bu.edu 'cat ", shQuote(filepath), "'")
    system(paste(ssh_cmd, ">", temp_file))
    if (file.exists(temp_file) && file.size(temp_file) > 0) {
      result <- data.table::fread(temp_file,
                                  col.names = NCBI_TAXDUMP_COLS,
                                  sep = "|") %>%
        mutate_all(~str_replace_all(., "\t", "")) %>%
        mutate_all(~na_if(., "")) %>%
        select(-`NA`)
      unlink(temp_file)
      return(result)
    } else {
      stop(paste("Failed to read remote NCBI taxdump:", filepath))
    }
  } else {
    stop(paste("NCBI taxdump file not found:", filepath))
  }
}

# Read a single Excel sheet for a specific taxonomic rank
read_rank_sheet <- function(filepath, sheet_index, rank_name, col_names) {
  read_xlsx(filepath, sheet = sheet_index, col_names = col_names, skip = 1) %>%
    mutate(rank = rank_name)
}

# Apply domain-specific filters to a rank mapping
apply_domain_filters <- function(mapping_df, rank_name, domain) {
  result <- mapping_df
  
  # Apply minimum genome filters
  if (rank_name == "phylum") {
    result <- result %>% filter(n_genomes > MIN_GENOMES_PHYLUM)
  } else if (rank_name == "class") {
    result <- result %>% filter(n_genomes > MIN_GENOMES_CLASS)
  } else if (rank_name == "order") {
    result <- result %>% filter(n_genomes > MIN_GENOMES_ORDER)
  } else if (rank_name == "family") {
    result <- result %>% filter(n_genomes > MIN_GENOMES_FAMILY)
  }
  
  # Apply domain-specific exclusions
  excluded <- EXCLUDED_TAXA[[domain]]
  if (!is.null(excluded) && length(excluded) > 0) {
    for (excluded_taxon in excluded) {
      result <- result %>% filter(GTDB_taxon != excluded_taxon)
    }
  }
  
  return(result)
}

# Process GTDB mapping Excel file for a domain (bacteria or archaea)
process_gtdb_mapping <- function(mapping_filepath, domain = "bacteria", source_name = "GTDB_mapping") {
  log_message(paste("Processing", domain, "mappings from:", basename(mapping_filepath)), source_name = source_name)
  
  # Read Excel file (local, HTTP URL, or remote SSH)
  excel_path <- read_excel_file(mapping_filepath, source_name = source_name)
  
  # Read all rank sheets
  rank_mappings <- list()
  for (i in seq_along(TAXONOMIC_RANKS)) {
    rank_name <- TAXONOMIC_RANKS[i]
    sheet_index <- EXCEL_SHEET_NAMES[i]
    
    rank_mapping <- read_rank_sheet(excel_path, sheet_index, rank_name, EXCEL_COL_NAMES)
    rank_mapping <- apply_domain_filters(rank_mapping, rank_name, domain)
    rank_mappings[[i]] <- rank_mapping
  }
  
  # Combine all ranks
  combined_mapping <- do.call(rbind, rank_mappings)
  
  return(combined_mapping)
}

# Process GTDB 207 mappings (slightly different filter logic)
process_gtdb_207_mapping <- function(mapping_filepath, domain = "bacteria", source_name = "GTDB_mapping") {
  log_message(paste("Processing", domain, "mappings from:", basename(mapping_filepath)), source_name = source_name)
  
  excel_path <- read_excel_file(mapping_filepath, source_name = source_name)
  
  # Define excluded taxa based on domain
  excluded_taxa <- if (domain == "bacteria") EXCLUDED_TAXA_BACTERIA else EXCLUDED_TAXA_ARCHAEA
  
  # Read each rank sheet
  rank_mappings <- list()
  for (i in seq_along(TAXONOMIC_RANKS)) {
    rank_name <- TAXONOMIC_RANKS[i]
    sheet_index <- EXCEL_SHEET_NAMES[i]
    
    rank_mapping <- read_xlsx(excel_path, sheet = sheet_index, 
                             col_names = EXCEL_COL_NAMES, skip = 1) %>%
      mutate(rank = rank_name)
    
    # Apply filters
    if (rank_name == "phylum") {
      rank_mapping <- rank_mapping %>% filter(n_genomes > MIN_GENOMES_PHYLUM)
    }
    
    # Apply excluded taxa filter
    for (excluded_taxon in excluded_taxa) {
      rank_mapping <- rank_mapping %>% filter(GTDB_taxon != excluded_taxon)
    }
    
    rank_mappings[[i]] <- rank_mapping
  }
  
  return(do.call(rbind, rank_mappings))
}

# Clean and extract NCBI taxon names from comma-separated list
extract_ncbi_taxon_name <- function(mapping_df) {
  mapping_df %>%
    separate(NCBI_taxa, sep = ", ", into = NCBI_TAXON_COLS,
           fill = "right", remove = FALSE) %>%
    mutate(
      NCBI_taxon1 = clean_ncbi(NCBI_taxon1),
      NCBI_taxon2 = clean_ncbi(NCBI_taxon2),
      NCBI_taxon3 = clean_ncbi(NCBI_taxon3),
      NCBI_taxon4 = clean_ncbi(NCBI_taxon4),
      NCBI_taxon5 = clean_ncbi(NCBI_taxon5),
      NCBI_taxon_name = ifelse(NCBI_taxon1 != "", NCBI_taxon1,
                             ifelse(NCBI_taxon2 != "", NCBI_taxon2,
                                    ifelse(NCBI_taxon3 != "", NCBI_taxon3, NA)))
    ) %>%
    filter(!is.na(NCBI_taxon_name))
}

# Transform mapping to long format for merging with NCBI taxdump
transform_to_long_format <- function(mapping_df, gtdb_version) {
  mapping_df %>%
    select(-c(n_genomes, n_ncbi_taxa, all_of(NCBI_TAXON_COLS))) %>%
    distinct(GTDB_taxon, rank, NCBI_taxon_name, .keep_all = TRUE) %>%
    pivot_wider(names_from = "rank", values_from = "NCBI_taxon_name") %>%
    select(-NCBI_taxa) %>%
    distinct(.keep_all = TRUE) %>%
    pivot_longer(cols = all_of(TAXONOMIC_RANKS),
                 values_to = "tax_name", names_to = "rank") %>%
    filter(!is.na(GTDB_taxon) & GTDB_taxon != "NA") %>%
    filter(!is.na(tax_name)) %>%
    mutate(GTDB_version = gtdb_version)
}

# Note: process_gtdb_release function removed - files are now downloaded upfront
# and processed directly in the main script for better efficiency

# Merge mapping with NCBI taxdump
merge_with_ncbi_taxdump <- function(mapping_long, ncbi_taxdump, source_name = "GTDB_mapping") {
  merged <- merge(mapping_long, ncbi_taxdump, all.x = TRUE, by = "tax_name") %>%
    filter(!is.na(tax_id)) %>%
    distinct(.keep_all = TRUE)
  
  return(merged)
}

# Format final mapping key
format_mapping_key <- function(mapping_key) {
  mapping_key %>%
    mutate(
      GTDB_taxon_prefix = GTDB_taxon,
      GTDB_taxon = gsub("s__|g__|f__|o__|c__|p__", "", GTDB_taxon_prefix)
    ) %>%
    select(GTDB_taxon, GTDB_version, GTDB_taxon_prefix, rank, 
           ncbi_tax_id = tax_id, NCBI_tax_name = tax_name, NCBI_version)
}

# Save output file (local or remote)
save_output_file <- function(data, output_file, source_name = "GTDB_mapping") {
  log_message(paste("Saving mapping keys to:", output_file), source_name = source_name)
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir) && !grepl("^/projectnb", output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Save locally or to remote via SSH
  if (grepl("^/projectnb", output_file)) {
    log_message("Saving to remote server via SSH", source_name = source_name)
    temp_file <- tempfile(fileext = ".rds")
    saveRDS(data, temp_file)
    ssh_cmd <- paste0("ssh zrwerbin@scc2.bu.edu 'cat > ", shQuote(output_file), "' < ", temp_file)
    system(ssh_cmd)
    unlink(temp_file)
  } else {
    saveRDS(data, output_file)
  }
}

# ============================================================================
# Main Processing
# ============================================================================

# Set up paths
GTDB_MAPPING_LOCAL_DIR <- file.path(GENOME_DB_DIR, "gtdb_mapping")
dir.create(GTDB_MAPPING_LOCAL_DIR, showWarnings = FALSE, recursive = TRUE)

OUTPUT_FILE <- Sys.getenv("GTDB_MAPPING_OUTPUT", unset = NA)
if (is.na(OUTPUT_FILE)) {
  # Unified mapping file containing all GTDB versions (214, 95, and 207)
  OUTPUT_FILE <- file.path(GTDB_MAPPING_LOCAL_DIR, "GTDB_NCBI_key.rds")
}

# ============================================================================
# Step 1: Download all Excel files upfront
# ============================================================================

log_message("=== Step 1: Downloading all GTDB Excel files ===", source_name = source_name)

# Download all Excel files upfront (more efficient than scattered downloads)
log_message("Downloading GTDB 214 files...", source_name = source_name)
file_214_bac <- read_excel_file(file.path(GTDB_214_URL, "gtdb_vs_ncbi_r214_bacteria.xlsx"), source_name = source_name)
file_214_ar <- read_excel_file(file.path(GTDB_214_URL, "gtdb_vs_ncbi_r214_archaea.xlsx"), source_name = source_name)

log_message("Downloading GTDB 95 files...", source_name = source_name)
file_95_bac <- read_excel_file(file.path(GTDB_95_URL, "gtdb_vs_ncbi_r95_bacteria.xlsx"), source_name = source_name)
file_95_ar <- read_excel_file(file.path(GTDB_95_URL, "gtdb_vs_ncbi_r95_archaea.xlsx"), source_name = source_name)

log_message("Downloading GTDB 207 files...", source_name = source_name)
file_207_bac <- read_excel_file(file.path(GTDB_207_URL, "gtdb_vs_ncbi_r207_bacteria.xlsx"), source_name = source_name)
file_207_ar <- read_excel_file(file.path(GTDB_207_URL, "gtdb_vs_ncbi_r207_archaea.xlsx"), source_name = source_name)

log_message("All Excel files downloaded successfully", source_name = source_name)

# ============================================================================
# Step 2: Process GTDB 214 and 95 mappings
# ============================================================================

log_message("=== Step 2: Processing GTDB 214 and 95 mappings ===", source_name = source_name)

# Process GTDB 214
log_message("Processing GTDB 214 mappings", source_name = source_name)
mapping_214_bacteria <- process_gtdb_mapping(file_214_bac, domain = "bacteria", source_name = source_name)
mapping_214_archaea <- process_gtdb_mapping(file_214_ar, domain = "archaea", source_name = source_name)
gtdb214_mapping <- rbind(mapping_214_bacteria, mapping_214_archaea) %>%
  extract_ncbi_taxon_name() %>%
  mutate(GTDB_version = "GTDB214")
gtdb214_mapping_long <- transform_to_long_format(gtdb214_mapping, "GTDB214")

# Process GTDB 95
log_message("Processing GTDB 95 mappings", source_name = source_name)
mapping_95_bacteria <- process_gtdb_mapping(file_95_bac, domain = "bacteria", source_name = source_name)
mapping_95_archaea <- process_gtdb_mapping(file_95_ar, domain = "archaea", source_name = source_name)
gtdb95_mapping <- rbind(mapping_95_bacteria, mapping_95_archaea) %>%
  extract_ncbi_taxon_name() %>%
  mutate(GTDB_version = "GTDB95")
gtdb95_mapping_long <- transform_to_long_format(gtdb95_mapping, "GTDB95")

# ============================================================================
# Step 3: Load NCBI taxdump files
# ============================================================================

log_message("=== Step 3: Loading NCBI taxdump files ===", source_name = source_name)

ncbi_taxdump_2021_path <- file.path(NCBI_TAX_DIR, "new_taxdump_2021-01-01", "rankedlineage.dmp")
ncbi_taxdump_2023_path <- file.path(NCBI_TAX_DIR, "new_taxdump_2023-12-01", "rankedlineage.dmp")

# Try remote paths if local not found
if (!file.exists(ncbi_taxdump_2021_path)) {
  ncbi_taxdump_2021_path <- file.path(GTDB_REMOTE_BASE, "new_taxdump_2021-01-01", "rankedlineage.dmp")
}
if (!file.exists(ncbi_taxdump_2023_path)) {
  ncbi_taxdump_2023_path <- file.path(GTDB_REMOTE_BASE, "new_taxdump_2023-12-01", "rankedlineage.dmp")
}

ncbi_taxdump_2021 <- read_ncbi_taxdump(ncbi_taxdump_2021_path) %>%
  mutate(NCBI_version = "NCBI_2021-01-01")

ncbi_taxdump_2023 <- read_ncbi_taxdump(ncbi_taxdump_2023_path) %>%
  mutate(NCBI_version = "NCBI_2023-12-01")

# ============================================================================
# Step 4: Merge mappings with NCBI taxdump
# ============================================================================

log_message("=== Step 4: Merging GTDB 214/95 mappings with NCBI taxdump ===", source_name = source_name)

mapping_key_2021_214 <- merge_with_ncbi_taxdump(gtdb214_mapping_long, ncbi_taxdump_2021)
mapping_key_2023_214 <- merge_with_ncbi_taxdump(gtdb214_mapping_long, ncbi_taxdump_2023)
mapping_key_2021_95 <- merge_with_ncbi_taxdump(gtdb95_mapping_long, ncbi_taxdump_2021)
mapping_key_2023_95 <- merge_with_ncbi_taxdump(gtdb95_mapping_long, ncbi_taxdump_2023)

log_message(paste("GTDB 214 mappings - 2021:", nrow(mapping_key_2021_214), "matched"), source_name = source_name)
log_message(paste("GTDB 214 mappings - 2023:", nrow(mapping_key_2023_214), "matched"), source_name = source_name)
log_message(paste("GTDB 95 mappings - 2021:", nrow(mapping_key_2021_95), "matched"), source_name = source_name)
log_message(paste("GTDB 95 mappings - 2023:", nrow(mapping_key_2023_95), "matched"), source_name = source_name)

# Prioritize 2023 over 2021 (remove duplicates from 2021)
mapping_key_2021_214 <- mapping_key_2021_214 %>% 
  filter(!GTDB_taxon %in% mapping_key_2023_214$GTDB_taxon)
mapping_key_2021_95 <- mapping_key_2021_95 %>% 
  filter(!GTDB_taxon %in% mapping_key_2023_95$GTDB_taxon)

# Combine and format final mapping keys
mapping_key_214 <- rbind(mapping_key_2023_214, mapping_key_2021_214)
mapping_key_95 <- rbind(mapping_key_2023_95, mapping_key_2021_95)

gtdb214_mapping_key <- format_mapping_key(mapping_key_214)
gtdb95_mapping_key <- format_mapping_key(mapping_key_95)

# ============================================================================
# Step 5: Process GTDB 207 mappings and combine with 214/95
# ============================================================================

log_message("=== Step 5: Processing GTDB 207 mappings ===", source_name = source_name)

# Process GTDB 207 mappings (using already downloaded files)
gtdb207_ncbi_mapping_bacteria <- process_gtdb_207_mapping(file_207_bac, domain = "bacteria")
gtdb207_ncbi_mapping_archaea <- process_gtdb_207_mapping(file_207_ar, domain = "archaea")

# Combine and process
gtdb207_ncbi_mapping <- rbind(gtdb207_ncbi_mapping_bacteria, gtdb207_ncbi_mapping_archaea) %>%
  extract_ncbi_taxon_name() %>%
  mutate(GTDB_version = "GTDB207")

gtdb207_ncbi_mapping_long <- transform_to_long_format(gtdb207_ncbi_mapping, "GTDB207")

# Merge GTDB 207 mappings with NCBI taxdump (already loaded in Step 3)
log_message("Merging GTDB 207 mappings with NCBI taxdump", source_name = source_name)

gtdb207_ncbi2021_mapping_key <- merge_with_ncbi_taxdump(gtdb207_ncbi_mapping_long, ncbi_taxdump_2021)
gtdb207_ncbi2023_mapping_key <- merge_with_ncbi_taxdump(gtdb207_ncbi_mapping_long, ncbi_taxdump_2023)

log_message(paste("GTDB 207 mappings - 2021:", nrow(gtdb207_ncbi2021_mapping_key), "matched"),
            source_name = source_name)
log_message(paste("GTDB 207 mappings - 2023:", nrow(gtdb207_ncbi2023_mapping_key), "matched"),
            source_name = source_name)

# Format mapping keys
gtdb207_ncbi2021_mapping_simple <- format_mapping_key(gtdb207_ncbi2021_mapping_key)
gtdb207_ncbi2023_mapping_simple <- format_mapping_key(gtdb207_ncbi2023_mapping_key)

# Combine all mappings into unified file
log_message("Combining all GTDB mappings (214, 95, and 207) into unified file", source_name = source_name)

# Combine GTDB 207 mappings (prioritize 2023 over 2021)
gtdb207_mapping_key <- rbind(gtdb207_ncbi2023_mapping_simple, 
                             gtdb207_ncbi2021_mapping_simple %>% 
                               filter(!GTDB_taxon_prefix %in% gtdb207_ncbi2023_mapping_simple$GTDB_taxon_prefix))

# Combine all three GTDB versions
gtdb_ncbi_mapping <- rbind(gtdb214_mapping_key, gtdb95_mapping_key, gtdb207_mapping_key) %>%
  distinct(GTDB_taxon_prefix, NCBI_tax_name, .keep_all = TRUE)

log_message(paste("Total unified mappings:", nrow(gtdb_ncbi_mapping)), source_name = source_name)
log_message(paste("  - GTDB 214:", sum(gtdb_ncbi_mapping$GTDB_version == "GTDB214"), "mappings"), source_name = source_name)
log_message(paste("  - GTDB 95:", sum(gtdb_ncbi_mapping$GTDB_version == "GTDB95"), "mappings"), source_name = source_name)
log_message(paste("  - GTDB 207:", sum(gtdb_ncbi_mapping$GTDB_version == "GTDB207"), "mappings"), source_name = source_name)

# Save unified mapping file
save_output_file(gtdb_ncbi_mapping, OUTPUT_FILE)

log_message("=== GTDB-NCBI mapping key creation completed ===", source_name = source_name)
log_message(paste("Unified mapping file saved to:", OUTPUT_FILE), source_name = source_name)
log_message("This file contains all GTDB versions (214, 95, 207) and is used by both SPIRE and SMAG", source_name = source_name)

