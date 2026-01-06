# Process GEM (Genomes from Earth's Microbiomes) / Nayfach catalog genomes for Soil Microbe Database
# Downloads genomes, processes metadata, and creates Struo2 input file
#
# Usage:
#   Rscript download_nayfach_MAGs.r
#   TEST_MODE=TRUE MAX_GENOMES=10 Rscript download_nayfach_MAGs.r  # Test with 10 genomes

library(tidyverse)
library(data.table)
library(readxl)
library(phyloseq)

# Load centralized configuration
# Find the create_database directory (go up from process_genomes/GEM)
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
# Go up to create_database directory (from process_genomes/GEM)
config_path <- file.path(script_dir, "../..", "config.R")
if (!file.exists(config_path)) {
  # Try alternative path from project root
  config_path <- file.path(script_dir, "..", "..", "..", "..", "scripts", "create_database", "config.R")
}
if (!file.exists(config_path)) {
  stop("Cannot find config.R. Please ensure you're running from the project directory.")
}
source(config_path)

# Get GEM-specific configuration
gem_config <- get_source_config("GEM")
source_name <- "GEM"

# Load helper functions
if (!file.exists(HELPER_FUNCTIONS)) {
  stop(paste("Helper functions not found at:", HELPER_FUNCTIONS))
}
source(HELPER_FUNCTIONS)

log_message("Starting GEM catalog processing", source_name = source_name)

# Create local directory if it doesn't exist
GEM_LOCAL_DIR <- get_source_dir("GEM")
dir.create(GEM_LOCAL_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Step 1: Download and filter GEM metadata
# ============================================================================

log_message("Loading GEM catalog metadata", source_name = source_name)

# Check for local metadata file first, then try to download
GEM_metadata_file <- file.path(GEM_LOCAL_DIR, "genome_metadata.tsv")
GEM_metadata <- NULL

if (file.exists(GEM_metadata_file)) {
  log_message(paste("Using local GEM metadata file:", GEM_metadata_file), source_name = source_name)
  GEM_metadata <- fread(GEM_metadata_file)
} else {
  log_message("Local metadata file not found, attempting to download from portal", source_name = source_name)
  tryCatch({
    GEM_metadata <- fread("https://portal.nersc.gov/GEM/genomes/genome_metadata.tsv")
    # Save a local copy for future use
    log_message(paste("Saving metadata to local file:", GEM_metadata_file), source_name = source_name)
    write_tsv(GEM_metadata, GEM_metadata_file)
  }, error = function(e) {
    stop(paste("Failed to download GEM metadata and local file not found.\n",
               "Please either:\n",
               "  1. Download the metadata file manually and place it at:", GEM_metadata_file, "\n",
               "  2. Ensure network access to https://portal.nersc.gov/GEM/genomes/genome_metadata.tsv\n",
               "Error:", e$message))
  })
}

if (is.null(GEM_metadata) || nrow(GEM_metadata) == 0) {
  stop("GEM metadata is empty or could not be loaded")
}

log_message("Filtering for soil genomes with quality thresholds", source_name = source_name)
soil_GEM_metadata <- GEM_metadata %>%
  filter(ecosystem_type == "Soil" & completeness > MIN_COMPLETENESS & contamination < MAX_CONTAMINATION)

log_message(paste("Found", nrow(soil_GEM_metadata), "soil GEM genomes meeting quality criteria"), source_name = source_name)

# Apply test mode limit if enabled
if (is_test_mode() && !is.na(MAX_GENOMES)) {
  log_message(paste("TEST MODE: Limiting to", MAX_GENOMES, "genomes"), source_name = source_name)
  soil_GEM_metadata <- soil_GEM_metadata %>% head(MAX_GENOMES)
}

# Create download links
soil_GEM_metadata$download_links <- paste0(gem_config$download_base_url, soil_GEM_metadata$genome_id, ".fna.gz")

# ============================================================================
# Step 2: Check for existing downloaded genomes
# ============================================================================

log_message("Checking for existing downloaded genomes", source_name = source_name)
soil_GEM_downloaded <- data.frame()

# Check local directory
local_genomes <- read_in_genomes(gem_config$local_genome_dir, pattern = gem_config$genome_pattern)
if (nrow(local_genomes) > 0) {
  soil_GEM_downloaded <- rbind(soil_GEM_downloaded, local_genomes)
  log_message(paste("Found", nrow(local_genomes), "genomes in local directory"), source_name = source_name)
}

# Check genome directory (if different from local)
if ("genome_dir" %in% names(gem_config) && gem_config$genome_dir != gem_config$local_genome_dir) {
  remote_genomes <- read_in_genomes(gem_config$genome_dir, pattern = gem_config$genome_pattern)
  if (nrow(remote_genomes) > 0) {
    soil_GEM_downloaded <- rbind(soil_GEM_downloaded, remote_genomes)
    log_message(paste("Found", nrow(remote_genomes), "genomes in genome directory"), source_name = source_name)
  }
}

# Remove duplicates if genomes exist in both local and remote
if (nrow(soil_GEM_downloaded) > 0) {
  soil_GEM_downloaded <- soil_GEM_downloaded %>%
    distinct(filename, .keep_all = TRUE)
  log_message(paste("Total unique genomes found:", nrow(soil_GEM_downloaded)), source_name = source_name)
}

# ============================================================================
# Step 3: Download missing genomes (if not in test mode)
# ============================================================================

if (!is_test_mode()) {
  missing_genomes <- soil_GEM_metadata %>%
    filter(!genome_id %in% soil_GEM_downloaded$user_genome)
  
  if (nrow(missing_genomes) > 0) {
    log_message(paste("Downloading", nrow(missing_genomes), "missing genomes"), source_name = source_name)
    dl_links <- missing_genomes$download_links %>% na.omit() %>% unique()
    
    for (i in seq_along(dl_links)) {
      if (basename(dl_links[i]) == "NA_genomic.fna.gz") next()
      dest_file <- file.path(gem_config$local_genome_dir, basename(dl_links[i]))
      if (!file.exists(dest_file)) {
        log_message(paste("Downloading:", basename(dl_links[i])), source_name = source_name)
        tryCatch({
          download.file(dl_links[i], dest_file, mode = "wb")
        }, error = function(e) {
          log_message(paste("Failed to download", basename(dl_links[i]), ":", e$message), source_name = source_name)
        })
      }
    }
    
    # Re-read downloaded genomes
    soil_GEM_downloaded <- read_in_genomes(gem_config$local_genome_dir, pattern = gem_config$genome_pattern)
  }
} else {
  log_message("TEST MODE: Skipping genome downloads", source_name = source_name)
  # Create placeholder filepaths for test mode
  soil_GEM_downloaded <- data.frame(
    user_genome = soil_GEM_metadata$genome_id,
    filepath = file.path(gem_config$local_genome_dir, paste0(soil_GEM_metadata$genome_id, ".fna.gz")),
    filename = paste0(soil_GEM_metadata$genome_id, ".fna.gz"),
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# Step 4: Load supplementary material from Nayfach et al.
# ============================================================================

log_message("Loading supplementary material from Nayfach et al.", source_name = source_name)

# Get supplementary file path
supplementary_file_config <- gem_config$supplementary_file
local_supplementary <- file.path(GEM_LOCAL_DIR, "41587_2020_718_MOESM3_ESM.xlsx")

# Use config file path if it exists, otherwise try local path
if (file.exists(supplementary_file_config)) {
  SUPPLEMENTARY_FILE <- supplementary_file_config
} else if (file.exists(local_supplementary)) {
  SUPPLEMENTARY_FILE <- local_supplementary
} else {
  stop(paste("Nayfach supplementary file not found. Tried:",
             "\n  Config path:", supplementary_file_config,
             "\n  Local path:", local_supplementary,
             "\nPlease ensure file is available locally or mount remote server."))
}

# Read supplementary material sheets
nayfach_s2 <- readxl::read_xlsx(SUPPLEMENTARY_FILE, sheet = 3)
nayfach_s5 <- readxl::read_xlsx(SUPPLEMENTARY_FILE, sheet = 6)

# ============================================================================
# Step 5: Merge metadata with downloaded genomes
# ============================================================================

log_message("Merging metadata with downloaded genomes", source_name = source_name)
# Match genomes by filename (user_genome from read_in_genomes matches genome_id)
soil_GEM_metadata <- soil_GEM_metadata %>%
  mutate(filename = paste0(genome_id, ".fna.gz")) %>%
  left_join(soil_GEM_downloaded %>% select(filename, filepath), by = "filename")

# Ensure filepath exists after merge (create if missing or fill NA values)
# Fill any NA filepaths (use local directory as default for missing genomes)
na_filepath <- is.na(soil_GEM_metadata$filepath)
if (any(na_filepath)) {
  soil_GEM_metadata$filepath[na_filepath] <- file.path(gem_config$local_genome_dir, paste0(soil_GEM_metadata$genome_id[na_filepath], ".fna.gz"))
}

# Match to OTU and taxonomy
soil_GEM_metadata$OTU_ID <- nayfach_s2[match(soil_GEM_metadata$genome_id, nayfach_s2$genome_id),]$species_id
soil_GEM_metadata$taxonomy2 <- nayfach_s5[match(soil_GEM_metadata$OTU_ID, nayfach_s5$otu_id),]$gtdb_taxonomy
soil_GEM_metadata$classified_rank <- nayfach_s5[match(soil_GEM_metadata$OTU_ID, nayfach_s5$otu_id),]$classified_rank

# Parse GTDB taxonomy using phyloseq function
log_message("Parsing GTDB taxonomy", source_name = source_name)
tax_split <- lapply(soil_GEM_metadata$taxonomy2, phyloseq::parse_taxonomy_qiime) %>% do.call(rbind, .)
soil_GEM_metadata <- cbind.data.frame(soil_GEM_metadata, tax_split)

# ============================================================================
# Step 6: Prepare taxonomy data for mapping
# ============================================================================

log_message("Preparing taxonomy data for NCBI mapping", source_name = source_name)
soil_GEM_long <- soil_GEM_metadata %>%
  pivot_longer(cols = c("Rank1", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               names_to = "rank", values_to = "name",
               values_drop_na = FALSE) %>%
  mutate(rank = tolower(rank)) %>%
  filter(rank == classified_rank)

soil_GEM_simple <- soil_GEM_long %>%
  select(c("genome_id", "metagenome_id", "otu_id",
																								 "ecosystem", "ecosystem_category", "ecosystem_type", "download_links", "OTU_ID", "taxonomy2",
           "classified_rank", "name", "filepath"))

# Determine higher classification for fallback taxonomy matching
soil_GEM_metadata$higher_classification <- ifelse(soil_GEM_metadata$classified_rank == "species", soil_GEM_metadata$Genus,
                                                   ifelse(soil_GEM_metadata$classified_rank == "genus", soil_GEM_metadata$Family,
                                                          ifelse(soil_GEM_metadata$classified_rank == "family", soil_GEM_metadata$Order,
                                                                 ifelse(soil_GEM_metadata$classified_rank == "order", soil_GEM_metadata$Class,
                                                                        ifelse(soil_GEM_metadata$classified_rank == "class", soil_GEM_metadata$Phylum, NA)))))

# Merge additional columns
soil_GEM_simple <- merge(soil_GEM_simple, 
                         soil_GEM_metadata[,c("genome_id", "taxonomy2", "higher_classification", "Rank1", "Phylum", "Species", "filepath")], 
                         by = "genome_id", all.x = TRUE, suffixes = c("", ".y")) %>%
  distinct(.keep_all = TRUE)

# Remove duplicate filepath column if merge created one
if ("filepath.y" %in% names(soil_GEM_simple)) {
  soil_GEM_simple$filepath <- ifelse(is.na(soil_GEM_simple$filepath), soil_GEM_simple$filepath.y, soil_GEM_simple$filepath)
  soil_GEM_simple$filepath.y <- NULL
}

# Ensure filepath exists (should already be there, but just in case)
if (!"filepath" %in% names(soil_GEM_simple) || any(is.na(soil_GEM_simple$filepath))) {
  na_filepath <- is.na(soil_GEM_simple$filepath) | !"filepath" %in% names(soil_GEM_simple)
  if (any(na_filepath) || !"filepath" %in% names(soil_GEM_simple)) {
    if (!"filepath" %in% names(soil_GEM_simple)) {
      soil_GEM_simple$filepath <- file.path(gem_config$local_genome_dir, paste0(soil_GEM_simple$genome_id, ".fna.gz"))
    } else {
      soil_GEM_simple$filepath[is.na(soil_GEM_simple$filepath)] <- 
        file.path(gem_config$local_genome_dir, paste0(soil_GEM_simple$genome_id[is.na(soil_GEM_simple$filepath)], ".fna.gz"))
    }
  }
}

# Clean higher classification and name (remove GTDB prefixes)
# Note: name from pivot_longer is already clean (parse_taxonomy_qiime removes prefixes)
# But we clean it anyway to be safe
soil_GEM_simple$higher_classification <- gsub("^[a-z]__", "", soil_GEM_simple$higher_classification)
soil_GEM_simple$name_clean <- gsub("^[a-z]__", "", soil_GEM_simple$name)

# ============================================================================
# Step 7: Load GTDB metadata files (for Strategy 1 direct taxonomy mapping)
# ============================================================================

log_message("Loading GTDB 95 and 214 metadata for direct taxonomy mapping", source_name = source_name)

# Function to load GTDB metadata (local or remote)
# Loads both bacterial and archaeal metadata files for the specified GTDB version
# Note: The unified mapping key (GTDB_NCBI_key.rds) is sufficient for matching all genomes.
# Metadata files only provide additional precision for Strategy 1 (direct full-string matching).
# Strategy 2 (mapping key) handles all taxa.
load_gtdb_metadata <- function(version, source_name) {
  if (version == "95") {
    bac_file <- file.path(GENOME_DB_DIR, "bac120_metadata_r95.tsv")
    ar_file <- file.path(GENOME_DB_DIR, "ar53_metadata_r95.tsv")
    remote_bac <- file.path(GTDB_REMOTE_BASE, "bac120_metadata_r95.tsv")
    remote_ar <- file.path(GTDB_REMOTE_BASE, "ar53_metadata_r95.tsv")
  } else if (version == "214") {
    bac_file <- file.path(GENOME_DB_DIR, "bac120_metadata_r214.tsv")
    ar_file <- file.path(GENOME_DB_DIR, "ar53_metadata_r214.tsv")
    remote_bac <- file.path(GTDB_REMOTE_BASE, "bac120_metadata_r214.tsv")
    remote_ar <- file.path(GTDB_REMOTE_BASE, "ar53_metadata_r214.tsv")
  } else {
    return(NULL)
  }
  
  result_list <- list()
  
  # Try local first - bacterial
  if (file.exists(bac_file)) {
    log_message(paste("Loading GTDB", version, "bacterial metadata from local file"), source_name = source_name)
    gtdb_bac <- fread(bac_file) %>%
      select(accession, gtdb_taxonomy, ncbi_taxid) %>%
      mutate(source = paste0("GTDB_", version))
    result_list$bac <- gtdb_bac
  }
  
  # Try local - archaeal
  if (file.exists(ar_file)) {
    log_message(paste("Loading GTDB", version, "archaeal metadata from local file"), source_name = source_name)
    gtdb_ar <- fread(ar_file) %>%
      select(accession, gtdb_taxonomy, ncbi_taxid) %>%
      mutate(source = paste0("GTDB_", version))
    result_list$ar <- gtdb_ar
  }
  
  # Combine results
  if (length(result_list) > 0) {
    return(do.call(rbind, result_list))
  }
  
  return(NULL)
}

# Load GTDB 95 and 214 metadata
gtdb_95_metadata <- load_gtdb_metadata("95", source_name)
gtdb_214_metadata <- load_gtdb_metadata("214", source_name)

if (!is.null(gtdb_95_metadata)) {
  log_message(paste("Loaded GTDB 95 metadata:", nrow(gtdb_95_metadata), "genomes (Strategy 1 enabled)"), source_name = source_name)
} else {
  log_message("GTDB 95 metadata not available (will use mapping key only - this is sufficient)", source_name = source_name)
}

if (!is.null(gtdb_214_metadata)) {
  log_message(paste("Loaded GTDB 214 metadata:", nrow(gtdb_214_metadata), "genomes (Strategy 1 enabled)"), source_name = source_name)
} else {
  log_message("GTDB 214 metadata not available (will use mapping key only - this is sufficient)", source_name = source_name)
}

# ============================================================================
# Step 8: Load GTDB-NCBI mapping file
# ============================================================================

log_message("Loading GTDB-NCBI mapping file", source_name = source_name)
if (!file.exists(GTDB_NCBI_MAPPING_FILE)) {
  stop(paste("GTDB-NCBI mapping file not found:", GTDB_NCBI_MAPPING_FILE,
             "\nPlease ensure remote server is mounted or copy file locally."))
}
ncbi_gtdb_mapping_all <- readRDS(GTDB_NCBI_MAPPING_FILE)

# Filter to GTDB 95, 214, and 207 mappings
# Note: Some taxa may be in GTDB 214 or 207 but not in GTDB 95
mapping_key_2021_gtdb95 <- ncbi_gtdb_mapping_all[
  ncbi_gtdb_mapping_all$GTDB_version == "GTDB95" & 
  ncbi_gtdb_mapping_all$NCBI_version == "NCBI_2021-01-01", ]
mapping_key_2023_gtdb95 <- ncbi_gtdb_mapping_all[
  ncbi_gtdb_mapping_all$GTDB_version == "GTDB95" & 
  ncbi_gtdb_mapping_all$NCBI_version == "NCBI_2023-12-01", ]

mapping_key_2021_gtdb214 <- ncbi_gtdb_mapping_all[
  ncbi_gtdb_mapping_all$GTDB_version == "GTDB214" & 
  ncbi_gtdb_mapping_all$NCBI_version == "NCBI_2021-01-01", ]
mapping_key_2023_gtdb214 <- ncbi_gtdb_mapping_all[
  ncbi_gtdb_mapping_all$GTDB_version == "GTDB214" & 
  ncbi_gtdb_mapping_all$NCBI_version == "NCBI_2023-12-01", ]

mapping_key_2021_gtdb207 <- ncbi_gtdb_mapping_all[
  ncbi_gtdb_mapping_all$GTDB_version == "GTDB207" & 
  ncbi_gtdb_mapping_all$NCBI_version == "NCBI_2021-01-01", ]
mapping_key_2023_gtdb207 <- ncbi_gtdb_mapping_all[
  ncbi_gtdb_mapping_all$GTDB_version == "GTDB207" & 
  ncbi_gtdb_mapping_all$NCBI_version == "NCBI_2023-12-01", ]

log_message(paste("Loaded GTDB 95 mappings - 2021:", nrow(mapping_key_2021_gtdb95), "entries"), source_name = source_name)
log_message(paste("Loaded GTDB 95 mappings - 2023:", nrow(mapping_key_2023_gtdb95), "entries"), source_name = source_name)
log_message(paste("Loaded GTDB 214 mappings - 2021:", nrow(mapping_key_2021_gtdb214), "entries"), source_name = source_name)
log_message(paste("Loaded GTDB 214 mappings - 2023:", nrow(mapping_key_2023_gtdb214), "entries"), source_name = source_name)
log_message(paste("Loaded GTDB 207 mappings - 2021:", nrow(mapping_key_2021_gtdb207), "entries"), source_name = source_name)
log_message(paste("Loaded GTDB 207 mappings - 2023:", nrow(mapping_key_2023_gtdb207), "entries"), source_name = source_name)

# Clean GTDB taxon names (remove prefixes)
mapping_key_2021_gtdb95$GTDB_taxon_clean <- gsub("s__|g__|f__|o__|c__|p__", "", mapping_key_2021_gtdb95$GTDB_taxon)
mapping_key_2023_gtdb95$GTDB_taxon_clean <- gsub("s__|g__|f__|o__|c__|p__", "", mapping_key_2023_gtdb95$GTDB_taxon)
mapping_key_2021_gtdb214$GTDB_taxon_clean <- gsub("s__|g__|f__|o__|c__|p__", "", mapping_key_2021_gtdb214$GTDB_taxon)
mapping_key_2023_gtdb214$GTDB_taxon_clean <- gsub("s__|g__|f__|o__|c__|p__", "", mapping_key_2023_gtdb214$GTDB_taxon)
mapping_key_2021_gtdb207$GTDB_taxon_clean <- gsub("s__|g__|f__|o__|c__|p__", "", mapping_key_2021_gtdb207$GTDB_taxon)
mapping_key_2023_gtdb207$GTDB_taxon_clean <- gsub("s__|g__|f__|o__|c__|p__", "", mapping_key_2023_gtdb207$GTDB_taxon)

# ============================================================================
# Step 9: Map to NCBI taxids using multiple strategies
# ============================================================================

log_message("Mapping GTDB taxonomy to NCBI taxids", source_name = source_name)

# Strategy 1: Direct mapping from GTDB metadata using full taxonomy string (highest priority)
# Try exact match first, then try cleaned taxonomy (removing trailing empty fields)
if (!is.null(gtdb_95_metadata)) {
  # Exact match
  soil_GEM_simple$NCBI_TaxID_genome_95 <- gtdb_95_metadata[
    match(soil_GEM_simple$taxonomy2, gtdb_95_metadata$gtdb_taxonomy),]$ncbi_taxid
  
  # For unmatched, try cleaned taxonomy (remove trailing empty fields like ;s__ or ;g__;s__)
  unmatched_idx <- which(is.na(soil_GEM_simple$NCBI_TaxID_genome_95))
  if (length(unmatched_idx) > 0) {
    for (i in unmatched_idx) {
      tax_clean <- gsub(";s__$", "", soil_GEM_simple$taxonomy2[i])
      tax_clean <- gsub(";g__;s__$", "", tax_clean)
      tax_clean <- gsub(";g__$", "", tax_clean)
      
      matches <- gtdb_95_metadata[gtdb_95_metadata$gtdb_taxonomy == tax_clean,]
      if (nrow(matches) > 0) {
        # Take first match (should be unique)
        soil_GEM_simple$NCBI_TaxID_genome_95[i] <- matches$ncbi_taxid[1]
      }
    }
  }
} else {
  soil_GEM_simple$NCBI_TaxID_genome_95 <- NA
}

if (!is.null(gtdb_214_metadata)) {
  # Exact match
  soil_GEM_simple$NCBI_TaxID_genome_214 <- gtdb_214_metadata[
    match(soil_GEM_simple$taxonomy2, gtdb_214_metadata$gtdb_taxonomy),]$ncbi_taxid
  
  # For unmatched, try cleaned taxonomy
  unmatched_idx <- which(is.na(soil_GEM_simple$NCBI_TaxID_genome_214))
  if (length(unmatched_idx) > 0) {
    for (i in unmatched_idx) {
      tax_clean <- gsub(";s__$", "", soil_GEM_simple$taxonomy2[i])
      tax_clean <- gsub(";g__;s__$", "", tax_clean)
      tax_clean <- gsub(";g__$", "", tax_clean)
      
      matches <- gtdb_214_metadata[gtdb_214_metadata$gtdb_taxonomy == tax_clean,]
      if (nrow(matches) > 0) {
        soil_GEM_simple$NCBI_TaxID_genome_214[i] <- matches$ncbi_taxid[1]
      }
    }
  }
} else {
  soil_GEM_simple$NCBI_TaxID_genome_214 <- NA
}

# Strategy 2: Direct mapping from GTDB mapping keys (try 95, then 214, then 207)
# Match by both name and rank (the unified mapping key includes rank information)
# Priority: GTDB 95 (2023) > GTDB 95 (2021) > GTDB 214 (2023) > GTDB 214 (2021) > GTDB 207 (2023) > GTDB 207 (2021)
soil_GEM_simple$NCBI_TaxID_custom_95_NCBI2023 <- NA
soil_GEM_simple$NCBI_TaxID_custom_95_NCBI2021 <- NA
soil_GEM_simple$NCBI_TaxID_custom_214_NCBI2023 <- NA
soil_GEM_simple$NCBI_TaxID_custom_214_NCBI2021 <- NA
soil_GEM_simple$NCBI_TaxID_custom_207_NCBI2023 <- NA
soil_GEM_simple$NCBI_TaxID_custom_207_NCBI2021 <- NA

for (i in 1:nrow(soil_GEM_simple)) {
  name_clean <- soil_GEM_simple$name_clean[i]
  rank <- soil_GEM_simple$classified_rank[i]
  
  # Helper function to match by name and rank
  match_by_name_rank <- function(mapping_key, name, rank_val) {
    # Try exact match by name and rank first
    matches <- mapping_key[
      mapping_key$GTDB_taxon_clean == name & 
      tolower(mapping_key$rank) == tolower(rank_val), ]
    if (nrow(matches) > 0) {
      return(matches$ncbi_tax_id[1])
    }
    # Fallback: match by name only
    matches <- mapping_key[mapping_key$GTDB_taxon_clean == name, ]
    if (nrow(matches) > 0) {
      return(matches$ncbi_tax_id[1])
    }
    return(NA)
  }
  
  # Try GTDB 95 (2023)
  if (is.na(soil_GEM_simple$NCBI_TaxID_custom_95_NCBI2023[i])) {
    soil_GEM_simple$NCBI_TaxID_custom_95_NCBI2023[i] <- match_by_name_rank(mapping_key_2023_gtdb95, name_clean, rank)
  }
  
  # Try GTDB 95 (2021)
  if (is.na(soil_GEM_simple$NCBI_TaxID_custom_95_NCBI2021[i])) {
    soil_GEM_simple$NCBI_TaxID_custom_95_NCBI2021[i] <- match_by_name_rank(mapping_key_2021_gtdb95, name_clean, rank)
  }
  
  # Try GTDB 214 (2023)
  if (is.na(soil_GEM_simple$NCBI_TaxID_custom_214_NCBI2023[i])) {
    soil_GEM_simple$NCBI_TaxID_custom_214_NCBI2023[i] <- match_by_name_rank(mapping_key_2023_gtdb214, name_clean, rank)
  }
  
  # Try GTDB 214 (2021)
  if (is.na(soil_GEM_simple$NCBI_TaxID_custom_214_NCBI2021[i])) {
    soil_GEM_simple$NCBI_TaxID_custom_214_NCBI2021[i] <- match_by_name_rank(mapping_key_2021_gtdb214, name_clean, rank)
  }
  
  # Try GTDB 207 (2023)
  if (is.na(soil_GEM_simple$NCBI_TaxID_custom_207_NCBI2023[i])) {
    soil_GEM_simple$NCBI_TaxID_custom_207_NCBI2023[i] <- match_by_name_rank(mapping_key_2023_gtdb207, name_clean, rank)
  }
  
  # Try GTDB 207 (2021)
  if (is.na(soil_GEM_simple$NCBI_TaxID_custom_207_NCBI2021[i])) {
    soil_GEM_simple$NCBI_TaxID_custom_207_NCBI2021[i] <- match_by_name_rank(mapping_key_2021_gtdb207, name_clean, rank)
  }
}

# Strategy 3: Higher classification mapping (fallback)
# Match higher classification by name across all GTDB versions
soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2023 <- NA
soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2021 <- NA

for (i in 1:nrow(soil_GEM_simple)) {
  if (!is.na(soil_GEM_simple$higher_classification[i])) {
    higher_clean <- soil_GEM_simple$higher_classification[i]
    
    # Determine expected rank for higher classification
    rank <- soil_GEM_simple$classified_rank[i]
    expected_higher_rank <- ifelse(rank == "species", "genus",
                                  ifelse(rank == "genus", "family",
                                         ifelse(rank == "family", "order",
                                                ifelse(rank == "order", "class",
                                                       ifelse(rank == "class", "phylum", NA)))))
    
    # Helper function to match higher classification
    match_higher <- function(mapping_key, name, expected_rank) {
      if (!is.na(expected_rank)) {
        matches <- mapping_key[
          mapping_key$GTDB_taxon_clean == name & 
          tolower(mapping_key$rank) == tolower(expected_rank), ]
        if (nrow(matches) > 0) {
          return(matches$ncbi_tax_id[1])
        }
      }
      # Fallback: match by name only
      matches <- mapping_key[mapping_key$GTDB_taxon_clean == name, ]
      if (nrow(matches) > 0) {
        return(matches$ncbi_tax_id[1])
      }
      return(NA)
    }
    
    # Try GTDB 95 (2023)
    if (is.na(soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2023[i])) {
      soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2023[i] <- match_higher(mapping_key_2023_gtdb95, higher_clean, expected_higher_rank)
    }
    
    # Try GTDB 95 (2021)
    if (is.na(soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2021[i])) {
      soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2021[i] <- match_higher(mapping_key_2021_gtdb95, higher_clean, expected_higher_rank)
    }
    
    # Also try GTDB 214 and 207 for higher classification (in case they're not in GTDB 95)
    if (is.na(soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2023[i])) {
      # Try GTDB 214 (2023)
      higher_match <- match_higher(mapping_key_2023_gtdb214, higher_clean, expected_higher_rank)
      if (!is.na(higher_match)) {
        soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2023[i] <- higher_match
      } else {
        # Try GTDB 207 (2023)
        higher_match <- match_higher(mapping_key_2023_gtdb207, higher_clean, expected_higher_rank)
        if (!is.na(higher_match)) {
          soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2023[i] <- higher_match
        }
      }
    }
  }
}

# Combine all mapping strategies (prioritize in order)
# Priority: GTDB 95 metadata > GTDB 214 metadata > GTDB 95 mapping (2023) > GTDB 95 mapping (2021) > 
#           GTDB 214 mapping (2023) > GTDB 214 mapping (2021) > GTDB 207 mapping (2023) > GTDB 207 mapping (2021) >
#           higher classification 2023 > higher classification 2021
soil_GEM_simple$NCBI_TaxID <- ifelse(!is.na(soil_GEM_simple$NCBI_TaxID_genome_95), 
                                      soil_GEM_simple$NCBI_TaxID_genome_95,
                                      ifelse(!is.na(soil_GEM_simple$NCBI_TaxID_genome_214), 
                                             soil_GEM_simple$NCBI_TaxID_genome_214,
                                             ifelse(!is.na(soil_GEM_simple$NCBI_TaxID_custom_95_NCBI2023), 
                                                    soil_GEM_simple$NCBI_TaxID_custom_95_NCBI2023,
                                                    ifelse(!is.na(soil_GEM_simple$NCBI_TaxID_custom_95_NCBI2021), 
                                                           soil_GEM_simple$NCBI_TaxID_custom_95_NCBI2021,
                                                           ifelse(!is.na(soil_GEM_simple$NCBI_TaxID_custom_214_NCBI2023), 
                                                                  soil_GEM_simple$NCBI_TaxID_custom_214_NCBI2023,
                                                                  ifelse(!is.na(soil_GEM_simple$NCBI_TaxID_custom_214_NCBI2021), 
                                                                         soil_GEM_simple$NCBI_TaxID_custom_214_NCBI2021,
                                                                         ifelse(!is.na(soil_GEM_simple$NCBI_TaxID_custom_207_NCBI2023), 
                                                                                soil_GEM_simple$NCBI_TaxID_custom_207_NCBI2023,
                                                                                ifelse(!is.na(soil_GEM_simple$NCBI_TaxID_custom_207_NCBI2021), 
                                                                                       soil_GEM_simple$NCBI_TaxID_custom_207_NCBI2021,
                                                                                       ifelse(!is.na(soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2023), 
                                                                                              soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2023,
                                                                                              ifelse(!is.na(soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2021), 
                                                                                                     soil_GEM_simple$higher_NCBI_TaxID_custom_95_NCBI2021,
                                                                                                     NA))))))))))

# Report mapping success
mapped_count <- sum(!is.na(soil_GEM_simple$NCBI_TaxID))
log_message(paste("Successfully mapped", mapped_count, "of", nrow(soil_GEM_simple), "genomes to NCBI taxids"), source_name = source_name)

# ============================================================================
# Step 10: Prepare final output
# ============================================================================

log_message("Preparing final output", source_name = source_name)
ready_genomes_nayfach <- soil_GEM_simple %>%
	filter(!is.na(NCBI_TaxID)) %>%
  mutate(source = "Nayfach_MAGs", 
         is_published = "Y",
         kingdom = gsub("^[a-z]__", "", Rank1),
         phylum = gsub("^[a-z]__", "", Phylum),
         species = gsub("^[a-z]__", "", Species)) %>%
	select(ncbi_species_taxid = NCBI_TaxID,
				 accession = genome_id,
         ncbi_organism_name = name_clean,
         kingdom,
         phylum,
         species,
				 fasta_file_path = filepath,
         ncbi_taxonomy = taxonomy2,
				 source,
				 is_published)

# Clean up intermediate variables (not needed after final output is created)
# Only remove variables that exist to avoid errors
vars_to_remove <- c("gtdb_95_metadata", "gtdb_214_metadata", 
                    "mapping_key_2021_gtdb95", "mapping_key_2023_gtdb95",
                    "mapping_key_2021_gtdb214", "mapping_key_2023_gtdb214",
                    "mapping_key_2021_gtdb207", "mapping_key_2023_gtdb207",
                    "ncbi_gtdb_mapping_all")
existing_vars <- vars_to_remove[vars_to_remove %in% ls()]
if (length(existing_vars) > 0) {
  rm(list = existing_vars)
}

# ============================================================================
# Step 11: Write output file
# ============================================================================

log_message(paste("Writing output to:", gem_config$output_file), source_name = source_name)
write_tsv(ready_genomes_nayfach, gem_config$output_file)

log_message(paste("Successfully processed", nrow(ready_genomes_nayfach), "GEM genomes"), source_name = source_name)
log_message("GEM catalog processing complete", source_name = source_name)
