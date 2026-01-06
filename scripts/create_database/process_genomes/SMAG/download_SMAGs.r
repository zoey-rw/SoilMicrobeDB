# Process SMAG (Soil Metagenome-Assembled Genomes) for Soil Microbe Database
# Downloads genomes, processes metadata, and creates Struo2 input file
#
# Usage:
#   Rscript download_SMAGs.r
#   TEST_MODE=TRUE MAX_GENOMES=10 Rscript download_SMAGs.r  # Test with 10 genomes

library(tidyverse)
library(readxl)

# Find the create_database directory (go up from process_genomes/SMAG)
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
# Go up to create_database directory (from process_genomes/SMAG)
config_path <- file.path(script_dir, "../..", "config.R")
if (!file.exists(config_path)) {
  # Try alternative path from project root
  config_path <- file.path(script_dir, "..", "..", "..", "..", "scripts", "create_database", "config.R")
}
if (!file.exists(config_path)) {
  stop("Cannot find config.R. Please ensure you're running from the project directory.")
}
source(config_path)

# Get SMAG-specific configuration
smag_config <- get_source_config("SMAG")
source_name <- "SMAG"

# Load helper functions
if (!file.exists(HELPER_FUNCTIONS)) {
  stop(paste("Helper functions not found at:", HELPER_FUNCTIONS))
}
source(HELPER_FUNCTIONS)

log_message("Starting SMAG processing", source_name = source_name)

# Check for required files and download if missing
if (!file.exists(smag_config$metadata_file)) {
  log_message(paste("SMAG metadata file not found locally. Downloading from figshare..."), source_name = source_name)
  # Create directory if it doesn't exist
  metadata_dir <- dirname(smag_config$metadata_file)
  if (!dir.exists(metadata_dir)) {
    dir.create(metadata_dir, showWarnings = FALSE, recursive = TRUE)
  }
  # Download from figshare - the download is a ZIP file containing multiple Excel files
  download_url <- "https://springernature.figshare.com/ndownloader/articles/23298791/versions/1"
  zip_file <- file.path(metadata_dir, "smag_data.zip")
  tryCatch({
    # Download ZIP file using curl
    curl_cmd <- paste0("curl -L -o ", shQuote(zip_file), 
                       " -H 'User-Agent: Mozilla/5.0' ", shQuote(download_url))
    result <- system(curl_cmd)
    if (result != 0 || !file.exists(zip_file) || file.size(zip_file) == 0) {
      stop("Failed to download ZIP file")
    }
    log_message("Extracting Excel file from ZIP archive...", source_name = source_name)
    # Extract the specific Excel file we need (Supplementary Data 2.xlsx)
    unzip(zip_file, files = basename(smag_config$metadata_file), exdir = metadata_dir, overwrite = TRUE)
    # Clean up ZIP file
    unlink(zip_file)
    if (file.exists(smag_config$metadata_file) && file.size(smag_config$metadata_file) > 0) {
      log_message(paste("Successfully downloaded and extracted SMAG metadata file to:", smag_config$metadata_file), source_name = source_name)
    } else {
      stop("Extracted file is empty or does not exist")
    }
  }, error = function(e) {
    if (file.exists(zip_file)) unlink(zip_file)
    stop(paste("Failed to download SMAG metadata file:", e$message,
               "\nURL:", download_url,
               "\nPlease download manually from: https://springernature.figshare.com/articles/dataset/A_Genomic_Catalogue_of_Soil_Microbiomes_Boosts_Mining_of_Biodiversity_and_Genetic_Resources/23298791"))
  })
}

# Load GTDB-NCBI mapping (contains all GTDB versions)
log_message("Loading GTDB-NCBI mapping", source_name = source_name)
if (!file.exists(GTDB_NCBI_MAPPING_FILE)) {
  stop(paste("GTDB-NCBI mapping file not found:", GTDB_NCBI_MAPPING_FILE,
             "\nPlease ensure remote server is mounted or copy file locally."))
}
ncbi_gtdb_mapping_all <- readRDS(GTDB_NCBI_MAPPING_FILE)

# Filter to GTDB 214 and 95 (SMAG uses both versions)
mapping_key <- ncbi_gtdb_mapping_all[ncbi_gtdb_mapping_all$GTDB_version == "GTDB214", ]
mapping_key_gtdb95 <- ncbi_gtdb_mapping_all[ncbi_gtdb_mapping_all$GTDB_version == "GTDB95", ]

log_message(paste("Loaded GTDB 214 mappings:", nrow(mapping_key), "entries"), source_name = source_name)
log_message(paste("Loaded GTDB 95 mappings:", nrow(mapping_key_gtdb95), "entries"), source_name = source_name)

# Check for GTDB 214 metadata (used for better taxonomy matching)
gtdb_214_metadata <- NULL
if (exists("gtdb_214_metadata", envir = .GlobalEnv)) {
  gtdb_214_metadata <- get("gtdb_214_metadata", envir = .GlobalEnv)
  log_message("Using GTDB 214 metadata from global environment", source_name = source_name)
} else {
  log_message("GTDB 214 metadata not found - will use mapping keys only", source_name = source_name)
}

# Prepare mapping for higher classification fallback
# The mapping already contains both NCBI 2021 and 2023 versions
# We'll use it for higher rank fallback matching

# Read in already downloaded genomes (check both local and remote)
log_message("Scanning for downloaded genomes", source_name = source_name)
smag_downloaded <- data.frame()

# Check local directory
local_genomes <- read_in_genomes(smag_config$local_genome_dir, pattern = smag_config$genome_pattern)
if (nrow(local_genomes) > 0) {
  smag_downloaded <- rbind(smag_downloaded, local_genomes)
  log_message(paste("Found", nrow(local_genomes), "genomes in local directory"), source_name = source_name)
}

# Check genome directory (if different from local)
if ("genome_dir" %in% names(smag_config) && smag_config$genome_dir != smag_config$local_genome_dir) {
  remote_genomes <- read_in_genomes(smag_config$genome_dir, pattern = smag_config$genome_pattern)
  if (nrow(remote_genomes) > 0) {
    smag_downloaded <- rbind(smag_downloaded, remote_genomes)
    log_message(paste("Found", nrow(remote_genomes), "genomes in genome directory"), source_name = source_name)
  }
}

# Remove duplicates if genomes exist in both local and remote
if (nrow(smag_downloaded) > 0) {
  smag_downloaded <- smag_downloaded %>%
    distinct(filename, .keep_all = TRUE)
  log_message(paste("Total unique genomes found:", nrow(smag_downloaded)), source_name = source_name)
}

# Read SMAG metadata
log_message("Reading SMAG metadata", source_name = source_name)
smag_metadata <- readxl::read_xlsx(smag_config$metadata_file, skip = 1)

# Filter high-quality MAGs
log_message(paste("Filtering for high-quality MAGs (completeness >", MIN_COMPLETENESS, "%, contamination <", MAX_CONTAMINATION, "%)"), source_name = source_name)
smag_hq <- smag_metadata %>%
  filter(completeness > MIN_COMPLETENESS & contamination < MAX_CONTAMINATION)

# Apply test mode limit if specified
if (is_test_mode() && !is.na(MAX_GENOMES)) {
  log_message(paste("TEST MODE: Limiting to", MAX_GENOMES, "genomes"), source_name = source_name)
  smag_hq <- smag_hq %>% slice_head(n = MAX_GENOMES)
}

log_message(paste("Found", nrow(smag_hq), "high-quality SMAG genomes"), source_name = source_name)

# Find missing genomes
smag_hq_missing <- smag_hq %>%
  filter(!user_genome %in% smag_downloaded$user_genome)

if (nrow(smag_hq_missing) > 0) {
  log_message(paste("Downloading", nrow(smag_hq_missing), "missing high-quality genomes"), source_name = source_name)
  
  # In test mode, limit downloads to MAX_GENOMES
  if (is_test_mode() && !is.na(MAX_GENOMES)) {
    n_to_download <- min(nrow(smag_hq_missing), MAX_GENOMES - nrow(smag_downloaded))
    if (n_to_download > 0) {
      smag_hq_missing <- smag_hq_missing[1:n_to_download,]
      log_message(paste("TEST MODE: Downloading", nrow(smag_hq_missing), "genomes (limited by MAX_GENOMES)"), source_name = source_name)
    } else {
      log_message("TEST MODE: Already have MAX_GENOMES genomes, skipping downloads", source_name = source_name)
      smag_hq_missing <- smag_hq_missing[0,]
    }
  }
  
  if (nrow(smag_hq_missing) > 0) {
    for (i in seq_len(nrow(smag_hq_missing))) {
      genome_id <- smag_hq_missing$user_genome[i]
      dl_link <- paste0(smag_config$download_base_url, genome_id, ".fa")
      # Use local directory for downloads
      file_path <- file.path(smag_config$local_genome_dir, paste0(genome_id, ".fa"))
      
      if (!file.exists(file_path)) {
        tryCatch({
          download.file(dl_link, destfile = file_path, quiet = TRUE)
          log_message(paste("Downloaded:", genome_id), source_name = source_name)
        }, error = function(e) {
          log_message(paste("Failed to download", genome_id, ":", e$message), source_name = source_name)
        })
      }
    }
  }
} else {
  log_message("All high-quality genomes already downloaded", source_name = source_name)
}

# Re-scan downloaded genomes (check both local and remote)
smag_downloaded <- data.frame()
local_genomes <- read_in_genomes(smag_config$local_genome_dir, pattern = smag_config$genome_pattern)
if (nrow(local_genomes) > 0) {
  smag_downloaded <- rbind(smag_downloaded, local_genomes)
}
if ("genome_dir" %in% names(smag_config) && smag_config$genome_dir != smag_config$local_genome_dir) {
  remote_genomes <- read_in_genomes(smag_config$genome_dir, pattern = smag_config$genome_pattern)
  if (nrow(remote_genomes) > 0) {
    smag_downloaded <- rbind(smag_downloaded, remote_genomes)
  }
}
if (nrow(smag_downloaded) > 0) {
  smag_downloaded <- smag_downloaded %>%
    distinct(filename, .keep_all = TRUE)
}

# Merge metadata with downloaded files
# In test mode, process genomes even if files aren't downloaded yet
if (nrow(smag_downloaded) > 0) {
  smag_hq_downloaded <- merge(smag_hq, smag_downloaded, by.x = "user_genome", by.y = "user_genome", all = FALSE)
  if (nrow(smag_hq_downloaded) == 0) {
    stop("No high-quality genomes found after merging with downloaded files")
  }
} else {
  # No downloaded files - use metadata only (for test mode or when files will be downloaded separately)
  if (is_test_mode()) {
    log_message("TEST MODE: Processing metadata without downloaded genome files", source_name = source_name)
    smag_hq_downloaded <- smag_hq
    # Create placeholder filepath column
    smag_hq_downloaded$filepath <- file.path(smag_config$local_genome_dir, paste0(smag_hq_downloaded$user_genome, smag_config$genome_pattern))
  } else {
    stop("No downloaded genomes found. Please download genomes first or run in TEST_MODE.")
  }
}

log_message(paste("Processing", nrow(smag_hq_downloaded), "genomes"), source_name = source_name)

# Parse GTDB taxonomy
smag_hq_downloaded <- smag_hq_downloaded %>%
  separate(GTDB214, sep = ";", into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), remove = FALSE) %>%
  mutate(rank = ifelse(species != "s__", "species",
                       ifelse(genus != "g__", "genus",
                              ifelse(family != "f__", "family",
                                     ifelse(order != "o__", "order",
                                            ifelse(class != "c__", "class",
                                                   ifelse(phylum != "p__", "phylum",
                                                          ifelse(kingdom != "d__", "kingdom", NA))))))))

# Determine GTDB taxon based on rank (strip prefixes like d__, p__, c__, etc.)
smag_hq_downloaded <- smag_hq_downloaded %>%
  mutate(GTDB_taxon = ifelse(rank == "phylum", gsub("^p__", "", phylum),
                             ifelse(rank == "class", gsub("^c__", "", class),
                                    ifelse(rank == "order", gsub("^o__", "", order),
                                           ifelse(rank == "family", gsub("^f__", "", family),
                                                  ifelse(rank == "genus", gsub("^g__", "", genus),
                                                         ifelse(rank == "species", gsub("^s__", "", species), NA)))))))

# Determine higher classification for fallback taxonomy matching (strip prefixes)
smag_hq_downloaded$higher_classification <- ifelse(smag_hq_downloaded$rank == "species", gsub("^g__", "", smag_hq_downloaded$genus),
                                                    ifelse(smag_hq_downloaded$rank == "genus", gsub("^f__", "", smag_hq_downloaded$family),
                                                           ifelse(smag_hq_downloaded$rank == "family", gsub("^o__", "", smag_hq_downloaded$order),
                                                                  ifelse(smag_hq_downloaded$rank == "order", gsub("^c__", "", smag_hq_downloaded$class),
                                                                         ifelse(smag_hq_downloaded$rank == "class", gsub("^p__", "", smag_hq_downloaded$phylum), NA)))))

# Map to NCBI taxids using multiple strategies
log_message("Mapping GTDB taxonomy to NCBI taxids", source_name = source_name)

# Strategy 1: Direct mapping from GTDB 214 metadata (if available)
if (!is.null(gtdb_214_metadata)) {
  smag_hq_downloaded$NCBI_TaxID_genome_214 <- gtdb_214_metadata[match(smag_hq_downloaded$GTDB214, gtdb_214_metadata$gtdb_taxonomy),]$ncbi_taxid
} else {
  smag_hq_downloaded$NCBI_TaxID_genome_214 <- NA
}

# Strategy 2: Mapping from GTDB 214 key
smag_hq_downloaded$NCBI_TaxID_gtdb214 <- mapping_key[match(smag_hq_downloaded$GTDB_taxon, mapping_key$GTDB_taxon),]$ncbi_tax_id

# Strategy 3: Mapping from GTDB 95 key
smag_hq_downloaded$NCBI_TaxID_gtdb95 <- mapping_key_gtdb95[match(smag_hq_downloaded$GTDB_taxon, mapping_key_gtdb95$GTDB_taxon),]$ncbi_tax_id

# Strategy 4: Higher classification mapping using mapping file
# Try matching at higher rank (e.g., if species fails, try genus; if genus fails, try family, etc.)
# Use both NCBI versions from mapping (prioritize 2023, then 2021)
mapping_2023 <- ncbi_gtdb_mapping_all[ncbi_gtdb_mapping_all$NCBI_version == "NCBI_2023-12-01" & 
                                       ncbi_gtdb_mapping_all$GTDB_version %in% c("GTDB214", "GTDB95"), ]
mapping_2021 <- ncbi_gtdb_mapping_all[ncbi_gtdb_mapping_all$NCBI_version == "NCBI_2021-01-01" & 
                                       ncbi_gtdb_mapping_all$GTDB_version %in% c("GTDB214", "GTDB95"), ]

smag_hq_downloaded$higher_NCBI_TaxID_2023 <- mapping_2023[match(smag_hq_downloaded$higher_classification, mapping_2023$GTDB_taxon),]$ncbi_tax_id
smag_hq_downloaded$higher_NCBI_TaxID_2021 <- mapping_2021[match(smag_hq_downloaded$higher_classification, mapping_2021$GTDB_taxon),]$ncbi_tax_id

# Combine all mapping strategies (prioritize in order)
smag_hq_downloaded$NCBI_TaxID <- ifelse(!is.na(smag_hq_downloaded$NCBI_TaxID_genome_214), smag_hq_downloaded$NCBI_TaxID_genome_214,
                                         ifelse(!is.na(smag_hq_downloaded$NCBI_TaxID_gtdb214), smag_hq_downloaded$NCBI_TaxID_gtdb214,
                                                ifelse(!is.na(smag_hq_downloaded$NCBI_TaxID_gtdb95), smag_hq_downloaded$NCBI_TaxID_gtdb95,
                                                       ifelse(!is.na(smag_hq_downloaded$higher_NCBI_TaxID_2023), smag_hq_downloaded$higher_NCBI_TaxID_2023,
                                                              ifelse(!is.na(smag_hq_downloaded$higher_NCBI_TaxID_2021), smag_hq_downloaded$higher_NCBI_TaxID_2021,
                                                                     NA)))))

# Report mapping success
mapped_count <- sum(!is.na(smag_hq_downloaded$NCBI_TaxID))
log_message(paste("Successfully mapped", mapped_count, "of", nrow(smag_hq_downloaded), "genomes to NCBI taxids"), source_name = source_name)

# Filter to genomes with NCBI taxids
ready_genomes_smag <- smag_hq_downloaded %>%
  filter(!is.na(NCBI_TaxID)) %>%
  mutate(source = "SMAG",
         is_published = "Y",
         kingdom = gsub("d__", "", kingdom),
         phylum = gsub("p__", "", phylum),
         species = gsub("s__", "", species)) %>%
  select(ncbi_species_taxid = NCBI_TaxID,
         accession = user_genome,
         ncbi_organism_name = GTDB_taxon,
         kingdom,
         phylum,
         species,
         fasta_file_path = filepath,
         ncbi_taxonomy = GTDB214,
         source,
         is_published)

log_message(paste("Writing", nrow(ready_genomes_smag), "genomes to output file"), source_name = source_name)
write_tsv(ready_genomes_smag, smag_config$output_file)

log_message(paste("SMAG processing complete. Output written to:", smag_config$output_file), source_name = source_name)
