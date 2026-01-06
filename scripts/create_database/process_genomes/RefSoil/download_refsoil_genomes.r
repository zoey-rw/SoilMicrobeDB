# Process RefSoil genomes for Soil Microbe Database
# RefSoil is a curated list of soil organisms from JGI GOLD
# Downloads genomes, processes metadata, and creates Struo2 input file
#
# Usage:
#   Rscript download_refsoil_genomes.r
#   TEST_MODE=TRUE MAX_GENOMES=10 Rscript download_refsoil_genomes.r  # Test with 10 genomes

library(readxl)
library(tidyverse)
library(data.table)
library(jsonlite)
library(CHNOSZ)

# Load centralized configuration
# Find the create_database directory (go up from process_genomes/RefSoil)
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
# Go up to create_database directory (from process_genomes/RefSoil)
config_path <- file.path(script_dir, "../..", "config.R")
if (!file.exists(config_path)) {
  # Try alternative path from project root
  config_path <- file.path(script_dir, "..", "..", "..", "..", "scripts", "create_database", "config.R")
}
if (!file.exists(config_path)) {
  stop("Cannot find config.R. Please ensure you're running from the project directory.")
}
source(config_path)

# Get RefSoil-specific configuration
refsoil_config <- get_source_config("RefSoil")
source_name <- "RefSoil"

# Load helper functions
if (!file.exists(HELPER_FUNCTIONS)) {
  stop(paste("Helper functions not found at:", HELPER_FUNCTIONS))
}
source(HELPER_FUNCTIONS)

log_message("Starting RefSoil processing", source_name = source_name)

# Get RefSoil file paths from config
REFSOIL_DATA_FILE_CONFIG <- Sys.getenv("REFSOIL_DATA_FILE", unset = refsoil_config$data_file)

# Create local directory if it doesn't exist
REFSOIL_LOCAL_DIR <- get_source_dir("RefSoil")
REFSOIL_REF_DATA_DIR <- file.path(REFSOIL_LOCAL_DIR, "ref_data")
dir.create(REFSOIL_REF_DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# Local file path
local_refsoil_data <- file.path(REFSOIL_REF_DATA_DIR, "RefSoil_s1.xlsx")

# Get RefSoil data file
# Use config file path if it exists, otherwise try local path
if (file.exists(REFSOIL_DATA_FILE_CONFIG)) {
  REFSOIL_DATA_FILE <- REFSOIL_DATA_FILE_CONFIG
} else if (file.exists(local_refsoil_data)) {
  REFSOIL_DATA_FILE <- local_refsoil_data
} else {
  stop(paste("RefSoil data file not found. Tried:",
             "\n  Config path:", REFSOIL_DATA_FILE_CONFIG,
             "\n  Local path:", local_refsoil_data,
             "\nPlease ensure file is available locally or mount remote server."))
}

# Read RefSoil data file
log_message("Reading RefSoil data file", source_name = source_name)
refsoil <- readxl::read_xlsx(REFSOIL_DATA_FILE)
log_message(paste("Found", nrow(refsoil), "RefSoil organisms"), source_name = source_name)

# Get JGI GOLD data files (RefSoil uses JGI GOLD to find matching genomes)
# Use the same approach as JGI_Gold script
JGI_GOLD_DATA_FILE_CONFIG <- Sys.getenv("JGI_GOLD_DATA_FILE", unset = NA)
if (is.na(JGI_GOLD_DATA_FILE_CONFIG)) {
  JGI_GOLD_LOCAL_DIR <- get_source_dir("JGI_Gold")
  local_gold_data <- file.path(JGI_GOLD_LOCAL_DIR, "ref_data", "goldData.xlsx")
  remote_gold_data <- file.path(Sys.getenv("JGI_GOLD_REMOTE_BASE", unset = "/projectnb/frpmars/soil_microbe_db/scripts/database_development/genome_downloads/jgi_gold/ref_data"), "goldData.xlsx")
  JGI_GOLD_DATA_FILE_CONFIG <- if (file.exists(local_gold_data)) local_gold_data else remote_gold_data
}

# Get GOLD data file
JGI_GOLD_LOCAL_DIR <- get_source_dir("JGI_Gold")
JGI_GOLD_REF_DATA_DIR <- file.path(JGI_GOLD_LOCAL_DIR, "ref_data")
dir.create(JGI_GOLD_REF_DATA_DIR, showWarnings = FALSE, recursive = TRUE)
local_gold_data <- file.path(JGI_GOLD_REF_DATA_DIR, "goldData.xlsx")

# Use config file path if it exists, otherwise try local path
if (file.exists(JGI_GOLD_DATA_FILE_CONFIG)) {
  GOLD_DATA_FILE <- JGI_GOLD_DATA_FILE_CONFIG
} else if (file.exists(local_gold_data)) {
  GOLD_DATA_FILE <- local_gold_data
} else {
  stop(paste("GOLD data file not found. Tried:",
             "\n  Config path:", JGI_GOLD_DATA_FILE_CONFIG,
             "\n  Local path:", local_gold_data,
             "\nPlease ensure file is available locally or mount remote server."))
}

# Read GOLD data files (only sheets we need)
log_message("Reading GOLD data files", source_name = source_name)
gold_data_organism <- readxl::read_xlsx(GOLD_DATA_FILE, sheet = 4)
gold_data_sequencing <- readxl::read_xlsx(GOLD_DATA_FILE, sheet = 5)
gold_data_analysis <- readxl::read_xlsx(GOLD_DATA_FILE, sheet = 6)

# Filter GOLD organisms to those matching RefSoil tax IDs
log_message("Filtering GOLD organisms to RefSoil tax IDs", source_name = source_name)
refsoil_organisms <- gold_data_organism %>% 
  filter(`ORGANISM NCBI TAX ID` %in% refsoil$`Taxon ID`)
log_message(paste("Found", nrow(refsoil_organisms), "GOLD organisms matching RefSoil tax IDs"), source_name = source_name)

refsoil_organism_ids <- refsoil_organisms %>% select(`ORGANISM GOLD ID`) %>% unique() %>% unlist()

# Filter to genomes (analysis projects) from RefSoil organisms
refsoil_sequencing <- gold_data_sequencing %>% filter(`ORGANISM GOLD ID` %in% refsoil_organism_ids)
refsoil_sequencing_ids <- refsoil_sequencing %>% select(`PROJECT GOLD ID`) %>% unlist()

# Filter to genomes (sequencing metadata) from RefSoil organisms
refsoil_ap <- gold_data_analysis %>% filter(`AP PROJECT GOLD IDS` %in% refsoil_sequencing_ids)

# Merge sequencing and analysis data
genome_info <- merge(refsoil_ap, refsoil_sequencing, by.x = "AP PROJECT GOLD IDS", by.y = "PROJECT GOLD ID")
log_message(paste("Found", nrow(genome_info), "genomes from RefSoil organisms"), source_name = source_name)

# Extract NCBI accessions from JSON
log_message("Extracting NCBI accessions from JSON", source_name = source_name)
get_accession <- function(accession_info_json) {
  if (is.na(accession_info_json) || accession_info_json == "") {
    return(NA)
  }
  tryCatch({
    accession_info <- jsonlite::parse_json(accession_info_json)
    if (length(accession_info) == 0) {
      return(NA)
    }
    out <- accession_info[[1]]
    to_return <- ifelse(is.null(out$assemblyAccession), NA, out$assemblyAccession)
    return(to_return)
  }, error = function(e) {
    return(NA)
  })
}

genome_info$NCBI_accession <- NA
for (i in seq_len(nrow(genome_info))) {
  if (i %% 1000 == 0) {
    log_message(paste("Processing accession", i, "of", nrow(genome_info)), source_name = source_name)
  }
  genome_info$NCBI_accession[i] <- get_accession(genome_info$`AP GENBANK`[i])
}

# Remove genomes without NCBI accessions
genome_info <- genome_info %>% filter(!is.na(NCBI_accession))
log_message(paste("Found", nrow(genome_info), "genomes with NCBI accessions"), source_name = source_name)

# Apply test mode limit if specified
if (is_test_mode() && !is.na(MAX_GENOMES)) {
  log_message(paste("TEST MODE: Limiting to", MAX_GENOMES, "genomes"), source_name = source_name)
  genome_info <- genome_info %>% slice_head(n = MAX_GENOMES)
}

# Ensure accessions have version numbers
genome_info <- genome_info %>%
  mutate(NCBI_accession = ifelse(!grepl("\\.\\d+$", NCBI_accession), 
                                  paste0(NCBI_accession, ".1"), 
                                  NCBI_accession))

# Load GTDB-NCBI mapping (contains all GTDB versions)
log_message("Loading GTDB-NCBI mapping", source_name = source_name)
if (!file.exists(GTDB_NCBI_MAPPING_FILE)) {
  stop(paste("GTDB-NCBI mapping file not found:", GTDB_NCBI_MAPPING_FILE,
             "\nPlease ensure remote server is mounted or copy file locally."))
}
ncbi_gtdb_mapping_all <- readRDS(GTDB_NCBI_MAPPING_FILE)

# Filter to GTDB 207 for RefSoil (can be adjusted if needed)
gtdb_ncbi_mapping <- ncbi_gtdb_mapping_all[ncbi_gtdb_mapping_all$GTDB_version == "GTDB207", ]
log_message(paste("Filtered to GTDB 207 mappings:", nrow(gtdb_ncbi_mapping), "entries"), source_name = source_name)

# Load GTDB 207 metadata (required for taxonomy matching)
if (!exists("gtdb_207_metadata", envir = .GlobalEnv) || is.null(get("gtdb_207_metadata", envir = .GlobalEnv))) {
  gtdb_207_bac_file <- file.path(GENOME_DB_DIR, "bac120_metadata_r207.tsv")
  gtdb_207_ar_file <- file.path(GENOME_DB_DIR, "ar53_metadata_r207.tsv")
  
  if (file.exists(gtdb_207_bac_file) && file.exists(gtdb_207_ar_file)) {
    log_message("Loading GTDB 207 metadata from local files", source_name = source_name)
    gtdb_207_metadata_bac <- fread(gtdb_207_bac_file) %>%
      select(accession, checkm_completeness, checkm_contamination, gtdb_genome_representative, gtdb_taxonomy,
             ncbi_date, ncbi_genbank_assembly_accession, ncbi_organism_name, ncbi_taxid, ncbi_taxonomy, ncbi_genome_category) %>%
      mutate(source = "GTDB_207")
    gtdb_207_metadata_ar <- fread(gtdb_207_ar_file) %>%
      select(accession, checkm_completeness, checkm_contamination, gtdb_genome_representative, gtdb_taxonomy,
             ncbi_date, ncbi_genbank_assembly_accession, ncbi_organism_name, ncbi_taxid, ncbi_taxonomy, ncbi_genome_category) %>%
      mutate(source = "GTDB_207")
    gtdb_207_metadata <- rbind(gtdb_207_metadata_bac, gtdb_207_metadata_ar) %>%
      mutate(is_MAG = ifelse(ncbi_genome_category == "derived from metagenome", TRUE, FALSE)) %>%
      select(-ncbi_genome_category)
    assign("gtdb_207_metadata", gtdb_207_metadata, envir = .GlobalEnv)
  } else {
    stop(paste("GTDB 207 metadata files are required but not found. Tried:",
               "\n  Local:", gtdb_207_bac_file, "and", gtdb_207_ar_file,
               "\nPlease ensure files are available locally or mount remote server."))
  }
} else {
  gtdb_207_metadata <- get("gtdb_207_metadata", envir = .GlobalEnv)
  if (is.null(gtdb_207_metadata) || nrow(gtdb_207_metadata) == 0) {
    stop("GTDB 207 metadata is required but is NULL or empty.")
  }
  log_message(paste("Using GTDB 207 metadata from global environment (", nrow(gtdb_207_metadata), "genomes)"), source_name = source_name)
}

# Map NCBI accessions to NCBI taxids using GTDB 207 metadata
log_message("Mapping NCBI accessions to NCBI taxids", source_name = source_name)
genome_info$base_accession <- gsub("\\.\\d+$", "", genome_info$NCBI_accession)
gtdb_207_metadata$base_accession <- gsub("\\.\\d+$", "", gtdb_207_metadata$ncbi_genbank_assembly_accession)

genome_info <- genome_info %>%
  left_join(gtdb_207_metadata %>% 
              select(base_accession, ncbi_taxid, ncbi_organism_name, gtdb_taxonomy),
            by = "base_accession")

# Filter to genomes with NCBI taxids
genome_info <- genome_info %>% filter(!is.na(ncbi_taxid))
log_message(paste("Found", nrow(genome_info), "genomes with NCBI taxids"), source_name = source_name)

# Check for already downloaded genomes
# RefSoil genomes are typically downloaded from NCBI, so check NCBI genome directory
NCBI_GENOME_DIR <- refsoil_config$genome_dir
NCBI_GENOME_LOCAL_DIR <- refsoil_config$local_genome_dir

# Use genome_dir as the primary directory (may be same as local)
if (is.null(NCBI_GENOME_DIR) || is.na(NCBI_GENOME_DIR)) {
  NCBI_GENOME_DIR <- NCBI_GENOME_LOCAL_DIR
}

# Check both local and remote directories
already_downloaded <- data.frame()
if (dir.exists(NCBI_GENOME_LOCAL_DIR)) {
  already_downloaded <- read_in_genomes(NCBI_GENOME_LOCAL_DIR, pattern = "\\.fna\\.gz$")
}

# Also check genome directory if different from local
if (NCBI_GENOME_DIR != NCBI_GENOME_LOCAL_DIR && dir.exists(NCBI_GENOME_DIR)) {
  remote_downloaded <- read_in_genomes(NCBI_GENOME_DIR, pattern = "\\.fna\\.gz$")
  already_downloaded <- rbind(already_downloaded, remote_downloaded)
}

if (nrow(already_downloaded) > 0) {
  already_downloaded <- already_downloaded %>%
    mutate(NCBI_accession = sapply(filename, function(x) {
      x1 <- strsplit(x, "_") %>% unlist
      return(paste0(x1[[1]], "_", x1[[2]]))
    })) %>%
    mutate(base_accession = gsub("\\.\\d+$", "", NCBI_accession)) %>%
    distinct(base_accession, .keep_all = TRUE)
  
  genome_info$already_downloaded <- ifelse(genome_info$base_accession %in% already_downloaded$base_accession, TRUE, FALSE)
  log_message(paste("Found", sum(genome_info$already_downloaded), "genomes already downloaded"), source_name = source_name)
} else {
  genome_info$already_downloaded <- FALSE
}

# Download missing genomes from NCBI
to_download <- genome_info %>% filter(!already_downloaded)

if (nrow(to_download) > 0) {
  # In test mode, limit downloads to MAX_GENOMES
  if (is_test_mode() && !is.na(MAX_GENOMES)) {
    n_already <- sum(genome_info$already_downloaded)
    n_to_download <- min(nrow(to_download), MAX_GENOMES - n_already)
    if (n_to_download > 0) {
      to_download <- to_download[1:n_to_download,]
      log_message(paste("TEST MODE: Downloading", nrow(to_download), "genomes (limited by MAX_GENOMES)"), source_name = source_name)
    } else {
      log_message("TEST MODE: Already have MAX_GENOMES genomes, skipping downloads", source_name = source_name)
      to_download <- to_download[0,]
    }
  } else {
    log_message(paste("Downloading", nrow(to_download), "missing genomes from NCBI"), source_name = source_name)
  }
  
  if (nrow(to_download) > 0) {
    # Create local directory if it doesn't exist
    dir.create(NCBI_GENOME_LOCAL_DIR, showWarnings = FALSE, recursive = TRUE)
    
    # Construct NCBI FTP paths from accessions
    # Pattern: GCA_000123456.1 -> ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/123/456/GCA_000123456.1/GCA_000123456.1_genomic.fna.gz
    construct_ncbi_ftp_path <- function(accession) {
      base_acc <- gsub("\\.\\d+$", "", accession)
      # Extract the numeric part (e.g., "000123456" from "GCA_000123456")
      numeric_part <- gsub("^[A-Z]+_", "", base_acc)
      # Split into groups of 3 digits for directory structure
      dir_parts <- c(substr(numeric_part, 1, 3),
                     substr(numeric_part, 4, 6),
                     substr(numeric_part, 7, 9))
      dir_parts <- dir_parts[dir_parts != ""]
      ftp_path <- paste0("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/",
                         gsub("_.*", "", base_acc), "/",
                         paste(dir_parts, collapse = "/"), "/",
                         accession, "/",
                         accession, "_genomic.fna.gz")
      return(ftp_path)
    }
    
    to_download$ftp_path <- sapply(to_download$NCBI_accession, construct_ncbi_ftp_path)
    to_download$local_path <- file.path(NCBI_GENOME_LOCAL_DIR, paste0(gsub("\\.", "_", to_download$base_accession), "_genomic.fna.gz"))
    
    for (i in seq_len(nrow(to_download))) {
      local_path <- to_download$local_path[i]
      if (!file.exists(local_path)) {
        tryCatch({
          download.file(to_download$ftp_path[i], destfile = local_path, quiet = TRUE, mode = "wb")
          log_message(paste("Downloaded:", basename(local_path)), source_name = source_name)
        }, error = function(e) {
          log_message(paste("Failed to download", to_download$NCBI_accession[i], ":", e$message), source_name = source_name)
        })
      }
    }
  }
} else {
  log_message("All selected genomes already downloaded", source_name = source_name)
}

# Prepare output with NCBI taxonomy information
log_message("Adding NCBI taxonomy information", source_name = source_name)

# Use local NCBI taxonomy directory (files should be copied via copy_ncbi_taxonomy_local.sh)
if (!file.exists(file.path(NCBI_TAX_DIR, "nodes.dmp"))) {
  stop(paste("NCBI taxonomy directory not found: ", NCBI_TAX_DIR,
             "\nPlease copy NCBI taxonomy files using: bash scripts/create_database/copy_ncbi_taxonomy_local.sh"))
}

log_message(paste("Using NCBI taxonomy directory:", NCBI_TAX_DIR), source_name = source_name)

# Update file paths to check if files actually exist (some downloads may have failed)
genome_info <- genome_info %>%
  mutate(expected_file_path = file.path(NCBI_GENOME_LOCAL_DIR, paste0(gsub("\\.", "_", base_accession), "_genomic.fna.gz"))) %>%
  mutate(file_exists = file.exists(expected_file_path))

# For test mode, we'll include genomes even if download failed (they may exist remotely)
# For production, you might want to filter to only genomes with existing files
log_message(paste("Genomes with local files:", sum(genome_info$file_exists), "of", nrow(genome_info)), source_name = source_name)

nodes <- getnodes(taxdir = NCBI_TAX_DIR)
genome_info$rank <- CHNOSZ::getrank(genome_info$ncbi_taxid, NCBI_TAX_DIR, nodes = nodes)

# Optimized taxonomy reading: use grep for small numbers of tax IDs
unique_tax_ids <- unique(genome_info$ncbi_taxid)
n_tax_ids <- length(unique_tax_ids)

log_message(paste("Reading NCBI taxonomy for", n_tax_ids, "unique tax IDs"), source_name = source_name)

# For small numbers of tax IDs (test mode), use grep to extract only needed lines
# This is much faster than reading the entire 315MB file
if (n_tax_ids <= 100 || is_test_mode()) {
  log_message("Using optimized grep-based extraction (fast for small tax ID sets)", source_name = source_name)
  
  ranked_lineage_file <- file.path(NCBI_TAX_DIR, "rankedlineage.dmp")
  
  if (!file.exists(ranked_lineage_file)) {
    stop(paste("rankedlineage.dmp not found at:", ranked_lineage_file))
  }
  
  # Use awk to match first field exactly (tax_id is first field in rankedlineage.dmp)
  # This is much faster than reading the entire 315MB file
  log_message("Extracting relevant lines from taxonomy file using awk...", source_name = source_name)
  
  # Create a temporary file with tax IDs (one per line) for awk to read
  temp_taxid_file <- tempfile()
  writeLines(as.character(unique_tax_ids), temp_taxid_file)
  
  # Use awk to match first field exactly against the list of tax IDs
  # Read tax IDs from file and match first field ($1) exactly
  awk_cmd <- paste0("awk -F'\t' 'FNR==NR {ids[$1]; next} $1 in ids' ", 
                    shQuote(temp_taxid_file), " ", shQuote(ranked_lineage_file))
  tax_lines <- system(awk_cmd, intern = TRUE)
  
  # Clean up temp file
  unlink(temp_taxid_file)
  
  if (length(tax_lines) == 0) {
    # Fallback: try simpler grep approach (matches anywhere in line, less precise)
    log_message("Awk extraction returned no results, trying grep fallback...", source_name = source_name)
    temp_taxid_file <- tempfile()
    writeLines(as.character(unique_tax_ids), temp_taxid_file)
    grep_cmd <- paste0("grep -F -f ", shQuote(temp_taxid_file), " ", shQuote(ranked_lineage_file))
    tax_lines <- system(grep_cmd, intern = TRUE)
    unlink(temp_taxid_file)
  }
  
  if (length(tax_lines) == 0) {
    stop("No matching taxonomy entries found. Check if tax IDs are valid.")
  }
  
  log_message(paste("Found", length(tax_lines), "matching taxonomy entries"), source_name = source_name)
  
  # Parse the extracted lines
  # rankedlineage.dmp format: tax_id\t|\ttax_name\t|\tspecies\t|\t... (tab-pipe-tab delimited)
  # Use data.table::fread for efficient parsing, but only on the small subset
  temp_tax_file <- tempfile()
  writeLines(tax_lines, temp_tax_file)
  
  # Read the subset using fread (fast even with pipe delimiters)
  ncbi_taxonomy <- fread(temp_tax_file, sep = "\t", data.table = FALSE, 
                         col.names = c("tax_id", "tax_name", "species", "genus", "family", 
                                      "order", "class", "phylum", "kingdom", "superkingdom", "extra"),
                         fill = TRUE) %>%
    as_tibble() %>%
    select(-extra) %>%
    # Clean up fields (remove leading/trailing whitespace and empty strings)
    mutate(across(everything(), ~trimws(.x))) %>%
    mutate(across(c(species, genus, family, order, class, phylum, kingdom, superkingdom), 
                  ~ifelse(.x == "" | is.na(.x), NA_character_, .x))) %>%
    mutate(tax_id = as.numeric(tax_id))
  
  unlink(temp_tax_file)
  
  # Join with requested tax IDs (some may not be found)
  ncbi_tax_ids_df <- data.frame(query_tax_id = as.numeric(unique_tax_ids))
  ncbi_taxonomy <- ncbi_taxonomy %>%
    right_join(ncbi_tax_ids_df, by = c("tax_id" = "query_tax_id")) %>%
    arrange(tax_id)
  
  log_message("NCBI taxonomy extraction complete", source_name = source_name)
} else {
  # For large numbers, use the standard function (reads entire file but handles all cases)
  log_message("Using standard taxonomy reading (for large tax ID sets)", source_name = source_name)
  ncbi_taxonomy <- tryCatch({
    tax_id_to_ranked_lineage(genome_info$ncbi_taxid, NCBI_TAX_DIR) %>%
      arrange(tax_id)
  }, error = function(e) {
    log_message(paste("Error reading taxonomy file:", e$message), source_name = source_name)
    stop(paste("Failed to read NCBI taxonomy file. Error:", e$message))
  })
  log_message("NCBI taxonomy read complete", source_name = source_name)
}

genome_info <- genome_info %>% arrange(ncbi_taxid) %>% mutate(ncbi_taxid = as.numeric(ncbi_taxid))

if (identical(genome_info$ncbi_taxid, ncbi_taxonomy$tax_id)) {
  genome_info <- cbind.data.frame(genome_info, ncbi_taxonomy) %>%
    mutate(species = ifelse(rank == "species", tax_name, species)) %>%
    mutate(phylum = ifelse(rank == "phylum", tax_name, phylum)) %>%
    mutate(class = ifelse(rank == "class", tax_name, class)) %>%
    mutate(genus = ifelse(rank == "genus", tax_name, genus)) %>%
    mutate(family = ifelse(rank == "family", tax_name, family)) %>%
    mutate(kingdom = ifelse(rank == "kingdom", tax_name, kingdom)) %>%
    mutate(order = ifelse(rank == "order", tax_name, order)) %>%
    mutate(ncbi_taxonomy = paste0("k__", ifelse(is.na(kingdom) | kingdom == "", superkingdom, kingdom),
                                   ";p__", ifelse(is.na(phylum), "", phylum),
                                   ";c__", ifelse(is.na(class), "", class),
                                   ";o__", ifelse(is.na(order), "", order),
                                   ";f__", ifelse(is.na(family), "", family),
                                   ";g__", ifelse(is.na(genus), "", genus),
                                   ";s__", ifelse(is.na(species), "", species)))
  
  # Use NCBI organism name from taxonomy if available
  genome_info$ncbi_organism_name <- ifelse(!is.na(genome_info$tax_name), 
                                           genome_info$tax_name, 
                                           genome_info$ncbi_organism_name)
  log_message("NCBI taxonomy information added", source_name = source_name)
} else {
  stop("Taxonomy ID mismatch - cannot proceed without valid NCBI taxonomy mapping")
}

# Prepare Struo2 input file
ready_genomes_refsoil <- genome_info %>%
  mutate(source = "RefSoil",
         is_published = ifelse(!is.na(`AP RESTRICTION STATUS`) && `AP RESTRICTION STATUS` == "Unrestricted", "Y", "NA"),
         accession = NCBI_accession,
         ncbi_species_taxid = as.numeric(ncbi_taxid),
         fasta_file_path = expected_file_path) %>%
  select(ncbi_species_taxid,
         accession,
         ncbi_organism_name,
         kingdom,
         phylum,
         species,
         fasta_file_path,
         ncbi_taxonomy,
         source,
         is_published) %>%
  distinct(accession, .keep_all = TRUE) %>%
  arrange(ncbi_species_taxid)

log_message(paste("Writing", nrow(ready_genomes_refsoil), "genomes to output file"), source_name = source_name)
write_tsv(ready_genomes_refsoil, refsoil_config$output_file)

log_message(paste("RefSoil processing complete. Output written to:", refsoil_config$output_file), source_name = source_name)