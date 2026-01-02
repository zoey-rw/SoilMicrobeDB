# Process SPIRE (Soil Profile Inferred Resource Economics) MAGs for Soil Microbe Database
# Downloads genomes, processes metadata, and creates Struo2 input file
#
# Usage:
#   Rscript download_spire_MAGs.r
#   TEST_MODE=TRUE MAX_GENOMES=10 Rscript download_spire_MAGs.r  # Test with 10 genomes

library(data.table)
library(tidyverse)
library(CHNOSZ)

# Load centralized configuration
# Find the create_database directory (go up from process_genomes/SPIRE)
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
# Go up to create_database directory (from process_genomes/SPIRE)
config_path <- file.path(script_dir, "../..", "config.R")
if (!file.exists(config_path)) {
  # Try alternative path from project root
  config_path <- file.path(script_dir, "..", "..", "..", "..", "scripts", "create_database", "config.R")
}
if (!file.exists(config_path)) {
  stop("Cannot find config.R. Please ensure you're running from the project directory.")
}
source(config_path)

# Get SPIRE-specific configuration
spire_config <- get_source_config("SPIRE")
source_name <- "SPIRE"

# Load helper functions
if (!file.exists(HELPER_FUNCTIONS)) {
  stop(paste("Helper functions not found at:", HELPER_FUNCTIONS))
}
source(HELPER_FUNCTIONS)

log_message("Starting SPIRE processing", source_name = source_name)

# Check for required files
# Try to access files - if remote, check if accessible via SSH or mount
check_file_access <- function(file_path) {
  if (file.exists(file_path)) {
    return(TRUE)
  }
  # If remote path, try to check via SSH
  if (grepl("^/projectnb", file_path)) {
    # Try to read first few bytes to check accessibility
    result <- tryCatch({
      con <- file(file_path, "rb")
      readBin(con, "raw", n = 1)
      close(con)
      TRUE
    }, error = function(e) FALSE)
    return(result)
  }
  return(FALSE)
}

if (!check_file_access(spire_config$metadata_file)) {
  stop(paste("SPIRE metadata file not accessible:", spire_config$metadata_file,
             "\nPlease ensure remote server is mounted or copy metadata files locally."))
}

if (!check_file_access(spire_config$taxonomy_file)) {
  stop(paste("SPIRE taxonomy file not accessible:", spire_config$taxonomy_file,
             "\nPlease ensure remote server is mounted or copy metadata files locally."))
}

if (!check_file_access(spire_config$microntology_file)) {
  stop(paste("SPIRE microntology file not accessible:", spire_config$microntology_file,
             "\nPlease ensure remote server is mounted or copy metadata files locally."))
}

# Load GTDB-NCBI mapping (file containing GTDB versions 95, 207, 214)
# Filter to GTDB 207 for SPIRE processing
log_message("Loading GTDB-NCBI mapping (will filter to GTDB 207)", source_name = source_name)
if (file.exists(GTDB_NCBI_MAPPING_FILE)) {
  gtdb_ncbi_mapping_all <- readRDS(GTDB_NCBI_MAPPING_FILE)
} else if (grepl("^/projectnb", GTDB_NCBI_MAPPING_FILE)) {
  # Read remote file via SSH
  log_message("Reading GTDB mapping file from remote server via SSH", source_name = source_name)
  remote_path <- GTDB_NCBI_MAPPING_FILE
  ssh_cmd <- paste0("ssh zrwerbin@scc2.bu.edu 'cat ", shQuote(remote_path), "'")
  temp_file <- tempfile(fileext = ".rds")
  system(paste(ssh_cmd, ">", temp_file))
  if (file.exists(temp_file) && file.size(temp_file) > 0) {
    gtdb_ncbi_mapping_all <- readRDS(temp_file)
    unlink(temp_file)
    log_message("Successfully loaded GTDB mapping from remote", source_name = source_name)
  } else {
    stop(paste("Failed to read GTDB-NCBI mapping file from remote:", GTDB_NCBI_MAPPING_FILE))
  }
} else {
  stop(paste("GTDB-NCBI mapping file not found:", GTDB_NCBI_MAPPING_FILE))
}

# Filter to GTDB 207 for SPIRE
gtdb_ncbi_mapping <- gtdb_ncbi_mapping_all[gtdb_ncbi_mapping_all$GTDB_version == "GTDB207", ]
log_message(paste("Filtered to GTDB 207 mappings:", nrow(gtdb_ncbi_mapping), "entries"), source_name = source_name)

# Load GTDB 207 metadata (required for taxonomy matching)
# This should be loaded by helper_functions.r, but verify it exists and load if needed
if (!exists("gtdb_207_metadata", envir = .GlobalEnv) || is.null(get("gtdb_207_metadata", envir = .GlobalEnv))) {
  # Try to load from local files if available
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
    # Try remote paths via SSH
    remote_bac <- file.path(GTDB_REMOTE_BASE, "bac120_metadata_r207.tsv")
    remote_ar <- file.path(GTDB_REMOTE_BASE, "ar53_metadata_r207.tsv")
    log_message("GTDB 207 metadata not found locally, attempting to load from remote via SSH", source_name = source_name)
    
    # Read remote files via SSH
    temp_bac <- tempfile(fileext = ".tsv")
    temp_ar <- tempfile(fileext = ".tsv")
    
    ssh_cmd_bac <- paste0("ssh zrwerbin@scc2.bu.edu 'cat ", shQuote(remote_bac), "'")
    ssh_cmd_ar <- paste0("ssh zrwerbin@scc2.bu.edu 'cat ", shQuote(remote_ar), "'")
    
    system(paste(ssh_cmd_bac, ">", temp_bac))
    system(paste(ssh_cmd_ar, ">", temp_ar))
    
    if (file.exists(temp_bac) && file.size(temp_bac) > 0 && file.exists(temp_ar) && file.size(temp_ar) > 0) {
      log_message("Loading GTDB 207 metadata from remote files", source_name = source_name)
      gtdb_207_metadata_bac <- fread(temp_bac) %>%
        select(accession, checkm_completeness, checkm_contamination, gtdb_genome_representative, gtdb_taxonomy,
               ncbi_date, ncbi_genbank_assembly_accession, ncbi_organism_name, ncbi_taxid, ncbi_taxonomy, ncbi_genome_category) %>%
        mutate(source = "GTDB_207")
      gtdb_207_metadata_ar <- fread(temp_ar) %>%
        select(accession, checkm_completeness, checkm_contamination, gtdb_genome_representative, gtdb_taxonomy,
               ncbi_date, ncbi_genbank_assembly_accession, ncbi_organism_name, ncbi_taxid, ncbi_taxonomy, ncbi_genome_category) %>%
        mutate(source = "GTDB_207")
      gtdb_207_metadata <- rbind(gtdb_207_metadata_bac, gtdb_207_metadata_ar) %>%
        mutate(is_MAG = ifelse(ncbi_genome_category == "derived from metagenome", TRUE, FALSE)) %>%
        select(-ncbi_genome_category)
      assign("gtdb_207_metadata", gtdb_207_metadata, envir = .GlobalEnv)
      unlink(c(temp_bac, temp_ar))
      log_message("Successfully loaded GTDB 207 metadata from remote", source_name = source_name)
    } else {
      unlink(c(temp_bac, temp_ar))
      stop(paste("GTDB 207 metadata files are required but not accessible. Tried:",
                 "\n  Local:", gtdb_207_bac_file, "and", gtdb_207_ar_file,
                 "\n  Remote:", remote_bac, "and", remote_ar,
                 "\nPlease ensure files are available locally or remotely accessible via SSH."))
    }
  }
} else {
  gtdb_207_metadata <- get("gtdb_207_metadata", envir = .GlobalEnv)
  if (is.null(gtdb_207_metadata) || nrow(gtdb_207_metadata) == 0) {
    stop("GTDB 207 metadata is required but is NULL or empty. Please ensure bac120_metadata_r207.tsv and ar53_metadata_r207.tsv are accessible.")
  }
  log_message(paste("Using GTDB 207 metadata from global environment (", nrow(gtdb_207_metadata), "genomes)"), source_name = source_name)
}

# Read SPIRE metadata
log_message("Reading SPIRE metadata", source_name = source_name)
spire_metadata <- fread(spire_config$metadata_file)

# Create download paths
spire_genome_dir <- get_source_dir("SPIRE")
spire_metadata$download_links <- paste0(spire_config$download_base_url, spire_metadata$genome_id)
spire_metadata$download_path <- file.path(spire_genome_dir, paste0(spire_metadata$genome_id, ".fna.gz"))

log_message(paste("Total genomes in metadata:", nrow(spire_metadata)), source_name = source_name)
log_message(paste("Domain distribution:"), source_name = source_name)
print(table(spire_metadata$domain))

# Filter high-quality MAGs
log_message(paste("Filtering for high-quality MAGs (completeness >", MIN_COMPLETENESS, "%, contamination <", MAX_CONTAMINATION, "%)"), source_name = source_name)
spire_hq <- spire_metadata %>%
  filter(completeness > MIN_COMPLETENESS & contamination < MAX_CONTAMINATION)

log_message(paste("Found", nrow(spire_hq), "high-quality SPIRE genomes"), source_name = source_name)

# Read microntology to filter for soil samples
log_message("Reading SPIRE microntology", source_name = source_name)
spire_microntology <- fread(spire_config$microntology_file, header = FALSE, col.names = c("sampleID", "studyID", "accession", "habitat"))

# Filter for soil-related habitats
soil_samples <- spire_microntology %>%
  filter(grepl(paste(spire_config$soil_habitat_keywords, collapse = "|"), habitat, ignore.case = TRUE))

log_message(paste("Found", nrow(soil_samples), "soil-related samples"), source_name = source_name)

# Filter genomes from soil samples
soil_genomes <- spire_hq %>%
  filter(derived_from_sample %in% soil_samples$sampleID)

log_message(paste("Found", nrow(soil_genomes), "high-quality genomes from soil samples"), source_name = source_name)

# Apply test mode limit if specified
if (is_test_mode() && !is.na(MAX_GENOMES)) {
  log_message(paste("TEST MODE: Limiting to", MAX_GENOMES, "genomes"), source_name = source_name)
  soil_genomes <- soil_genomes %>% slice_head(n = MAX_GENOMES)
}

# Check for already downloaded genomes locally
already_downloaded <- list.files(spire_genome_dir, full.names = TRUE, pattern = spire_config$genome_pattern)
log_message(paste("Found", length(already_downloaded), "genomes already downloaded locally"), source_name = source_name)

to_download <- soil_genomes[!basename(soil_genomes$download_path) %in% basename(already_downloaded),]

if (nrow(to_download) > 0) {
  # In test mode, limit downloads to MAX_GENOMES
  if (is_test_mode() && !is.na(MAX_GENOMES)) {
    n_to_download <- min(nrow(to_download), MAX_GENOMES - length(already_downloaded))
    if (n_to_download > 0) {
      to_download <- to_download[1:n_to_download,]
      log_message(paste("TEST MODE: Downloading", nrow(to_download), "genomes (limited by MAX_GENOMES)"), source_name = source_name)
    } else {
      log_message("TEST MODE: Already have MAX_GENOMES genomes, skipping downloads", source_name = source_name)
      to_download <- to_download[0,]
    }
  } else {
    log_message(paste("Downloading", nrow(to_download), "missing genomes"), source_name = source_name)
  }
  
  if (nrow(to_download) > 0) {
    for (i in seq_len(nrow(to_download))) {
      local_path <- to_download$download_path[i]
      if (!file.exists(local_path)) {
        tryCatch({
          download.file(to_download$download_links[i], destfile = local_path, quiet = TRUE)
          log_message(paste("Downloaded:", basename(local_path)), source_name = source_name)
        }, error = function(e) {
          log_message(paste("Failed to download", to_download$genome_id[i], ":", e$message), source_name = source_name)
        })
      }
    }
  }
} else {
  log_message("All selected genomes already downloaded locally", source_name = source_name)
}

# Map to NCBI taxids using GTDB 207 metadata (required)
# gtdb_207_metadata should be loaded above
if (is.null(gtdb_207_metadata) || nrow(gtdb_207_metadata) == 0) {
  stop("GTDB 207 metadata is required but not available. Please ensure bac120_metadata_r207.tsv and ar53_metadata_r207.tsv are accessible.")
}
soil_genomes$simple_NCBI_match <- gtdb_207_metadata[match(soil_genomes$classification, gtdb_207_metadata$gtdb_taxonomy),]$ncbi_taxid

# Clean up classification
soil_genomes <- soil_genomes %>%
  mutate(domain = ifelse(classification == "Unclassified Archaea", "Archaea",
                         ifelse(classification == "Unclassified Bacteria", "Bacteria", domain))) %>%
  mutate(classification = ifelse(classification == "", "Unclassified", classification))

# Determine taxonomic rank
soil_genomes <- soil_genomes %>%
  mutate(rank = ifelse(species != "", "species",
                       ifelse(genus != "", "genus",
                              ifelse(family != "", "family",
                                     ifelse(order != "", "order",
                                            ifelse(class != "", "class",
                                                   ifelse(phylum != "", "phylum",
                                                          ifelse(domain != "", "domain", NA))))))))

# Determine GTDB taxon based on rank
soil_genomes <- soil_genomes %>%
  mutate(GTDB_taxon = ifelse(rank == "domain", domain,
                             ifelse(rank == "phylum", phylum,
                                    ifelse(rank == "class", class,
                                           ifelse(rank == "order", order,
                                                  ifelse(rank == "family", family,
                                                         ifelse(rank == "genus", genus,
                                                                ifelse(rank == "species", species, NA))))))))

# Merge with GTDB-NCBI mapping (already filtered to GTDB207)
log_message("Mapping GTDB taxonomy to NCBI taxids", source_name = source_name)
soil_genomes_ncbi <- merge(soil_genomes,
                           gtdb_ncbi_mapping,
                           by.x = "GTDB_taxon",
                           by.y = "GTDB_taxon",
                           all.x = TRUE)

# Determine higher classification for fallback taxonomy matching
soil_genomes_ncbi$higher_classification <- ifelse(soil_genomes_ncbi$rank == "species", soil_genomes_ncbi$genus,
                                                   ifelse(soil_genomes_ncbi$rank == "genus", soil_genomes_ncbi$family,
                                                          ifelse(soil_genomes_ncbi$rank == "family", soil_genomes_ncbi$order,
                                                                 ifelse(soil_genomes_ncbi$rank == "order", soil_genomes_ncbi$class,
                                                                        ifelse(soil_genomes_ncbi$rank == "class", soil_genomes_ncbi$phylum,
                                                                               ifelse(soil_genomes_ncbi$rank == "phylum", soil_genomes_ncbi$domain,
                                                                                      ifelse(soil_genomes_ncbi$rank == "domain", "Unclassified", NA)))))))

# Map higher classification and phylum to NCBI taxids
soil_genomes_ncbi$higher_NCBI_TaxID <- gtdb_ncbi_mapping[match(soil_genomes_ncbi$higher_classification,
                                                                 gtdb_ncbi_mapping$GTDB_taxon),]$ncbi_tax_id

soil_genomes_ncbi$phylum_NCBI_TaxID <- gtdb_ncbi_mapping[match(soil_genomes_ncbi$phylum,
                                                                 gtdb_ncbi_mapping$GTDB_taxon),]$ncbi_tax_id

# Combine all mapping strategies (prioritize in order)
soil_genomes_ncbi$NCBI_TaxID <- ifelse(!is.na(soil_genomes_ncbi$simple_NCBI_match), soil_genomes_ncbi$simple_NCBI_match,
                                       ifelse(!is.na(soil_genomes_ncbi$ncbi_tax_id), soil_genomes_ncbi$ncbi_tax_id,
                                              ifelse(!is.na(soil_genomes_ncbi$higher_NCBI_TaxID), soil_genomes_ncbi$higher_NCBI_TaxID,
                                                     ifelse(!is.na(soil_genomes_ncbi$phylum_NCBI_TaxID), soil_genomes_ncbi$phylum_NCBI_TaxID,
                                                            NA))))

# Report mapping success
mapped_count <- sum(!is.na(soil_genomes_ncbi$NCBI_TaxID))
log_message(paste("Successfully mapped", mapped_count, "of", nrow(soil_genomes_ncbi), "genomes to NCBI taxids"), source_name = source_name)

# Filter to genomes with NCBI taxids
soil_genomes_ncbi <- soil_genomes_ncbi %>%
  filter(!is.na(NCBI_TaxID))

log_message(paste("Processing", nrow(soil_genomes_ncbi), "genomes with NCBI taxids"), source_name = source_name)

# Prepare output
ready_genomes_spire <- soil_genomes_ncbi %>%
  mutate(source = "SPIRE_MAGs",
         is_published = "Y",
         ncbi_organism_name = paste0(GTDB_taxon, "_", genome_id),
         ncbi_species_taxid = as.numeric(NCBI_TaxID)) %>%
  select(ncbi_species_taxid,
         accession = genome_id,
         ncbi_organism_name,
         fasta_file_path = download_path,
         source,
         is_published) %>%
  arrange(ncbi_species_taxid)

# Add NCBI taxonomy information (required)
if (!dir.exists(NCBI_TAX_DIR)) {
  stop(paste("NCBI taxonomy directory is required but not found:", NCBI_TAX_DIR,
             "\nPlease ensure NCBI taxonomy files are available."))
}

log_message("Adding NCBI taxonomy information", source_name = source_name)
nodes <- getnodes(taxdir = NCBI_TAX_DIR)

ready_genomes_spire$rank <- CHNOSZ::getrank(ready_genomes_spire$ncbi_species_taxid, NCBI_TAX_DIR, nodes = nodes)
ncbi_taxonomy <- tax_id_to_ranked_lineage(ready_genomes_spire$ncbi_species_taxid, NCBI_TAX_DIR) %>%
  arrange(tax_id)
ncbi_taxonomy$kingdom <- ifelse(ncbi_taxonomy$kingdom == "", ncbi_taxonomy$superkingdom, ncbi_taxonomy$kingdom)

if (identical(ready_genomes_spire$ncbi_species_taxid, ncbi_taxonomy$tax_id)) {
  ready_genomes_spire <- cbind(ready_genomes_spire, ncbi_taxonomy) %>%
    mutate(species = ifelse(rank == "species", tax_name, species)) %>%
    mutate(phylum = ifelse(rank == "phylum", tax_name, phylum)) %>%
    mutate(class = ifelse(rank == "class", tax_name, class)) %>%
    mutate(genus = ifelse(rank == "genus", tax_name, genus)) %>%
    mutate(family = ifelse(rank == "family", tax_name, family)) %>%
    mutate(kingdom = ifelse(rank == "kingdom", tax_name, kingdom)) %>%
    mutate(order = ifelse(rank == "order", tax_name, order)) %>%
    mutate(ncbi_taxonomy = paste0("k__", kingdom,
                                   ";p__", phylum,
                                   ";c__", class,
                                   ";o__", order,
                                   ";f__", family,
                                   ";g__", genus,
                                   ";s__", species))
  ready_genomes_spire$ncbi_organism_name <- ready_genomes_spire$tax_name
  log_message("NCBI taxonomy information added", source_name = source_name)
} else {
  stop("Taxonomy ID mismatch - cannot proceed without valid NCBI taxonomy mapping")
}

log_message(paste("Writing", nrow(ready_genomes_spire), "genomes to output file"), source_name = source_name)
write_tsv(ready_genomes_spire, spire_config$output_file)

log_message(paste("SPIRE processing complete. Output written to:", spire_config$output_file), source_name = source_name)