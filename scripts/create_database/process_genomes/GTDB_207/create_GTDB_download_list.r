# Process GTDB 207 genomes for Soil Microbe Database
# Filters GTDB 207 metadata, checks for existing downloads, subsamples genomes,
# and creates Struo2 input file
#
# Usage:
#   Rscript create_GTDB_download_list.r
#   TEST_MODE=TRUE MAX_GENOMES=10 Rscript create_GTDB_download_list.r  # Test with 10 genomes

library(data.table)
library(tidyverse)
library(CHNOSZ)

# Load configuration
# Find the create_database directory (go up from process_genomes/GTDB_207)
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
# Go up to create_database directory (from process_genomes/GTDB_207)
config_path <- file.path(script_dir, "../..", "config.R")
if (!file.exists(config_path)) {
  # Try alternative path from project root
  config_path <- file.path(script_dir, "..", "..", "..", "..", "scripts", "create_database", "config.R")
}
if (!file.exists(config_path)) {
  stop("Cannot find config.R. Please ensure you're running from the project directory.")
}
source(config_path)

# Get GTDB_207-specific configuration
gtdb_config <- get_source_config("GTDB_207")
source_name <- "GTDB_207"

# Load helper functions
if (!file.exists(HELPER_FUNCTIONS)) {
  stop(paste("Helper functions not found at:", HELPER_FUNCTIONS))
}
source(HELPER_FUNCTIONS)

log_message("Starting GTDB 207 processing", source_name = source_name)

# Load GTDB 207 metadata (required for taxonomy matching)
# This should be loaded by helper_functions.r, but verify it exists and load if needed
if (!exists("gtdb_207_metadata", envir = .GlobalEnv) || is.null(get("gtdb_207_metadata", envir = .GlobalEnv))) {
  # Try local paths first
  if (file.exists(GTDB_207_BAC_FILE) && file.exists(GTDB_207_AR_FILE)) {
    log_message("Loading GTDB 207 metadata from local files", source_name = source_name)
    gtdb_207_metadata_bac <- fread(GTDB_207_BAC_FILE) %>%
      select(accession, checkm_completeness, checkm_contamination, gtdb_genome_representative, gtdb_taxonomy,
             ncbi_date, ncbi_genbank_assembly_accession, ncbi_organism_name, ncbi_taxid, ncbi_taxonomy, ncbi_genome_category) %>%
      mutate(source = "GTDB_207")
    gtdb_207_metadata_ar <- fread(GTDB_207_AR_FILE) %>%
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
      stop("GTDB 207 metadata files are required but not found locally or remotely. Please ensure files are accessible.")
    }
  }
} else {
  gtdb_207_metadata <- get("gtdb_207_metadata", envir = .GlobalEnv)
  log_message("Using GTDB 207 metadata from global environment", source_name = source_name)
}

# Verify metadata is loaded and has data
if (is.null(gtdb_207_metadata) || nrow(gtdb_207_metadata) == 0) {
  stop("GTDB 207 metadata is empty or not loaded. Please ensure metadata files are accessible.")
}

log_message(paste("GTDB 207 metadata loaded:", nrow(gtdb_207_metadata), "genomes"), source_name = source_name)

# Check if required columns exist
required_cols <- c("gtdb_genome_representative", "checkm_completeness", "checkm_contamination", "ncbi_genbank_assembly_accession", "gtdb_taxonomy")
missing_cols <- setdiff(required_cols, colnames(gtdb_207_metadata))
if (length(missing_cols) > 0) {
  stop(paste("GTDB 207 metadata is missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Filter GTDB 207 metadata for high-quality representative genomes
# Filter criteria: gtdb_genome_representative == accession (genome is representative), completeness >= 90, contamination < 5
log_message("Filtering GTDB 207 metadata for high-quality representative genomes", source_name = source_name)
gtdb_filtered <- gtdb_207_metadata %>%
  filter(gtdb_genome_representative == accession & 
         checkm_completeness >= 90 & 
         checkm_contamination < 5) %>%
  mutate(base_accession = gsub(".1$|.2$", "", ncbi_genbank_assembly_accession))

# Extract kingdom from GTDB taxonomy
gtdb_filtered <- gtdb_filtered %>%
  mutate(kingdom = unlist(lapply(strsplit(gtdb_taxonomy, ";p"), "[[", 1)))
kingdom_table <- table(gtdb_filtered$kingdom)
log_message(paste("Filtered genomes by kingdom:", paste(names(kingdom_table), kingdom_table, sep = "=", collapse = ", ")), source_name = source_name)

# Create high-quality and medium-quality subsets
gtdb_hq <- gtdb_filtered %>% filter(checkm_completeness > 95)
gtdb_mq <- gtdb_filtered %>% filter(checkm_completeness > 50)

log_message(paste("High-quality genomes (>95% completeness):", nrow(gtdb_hq)), source_name = source_name)
log_message(paste("Medium-quality genomes (>50% completeness):", nrow(gtdb_mq)), source_name = source_name)

# Check for already downloaded genomes to avoid duplicates
# Check both local and remote genome directories
log_message("Checking for already downloaded genomes", source_name = source_name)

# Check local NCBI genomes directory
local_ncbi_dir <- JGI_GOLD_CONFIG$local_genome_dir
jgi_downloaded_local <- read_in_genomes(local_ncbi_dir, pattern = ".fna.gz")

# Check remote NCBI genomes directory (if accessible)
remote_ncbi_dir <- JGI_GOLD_CONFIG$remote_genome_dir
jgi_downloaded_remote <- read_in_genomes(remote_ncbi_dir, pattern = ".fna.gz")

# Combine downloaded genomes
jgi_downloaded <- rbind(jgi_downloaded_local, jgi_downloaded_remote)

# Process downloaded genomes if any exist
if (nrow(jgi_downloaded) > 0) {
  jgi_downloaded <- jgi_downloaded %>%
    mutate(NCBI_accession = lapply(filename, function(x) {
      x1 <- strsplit(x, "_") %>% unlist
      return(paste0(x1[[1]], "_", x1[[2]]))
    }) %>% unlist) %>%
    mutate(base_accession = gsub(".1$|.2$", "", NCBI_accession)) %>%
    distinct(NCBI_accession, .keep_all = TRUE)
} else {
  # Create empty data frame with required columns
  jgi_downloaded <- data.frame(
    filename = character(0),
    base_accession = character(0),
    NCBI_accession = character(0),
    stringsAsFactors = FALSE
  )
}

log_message(paste("Found", nrow(jgi_downloaded), "already downloaded genomes"), source_name = source_name)

# Filter out already downloaded genomes
gtdb_mq$already_downloaded <- ifelse(gtdb_mq$base_accession %in% jgi_downloaded$base_accession, 1, 0)
gtdb_mq_subset <- gtdb_mq %>%
  filter(!already_downloaded)

log_message(paste("Genomes remaining after filtering duplicates:", nrow(gtdb_mq_subset)), source_name = source_name)

if (nrow(gtdb_mq_subset) == 0) {
  stop("No genomes remaining after filtering duplicates. All genomes may already be downloaded.")
}

# Apply test mode limit if specified (before subsampling)
if (is_test_mode() && !is.na(MAX_GENOMES)) {
  log_message(paste("TEST MODE: Limiting to", MAX_GENOMES, "genomes before subsampling"), source_name = source_name)
  gtdb_mq_subset <- gtdb_mq_subset %>% slice_head(n = MAX_GENOMES)
}

# Subsample to max 5 random genomes per NCBI taxid
max_genome_count <- 5
set.seed(1)

ngenomes <- list()
multigenomes <- list()
sampled <- list()

log_message("Subsampling genomes (max 5 per taxid)", source_name = source_name)
unique_taxids <- unique(gtdb_mq_subset$ncbi_taxid)
log_message(paste("Processing", length(unique_taxids), "unique taxids"), source_name = source_name)

for (i in unique_taxids) {
  subset <- gtdb_mq_subset[gtdb_mq_subset$ncbi_taxid == i, ]
  ngenomes[[i]] <- nrow(subset)
  if (nrow(subset) > max_genome_count) {
    multigenomes[[i]] <- subset
    sampled[[i]] <- subset %>% sample_n(max_genome_count)
  }
}

ngenomes_count <- unlist(ngenomes)
if (length(multigenomes) > 0) {
  multigenomes <- data.table::rbindlist(multigenomes)
} else {
  multigenomes <- data.frame()
}

if (length(sampled) > 0) {
  sampled <- data.table::rbindlist(sampled)
  sampled_taxids <- unique(sampled$ncbi_taxid)
  orig <- gtdb_mq_subset %>% filter(!ncbi_taxid %in% sampled_taxids)
  combined_subsampled <- rbind(orig, sampled)
} else {
  # No subsampling needed - all taxids have <= max_genome_count genomes
  combined_subsampled <- gtdb_mq_subset
}

log_message(paste("After subsampling:", nrow(combined_subsampled), "genomes"), source_name = source_name)

# Get NCBI taxonomy information
log_message("Loading NCBI taxonomy information", source_name = source_name)

# Check if NCBI taxonomy directory exists
# If NCBI_TAX_DIR is a remote path (starts with /projectnb), assume it's accessible via mount/SSH
if (!dir.exists(NCBI_TAX_DIR) && !grepl("^/projectnb", NCBI_TAX_DIR)) {
  # Try remote path if local doesn't exist and it's not already a remote path
  remote_tax_dir <- file.path(GTDB_REMOTE_BASE, "ncbi_taxonomy")
  if (grepl("^/projectnb", remote_tax_dir)) {
    # For remote paths, we'll try to use them directly (assuming they're mounted or accessible)
    log_message(paste("Local NCBI taxonomy directory not found, will attempt to use remote:", remote_tax_dir), source_name = source_name)
    NCBI_TAX_DIR <- remote_tax_dir
  } else {
    stop(paste("NCBI taxonomy directory not found:", NCBI_TAX_DIR,
               "\nPlease set NCBI_TAX_DIR environment variable or ensure NCBI taxonomy files are available."))
  }
}

log_message(paste("Using NCBI taxonomy directory:", NCBI_TAX_DIR), source_name = source_name)

nodes <- getnodes(taxdir = NCBI_TAX_DIR)
combined_subsampled$rank <- CHNOSZ::getrank(combined_subsampled$ncbi_taxid, NCBI_TAX_DIR, nodes = nodes)
ncbi_taxonomy <- tax_id_to_ranked_lineage(combined_subsampled$ncbi_taxid, NCBI_TAX_DIR) %>%
  arrange(tax_id)

combined_subsampled <- combined_subsampled %>% 
  arrange(ncbi_taxid) %>% 
  mutate(ncbi_taxid = as.numeric(ncbi_taxid))

check_match <- identical(ncbi_taxonomy$tax_id, combined_subsampled$ncbi_taxid)
ncbi_taxonomy$kingdom <- ncbi_taxonomy$superkingdom

log_message(paste("Taxonomy good to join with TaxIDs?", check_match), source_name = source_name)

if (!check_match) {
  warning("Taxonomy tax_id mismatch - proceeding with join anyway", source_name = source_name)
}

ncbi_taxonomy <- cbind.data.frame(combined_subsampled, ncbi_taxonomy) %>%
  mutate(species = ifelse(rank == "species", tax_name, species)) %>%
  mutate(custom_taxonomy =
           paste0("k__", kingdom,
                  ";p__", phylum,
                  ";c__", class,
                  ";o__", order,
                  ";f__", family,
                  ";g__", genus,
                  ";s__", species))

# Prepare final output (filter to published only - all GTDB genomes are published)
log_message("Preparing final output", source_name = source_name)
to_write <- inner_join(combined_subsampled,
                       ncbi_taxonomy,
                       relationship = "many-to-many") %>%
  mutate(source = "GTDB_207",
         is_published = "Y") %>%
  select(ncbi_species_taxid = ncbi_taxid,
         accession = ncbi_genbank_assembly_accession,
         ncbi_genbank_assembly_accession = ncbi_genbank_assembly_accession,
         ncbi_organism_name,
         kingdom,
         phylum,
         species,
         ncbi_taxonomy,
         source, 
         is_published)

log_message(paste("Final output contains", nrow(to_write), "genomes"), source_name = source_name)

# Write output file
output_file <- gtdb_config$output_file
log_message(paste("Writing output to:", output_file), source_name = source_name)

# Ensure output directory exists
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

write_tsv(to_write, output_file)

log_message("GTDB 207 processing completed successfully", source_name = source_name)
log_message(paste("Output file:", output_file), source_name = source_name)
