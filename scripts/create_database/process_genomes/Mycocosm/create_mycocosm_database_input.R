# Process Mycocosm fungal genomes for Soil Microbe Database
# Processes manually downloaded genomes, adds taxonomy, and creates Struo2 input files
#
# Usage:
#   Rscript create_mycocosm_database_input.R
#   TEST_MODE=TRUE MAX_GENOMES=10 Rscript create_mycocosm_database_input.R  # Test with 10 genomes
#
# Note: Genome downloads from Mycocosm must be done manually. This script processes
# genomes that have already been downloaded to the genome directory.

library(tidyverse)
library(data.table)
library(CHNOSZ)

# Load centralized configuration
# Find the create_database directory (go up from process_genomes/Mycocosm)
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
# Go up to create_database directory (from process_genomes/Mycocosm)
config_path <- file.path(script_dir, "../..", "config.R")
if (!file.exists(config_path)) {
  # Try alternative path from project root
  config_path <- file.path(script_dir, "..", "..", "..", "..", "scripts", "create_database", "config.R")
}
if (!file.exists(config_path)) {
  stop("Cannot find config.R. Please ensure you're running from the project directory.")
}
source(config_path)

# Get Mycocosm-specific configuration
mycocosm_config <- get_source_config("Mycocosm")
source_name <- "Mycocosm"

# Load helper functions
if (!file.exists(HELPER_FUNCTIONS)) {
  stop(paste("Helper functions not found at:", HELPER_FUNCTIONS))
}
source(HELPER_FUNCTIONS)

log_message("Starting Mycocosm processing", source_name = source_name)

# Use local NCBI taxonomy directory (files should be copied via copy_ncbi_taxonomy_local.sh)
if (!file.exists(file.path(NCBI_TAX_DIR, "rankedlineage.dmp"))) {
  stop(paste("NCBI taxonomy directory not found: ", NCBI_TAX_DIR,
             "\nPlease copy NCBI taxonomy files using: bash scripts/create_database/copy_ncbi_taxonomy_local.sh"))
}

taxdir <- NCBI_TAX_DIR
log_message(paste("Using NCBI taxonomy directory:", taxdir), source_name = source_name)

# Read Mycocosm catalog from URL
log_message("Downloading Mycocosm catalog from JGI website", source_name = source_name)
tryCatch({
  # Read CSV with actual column names from the file
  myco_in <- read.csv(mycocosm_config$catalog_url, 
                      check.names = FALSE, 
                      stringsAsFactors = FALSE) %>%
    # Rename columns to match expected names (handle spaces and special characters)
    rename(row = `##`,
           Name = name,
           NCBI_TaxID = `NCBI Taxon`,
           gene_count = `#of genes`,
           is_public = `is public`,
           is_published = `is published`,
           is_superseded = `is superseded`,
           `superseded by` = `superseded by`,
           publications = `publication(s)`,
           pubmed_id = `pubmed id(s)`,
           doi_id = `doi id(s)`) %>%
    mutate(NCBI_TaxID = as.numeric(NCBI_TaxID)) %>%
    arrange(NCBI_TaxID)
  myco_in$Name <- gsub('\\"', "", myco_in$Name)
  myco_in$filename <- paste0(myco_in$portal, "_AssemblyScaffolds_Repeatmasked.fasta.gz")
  log_message(paste("Downloaded catalog with", nrow(myco_in), "genomes"), source_name = source_name)
}, error = function(e) {
  stop(paste("Failed to download Mycocosm catalog:", e$message,
             "\nURL:", mycocosm_config$catalog_url))
})

# Check for already downloaded genomes (local and remote)
log_message("Checking for downloaded genomes", source_name = source_name)
mycocosm_downloaded <- data.frame()

# Check local directory
if (dir.exists(mycocosm_config$local_genome_dir)) {
  local_genomes <- read_in_genomes(mycocosm_config$local_genome_dir, 
                                    pattern = mycocosm_config$genome_pattern)
  if (nrow(local_genomes) > 0) {
    mycocosm_downloaded <- rbind(mycocosm_downloaded, local_genomes)
    log_message(paste("Found", nrow(local_genomes), "genomes in local directory"), source_name = source_name)
  }
}

# Check genome directory (if different from local)
if ("genome_dir" %in% names(mycocosm_config) && mycocosm_config$genome_dir != mycocosm_config$local_genome_dir) {
  if (dir.exists(mycocosm_config$genome_dir)) {
    log_message(paste("Checking genome directory:", mycocosm_config$genome_dir), source_name = source_name)
    remote_genomes <- read_in_genomes(mycocosm_config$genome_dir, 
                                       pattern = mycocosm_config$genome_pattern)
    if (nrow(remote_genomes) > 0) {
      mycocosm_downloaded <- rbind(mycocosm_downloaded, remote_genomes)
      log_message(paste("Found", nrow(remote_genomes), "genomes in genome directory"), source_name = source_name)
    } else {
      log_message(paste("Genome directory exists but no genomes found matching pattern:", mycocosm_config$genome_pattern), source_name = source_name)
    }
  }
}

if (nrow(mycocosm_downloaded) == 0) {
  log_message("WARNING: No downloaded genomes found. Please download genomes manually from Mycocosm.", source_name = source_name)
  log_message("Genome directory locations checked:", source_name = source_name)
  log_message(paste("  Local:", mycocosm_config$local_genome_dir), source_name = source_name)
  if ("genome_dir" %in% names(mycocosm_config) && mycocosm_config$genome_dir != mycocosm_config$local_genome_dir) {
    log_message(paste("  Genome:", mycocosm_config$genome_dir), source_name = source_name)
  }
  stop("No genomes found. Please ensure genomes are downloaded to one of the above directories.")
}

# Remove duplicates if genomes exist in both local and remote
mycocosm_downloaded <- mycocosm_downloaded %>%
  distinct(filename, .keep_all = TRUE)

log_message(paste("Total unique genomes found:", nrow(mycocosm_downloaded)), source_name = source_name)

# Apply test mode limit if specified
if (is_test_mode() && !is.na(MAX_GENOMES)) {
  log_message(paste("TEST MODE: Limiting to", MAX_GENOMES, "genomes"), source_name = source_name)
  mycocosm_downloaded <- mycocosm_downloaded %>% slice_head(n = MAX_GENOMES)
}

# Match downloaded genomes with catalog
log_message("Matching downloaded genomes with catalog", source_name = source_name)
mycocosm_downloaded$portal <- myco_in[match(mycocosm_downloaded$filename, myco_in$filename), ]$portal
mycocosm_downloaded$is_published <- myco_in[match(mycocosm_downloaded$filename, myco_in$filename), ]$is_published
mycocosm_downloaded$NCBI_TaxID <- myco_in[match(mycocosm_downloaded$portal, myco_in$portal), ]$NCBI_TaxID %>% as.numeric()
mycocosm_downloaded$Name <- myco_in[match(mycocosm_downloaded$portal, myco_in$portal), ]$Name

# Filter to genomes with valid TaxIDs
genomes_with_taxid <- mycocosm_downloaded %>% filter(!is.na(NCBI_TaxID))
log_message(paste("Found", nrow(genomes_with_taxid), "genomes with valid NCBI TaxIDs"), source_name = source_name)

if (nrow(genomes_with_taxid) == 0) {
  stop("No genomes with valid NCBI TaxIDs found. Cannot proceed with taxonomy assignment.")
}

# Add NCBI taxonomy
log_message("Adding NCBI taxonomy information", source_name = source_name)

# Read ranked lineage file (should be local now after copying if needed)
ranked_lineage_file <- file.path(taxdir, "rankedlineage.dmp")
if (!file.exists(ranked_lineage_file)) {
  stop(paste("rankedlineage.dmp not found at:", ranked_lineage_file,
             "\nPlease ensure NCBI taxonomy files are available."))
}

log_message("Reading rankedlineage.dmp", source_name = source_name)
ncbi_taxdump <- data.table::fread(ranked_lineage_file,
                                   col.names = c("tax_id", "tax_name", "species", "genus", "family",
                                                 "order", "class", "phylum", "kingdom", "superkingdom", "NA"), 
                                   sep = "|")

# Get nodes for CHNOSZ
nodes <- getnodes(taxdir = taxdir)
genomes_with_taxid$rank <- CHNOSZ::getrank(genomes_with_taxid$NCBI_TaxID, taxdir, nodes = nodes)

# Get ranked lineage from already-read taxdump file
# rankedlineage.dmp has pipe-delimited format with trailing pipes, so we need to select every other column
# The file has 11 columns total (tax_id, tax_name, species, genus, family, order, class, phylum, kingdom, superkingdom, plus trailing empty)
ncbi_tax_ids_df <- data.frame(query_tax_id = as.numeric(genomes_with_taxid$NCBI_TaxID))
# Select every other column (1, 3, 5, 7, 9, 11) - this gives us the 10 data columns
ranked_lineage <- ncbi_taxdump %>%
  as_tibble() %>%
  select(seq(from = 1, to = min(11, ncol(ncbi_taxdump)), by = 2)) %>%
  `colnames<-`(c("tax_id", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"))

ncbi_taxonomy <- ranked_lineage %>%
  right_join(ncbi_tax_ids_df, by = c("tax_id" = "query_tax_id")) %>%
  arrange(tax_id)

genomes_with_taxid <- genomes_with_taxid %>% arrange(NCBI_TaxID) %>% mutate(NCBI_TaxID = as.numeric(NCBI_TaxID))

# Verify taxonomy matching (original approach)
check_match <- identical(ncbi_taxonomy$tax_id, genomes_with_taxid$NCBI_TaxID)
log_message(paste("Taxonomy good to join with TaxIDs?", check_match), source_name = source_name)

if (!check_match) {
  log_message("WARNING: Taxonomy ID order mismatch - attempting to join by TaxID", source_name = source_name)
  # Join by TaxID instead
  genomes_with_taxid <- genomes_with_taxid %>%
    left_join(ncbi_taxonomy, by = c("NCBI_TaxID" = "tax_id"))
} else {
  log_message("Taxonomy IDs match - joining directly", source_name = source_name)
  # Use cbind like the original script - exclude tax_id from ncbi_taxonomy to avoid duplicate
  ncbi_taxonomy_no_id <- ncbi_taxonomy %>% select(-tax_id)
  genomes_with_taxid <- cbind.data.frame(genomes_with_taxid, ncbi_taxonomy_no_id)
}

# Check which columns we have after join
log_message(paste("Columns after join:", paste(colnames(genomes_with_taxid), collapse = ", ")), source_name = source_name)

# Create custom taxonomy string (original approach)
# Handle case where kingdom might be empty - use superkingdom
if ("kingdom" %in% colnames(genomes_with_taxid)) {
  genomes_with_taxid$kingdom <- ifelse(is.na(genomes_with_taxid$kingdom) | genomes_with_taxid$kingdom == "", 
                                       genomes_with_taxid$superkingdom, 
                                       genomes_with_taxid$kingdom)
}

genomes_with_taxid <- genomes_with_taxid %>%
  mutate(species = ifelse(rank == "species", tax_name, species)) %>%
  mutate(custom_taxonomy = paste0("k__", ifelse(is.na(kingdom) | kingdom == "", 
                                                ifelse(is.na(superkingdom), "", superkingdom), 
                                                kingdom),
                                  ";p__", ifelse(is.na(phylum), "", phylum),
                                  ";c__", ifelse(is.na(class), "", class),
                                  ";o__", ifelse(is.na(order), "", order),
                                  ";f__", ifelse(is.na(family), "", family),
                                  ";g__", ifelse(is.na(genus), "", genus),
                                  ";s__", ifelse(is.na(species), "", species)))

# Prepare Struo2 input files
log_message("Preparing Struo2 input files", source_name = source_name)
ready_genomes_myco <- genomes_with_taxid %>%
  mutate(source = "Mycocosm") %>%
  select(ncbi_species_taxid = NCBI_TaxID,
         accession = Name,
         ncbi_organism_name = tax_name,
         species,
         fasta_file_path = filepath,
         ncbi_taxonomy = custom_taxonomy,
         source,
         is_published)

# Filter to published only
to_write_published <- ready_genomes_myco %>%
  filter(is_published == "Y")

# Write output files
log_message(paste("Writing", nrow(to_write_published), "published genomes to output file"), source_name = source_name)
write_tsv(to_write_published, mycocosm_config$output_file_published)

log_message(paste("Writing", nrow(ready_genomes_myco), "total genomes to output file"), source_name = source_name)
write_tsv(ready_genomes_myco, mycocosm_config$output_file_all)

# Also write to main output file (published only, matching other sources)
write_tsv(to_write_published, mycocosm_config$output_file)

log_message(paste("Mycocosm processing complete. Output files written to:", 
                  "\n  Published:", mycocosm_config$output_file_published,
                  "\n  All:", mycocosm_config$output_file_all,
                  "\n  Main:", mycocosm_config$output_file), source_name = source_name)
