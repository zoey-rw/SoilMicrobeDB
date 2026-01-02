# Process JGI GOLD (Genomes Online Database) genomes for Soil Microbe Database
# Downloads genomes, processes metadata, and creates Struo2 input file
#
# Usage:
#   Rscript download_jgi_gold_genomes.r
#   TEST_MODE=TRUE MAX_GENOMES=10 Rscript download_jgi_gold_genomes.r  # Test with 10 genomes

library(readxl)
library(tidyverse)
library(data.table)
library(jsonlite)
library(CHNOSZ)

# Load centralized configuration
# Find the create_database directory (go up from process_genomes/JGI_Gold)
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
# Go up to create_database directory (from process_genomes/JGI_Gold)
config_path <- file.path(script_dir, "../..", "config.R")
if (!file.exists(config_path)) {
  # Try alternative path from project root
  config_path <- file.path(script_dir, "..", "..", "..", "..", "scripts", "create_database", "config.R")
}
if (!file.exists(config_path)) {
  stop("Cannot find config.R. Please ensure you're running from the project directory.")
}
source(config_path)

# Get JGI GOLD-specific configuration
jgi_gold_config <- get_source_config("JGI_Gold")
source_name <- "JGI_Gold"

# Load helper functions
if (!file.exists(HELPER_FUNCTIONS)) {
  stop(paste("Helper functions not found at:", HELPER_FUNCTIONS))
}
source(HELPER_FUNCTIONS)

log_message("Starting JGI GOLD processing", source_name = source_name)

# Get JGI GOLD file paths from config (can be overridden by environment variables)
GOLD_DATA_FILE_CONFIG <- Sys.getenv("GOLD_DATA_FILE", unset = jgi_gold_config$data_file)
GOLD_ECOSYSTEM_FILE_CONFIG <- Sys.getenv("GOLD_ECOSYSTEM_FILE", unset = jgi_gold_config$ecosystem_file)

# Create local directory if it doesn't exist
JGI_GOLD_LOCAL_DIR <- get_source_dir("JGI_Gold")
JGI_GOLD_REF_DATA_DIR <- file.path(JGI_GOLD_LOCAL_DIR, "ref_data")
dir.create(JGI_GOLD_REF_DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# Local file paths
local_gold_data <- file.path(JGI_GOLD_REF_DATA_DIR, "goldData.xlsx")
local_ecosystem <- file.path(JGI_GOLD_REF_DATA_DIR, "GOLDs5levelEcosystemClassificationPaths.xlsx")

# Function to get GOLD data file (try local, then remote copy, then download)
get_gold_file <- function(remote_path, local_path, description, download_url = NULL) {
  # First, check if local file exists
  if (file.exists(local_path) && file.size(local_path) > 0) {
    log_message(paste("Using local", description), source_name = source_name)
    return(local_path)
  }
  
  # Try to copy from remote via SSH
  if (grepl("^/projectnb", remote_path)) {
    log_message(paste("Attempting to copy", description, "from remote server via SSH"), source_name = source_name)
    ssh_cmd <- paste0("scp zrwerbin@scc2.bu.edu:", shQuote(remote_path), " ", shQuote(local_path))
    result <- system(ssh_cmd, ignore.stderr = TRUE)
    if (result == 0 && file.exists(local_path) && file.size(local_path) > 0) {
      log_message(paste("Successfully copied", description, "from remote"), source_name = source_name)
      return(local_path)
    }
  }
  
  # Try to download from JGI GOLD website
  if (!is.null(download_url)) {
    log_message(paste("Downloading", description, "from JGI GOLD website..."), source_name = source_name)
    tryCatch({
      # Download using curl (more reliable for large files)
      curl_cmd <- paste0("curl -L -o ", shQuote(local_path), 
                         " -H 'User-Agent: Mozilla/5.0' ", shQuote(download_url))
      result <- system(curl_cmd)
      if (result == 0 && file.exists(local_path) && file.size(local_path) > 0) {
        log_message(paste("Successfully downloaded", description, "to:", local_path), source_name = source_name)
        return(local_path)
      } else {
        if (file.exists(local_path)) unlink(local_path)
        stop("Downloaded file is empty or download failed")
      }
    }, error = function(e) {
      if (file.exists(local_path)) unlink(local_path)
      stop(paste("Failed to download", description, ":", e$message,
                 "\nURL:", download_url,
                 "\nPlease download manually from: https://gold.jgi.doe.gov/download?mode=site_excel"))
    })
  }
  
  # Provide helpful error message
  error_msg <- paste0(description, " not found. Tried:\n",
                      "  Local: ", local_path, "\n",
                      "  Remote: ", remote_path)
  if (!is.null(download_url)) {
    error_msg <- paste0(error_msg, "\n  Download URL: ", download_url)
  }
  error_msg <- paste0(error_msg, 
                       "\n\nPlease download the GOLD data file manually from:",
                       "\n  https://gold.jgi.doe.gov/download?mode=site_excel",
                       "\n  Save it as: ", local_path,
                       "\n\nOr copy from remote server using:",
                       "\n  scp zrwerbin@scc2.bu.edu:", remote_path, " ", local_path)
  stop(error_msg)
}

# Get GOLD data files
# Note: JGI GOLD website requires manual download, so we'll try remote copy first
GOLD_DATA_FILE <- get_gold_file(GOLD_DATA_FILE_CONFIG, local_gold_data, "GOLD data file")
GOLD_ECOSYSTEM_FILE <- get_gold_file(GOLD_ECOSYSTEM_FILE_CONFIG, local_ecosystem, "GOLD ecosystem file")

# Read GOLD data files (only sheets we need)
log_message("Reading GOLD data files", source_name = source_name)
gold_data_organism <- readxl::read_xlsx(GOLD_DATA_FILE, sheet = 4)
gold_data_sequencing <- readxl::read_xlsx(GOLD_DATA_FILE, sheet = 5)
gold_data_analysis <- readxl::read_xlsx(GOLD_DATA_FILE, sheet = 6)

# Read ecosystem paths file
# Extract unique ecosystem paths from the Organism sheet (matching original approach)
log_message("Reading ecosystem paths file", source_name = source_name)
ecosystem_paths <- gold_data_organism %>%
  select(`ORGANISM ECOSYSTEM PATH ID`, `ORGANISM ECOSYSTEM TYPE`) %>%
  distinct(`ORGANISM ECOSYSTEM PATH ID`, .keep_all = TRUE)

# Filter for soil ecosystems (matching original script approach)
log_message("Filtering for soil ecosystems", source_name = source_name)
soil_paths <- ecosystem_paths %>% filter(grepl("Soil", `ORGANISM ECOSYSTEM TYPE`))
soil_path_ids <- soil_paths %>% select(`ORGANISM ECOSYSTEM PATH ID`) %>% unlist()

# Filter organisms to organisms found in soil
soil_organisms <- gold_data_organism %>% 
  filter(`ORGANISM ECOSYSTEM PATH ID` %in% soil_path_ids)

# Optional: Remove fungi (currently commented out)
# soil_organisms <- soil_organisms %>% filter(`ORGANISM NCBI SUPERKINGDOM` != "Eukaryota")

soil_organism_ids <- soil_organisms %>% select(`ORGANISM GOLD ID`) %>% unique() %>% unlist()
log_message(paste("Found", length(soil_organism_ids), "soil organisms"), source_name = source_name)

# Filter to genomes (analysis projects) from soil organisms
soil_sequencing <- gold_data_sequencing %>% filter(`ORGANISM GOLD ID` %in% soil_organism_ids)
soil_sequencing_ids <- soil_sequencing %>% select(`PROJECT GOLD ID`) %>% unlist()

# Filter to genomes (sequencing metadata) from soil organisms
soil_ap <- gold_data_analysis %>% filter(`AP PROJECT GOLD IDS` %in% soil_sequencing_ids)

# Merge sequencing and analysis data
genome_info <- merge(soil_ap, soil_sequencing, by.x = "AP PROJECT GOLD IDS", by.y = "PROJECT GOLD ID")

log_message(paste("Found", nrow(genome_info), "genomes from soil organisms"), source_name = source_name)

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
if (file.exists(GTDB_NCBI_MAPPING_FILE)) {
  ncbi_gtdb_mapping_all <- readRDS(GTDB_NCBI_MAPPING_FILE)
} else if (grepl("^/projectnb", GTDB_NCBI_MAPPING_FILE)) {
  # Read remote file via SSH
  log_message("Reading GTDB mapping file from remote server via SSH", source_name = source_name)
  remote_path <- GTDB_NCBI_MAPPING_FILE
  ssh_cmd <- paste0("ssh zrwerbin@scc2.bu.edu 'cat ", shQuote(remote_path), "'")
  temp_file <- tempfile(fileext = ".rds")
  system(paste(ssh_cmd, ">", temp_file))
  if (file.exists(temp_file) && file.size(temp_file) > 0) {
    ncbi_gtdb_mapping_all <- readRDS(temp_file)
    unlink(temp_file)
    log_message("Successfully loaded GTDB mapping from remote", source_name = source_name)
  } else {
    stop(paste("Failed to read GTDB-NCBI mapping file from remote:", GTDB_NCBI_MAPPING_FILE))
  }
} else {
  stop(paste("GTDB-NCBI mapping file not found:", GTDB_NCBI_MAPPING_FILE))
}

# Filter to GTDB 207 for JGI GOLD (can be adjusted if needed)
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
      stop(paste("GTDB 207 metadata files are required but not accessible. Tried:",
                 "\n  Local:", gtdb_207_bac_file, "and", gtdb_207_ar_file,
                 "\n  Remote:", remote_bac, "and", remote_ar,
                 "\nPlease ensure files are available locally or remotely accessible via SSH."))
    }
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
# JGI GOLD genomes are typically downloaded from NCBI, so check NCBI genome directory
NCBI_GENOME_REMOTE_BASE <- jgi_gold_config$remote_genome_dir
NCBI_GENOME_LOCAL_DIR <- jgi_gold_config$local_genome_dir

# Check both local and remote directories
already_downloaded <- data.frame()
if (dir.exists(NCBI_GENOME_LOCAL_DIR)) {
  already_downloaded <- read_in_genomes(NCBI_GENOME_LOCAL_DIR, pattern = "\\.fna\\.gz$")
}

# Also check remote if accessible
if (dir.exists(NCBI_GENOME_REMOTE_BASE)) {
  remote_downloaded <- read_in_genomes(NCBI_GENOME_REMOTE_BASE, pattern = "\\.fna\\.gz$")
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
if (!dir.exists(NCBI_TAX_DIR)) {
  stop(paste("NCBI taxonomy directory is required but not found:", NCBI_TAX_DIR,
             "\nPlease ensure NCBI taxonomy files are available."))
}

nodes <- getnodes(taxdir = NCBI_TAX_DIR)
genome_info$rank <- CHNOSZ::getrank(genome_info$ncbi_taxid, NCBI_TAX_DIR, nodes = nodes)
ncbi_taxonomy <- tax_id_to_ranked_lineage(genome_info$ncbi_taxid, NCBI_TAX_DIR) %>%
  arrange(tax_id)

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
ready_genomes_jgi <- genome_info %>%
  mutate(source = "JGI_GOLD",
         is_published = "Y",
         accession = NCBI_accession,
         ncbi_species_taxid = as.numeric(ncbi_taxid),
         fasta_file_path = file.path(NCBI_GENOME_LOCAL_DIR, paste0(gsub("\\.", "_", base_accession), "_genomic.fna.gz"))) %>%
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

log_message(paste("Writing", nrow(ready_genomes_jgi), "genomes to output file"), source_name = source_name)
write_tsv(ready_genomes_jgi, jgi_gold_config$output_file)

log_message(paste("JGI GOLD processing complete. Output written to:", jgi_gold_config$output_file), source_name = source_name)

