# Configuration file for database creation scripts
# Source this file at the beginning of all processing scripts
#
# Usage:
#   source("scripts/create_database/config.R")
#   # Or from within a subdirectory:
#   source(file.path(dirname(getwd()), "config.R"))
#
# Setup Instructions:
#   1. Set the GitHub repository URL (if hosting reference data on GitHub):
#      - Set environment variable: SOIL_MICROBE_DB_GITHUB
#      - Or edit GITHUB_REPO_BASE below
#   2. Place large genome files in the appropriate directories (see below)
#   3. Small reference data files will be downloaded from GitHub if not found locally

# Base directory - can be overridden by environment variable
# Automatically detects project root if running from within project
BASE_DIR <- Sys.getenv("BASE_DIR", unset = NA)
if (is.na(BASE_DIR)) {
  # Try to find project root by looking for common markers
  current_dir <- getwd()
  if (grepl("create_database", current_dir)) {
    # If we're in create_database or a subdirectory, go up to project root
    # From process_genomes/SOURCE: go up 4 levels (SPIRE -> process_genomes -> create_database -> scripts -> project_root)
    # From create_database: go up 2 levels (create_database -> scripts -> project_root)
    if (grepl("process_genomes", current_dir)) {
      BASE_DIR <- normalizePath(file.path(current_dir, "../../../.."))
    } else {
      BASE_DIR <- normalizePath(file.path(current_dir, "../.."))
    }
  } else if (grepl("scripts", current_dir)) {
    BASE_DIR <- normalizePath(file.path(current_dir, ".."))
  } else {
    BASE_DIR <- getwd()
  }
}
BASE_DIR <- normalizePath(BASE_DIR)

# Main directory structure
DATA_DIR <- file.path(BASE_DIR, "data")
GENOME_DB_DIR <- file.path(DATA_DIR, "genome_database")
STRUO2_INPUT_DIR <- file.path(GENOME_DB_DIR, "struo2_input_tables")
SCRIPTS_DIR <- file.path(BASE_DIR, "scripts")
CREATE_DB_DIR <- file.path(SCRIPTS_DIR, "create_database")
LOG_BASE_DIR <- file.path(CREATE_DB_DIR, "logs")

# Create base directories if they don't exist
dir.create(GENOME_DB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(STRUO2_INPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_BASE_DIR, showWarnings = FALSE, recursive = TRUE)

# Testing mode settings
# Set MAX_GENOMES to limit processing for testing (NULL = process all)
MAX_GENOMES <- as.numeric(Sys.getenv("MAX_GENOMES", unset = NA))
TEST_MODE <- as.logical(Sys.getenv("TEST_MODE", unset = "FALSE"))

# Default quality thresholds (can be overridden per source)
MIN_COMPLETENESS <- 95
MAX_CONTAMINATION <- 5

# Helper functions path - use absolute path to avoid issues
HELPER_FUNCTIONS <- normalizePath(file.path(CREATE_DB_DIR, "helper_functions.r"), mustWork = FALSE)

# Load helper functions (must be loaded after directories are defined)
source(HELPER_FUNCTIONS)

# ============================================================================
# Reference Data Files
# ============================================================================
# Small reference data files can be downloaded from GitHub if not found locally.
# Large genome files should be placed in the appropriate directories.

# GitHub repository for reference data
GITHUB_REPO_BASE <- Sys.getenv("SOIL_MICROBE_DB_GITHUB", 
                                unset = "https://raw.githubusercontent.com/zoey-rw/SoilMicrobeDB/master/data/genome_database")

# GTDB mapping files
GTDB_MAPPING_LOCAL_DIR <- file.path(GENOME_DB_DIR, "gtdb_mapping")
dir.create(GTDB_MAPPING_LOCAL_DIR, showWarnings = FALSE, recursive = TRUE)

# GTDB-NCBI mapping file (contains all GTDB versions: 214, 95, and 207)
# Scripts filter this file as needed (e.g., SPIRE filters to GTDB207, SMAG uses all versions)
GTDB_NCBI_MAPPING_FILE <- Sys.getenv("GTDB_NCBI_MAPPING_FILE", unset = NA)
if (is.na(GTDB_NCBI_MAPPING_FILE)) {
  local_mapping <- file.path(GTDB_MAPPING_LOCAL_DIR, "GTDB_NCBI_key.rds")
  # Try local first, then download from GitHub if available
  if (!file.exists(local_mapping)) {
    github_url <- paste0(GITHUB_REPO_BASE, "/gtdb_mapping/GTDB_NCBI_key.rds")
    tryCatch({
      download_from_github(github_url, local_mapping, "GTDB-NCBI mapping file")
    }, error = function(e) {
      message("Note: GTDB-NCBI mapping file not found locally or on GitHub.")
      message("Please download it manually and place at: ", local_mapping)
    })
  }
  GTDB_NCBI_MAPPING_FILE <- local_mapping
}

# GTDB 207 metadata files (required for taxonomy matching)
# These are large files - place them in data/genome_database/ or set via environment variable
GTDB_207_BAC_FILE <- Sys.getenv("GTDB_207_BAC_FILE", unset = file.path(GENOME_DB_DIR, "bac120_metadata_r207.tsv"))
GTDB_207_AR_FILE <- Sys.getenv("GTDB_207_AR_FILE", unset = file.path(GENOME_DB_DIR, "ar53_metadata_r207.tsv"))

# NCBI taxonomy directory
# Set via environment variable or use default local directory
NCBI_TAX_DIR <- Sys.getenv("NCBI_TAX_DIR", unset = file.path(GENOME_DB_DIR, "ncbi_taxonomy"))

# Load helper functions (must be loaded after directories are defined)
source(HELPER_FUNCTIONS)

# Source-specific configurations
# These can be overridden or extended by individual scripts

# SMAG-specific settings
SMAG_CONFIG <- list(
  metadata_file = file.path(get_source_dir("SMAG"), "Supplementary Data 2.xlsx"),
  download_base_url = "https://smag.microbmalab.cn/industrialshow/download/fileDownload?fileKey=SMAG/MAG/",
  output_file = file.path(STRUO2_INPUT_DIR, "SMAG_struo.tsv"),
  genome_pattern = ".fa"
)

# SPIRE-specific settings
SPIRE_LOCAL_DIR <- get_source_dir("SPIRE")
dir.create(SPIRE_LOCAL_DIR, showWarnings = FALSE, recursive = TRUE)

# Metadata files - try local first, then download from GitHub if available
SPIRE_METADATA_PATH <- Sys.getenv("SPIRE_METADATA_FILE", unset = NA)
if (is.na(SPIRE_METADATA_PATH)) {
  local_meta <- file.path(SPIRE_LOCAL_DIR, "spire_v1_genome_metadata.tsv.gz")
  if (!file.exists(local_meta)) {
    github_url <- paste0(GITHUB_REPO_BASE, "/SPIRE/spire_v1_genome_metadata.tsv.gz")
    tryCatch({
      download_from_github(github_url, local_meta, "SPIRE metadata file")
    }, error = function(e) {
      message("Note: SPIRE metadata file not found. Please download manually.")
    })
  }
  SPIRE_METADATA_PATH <- local_meta
}

SPIRE_TAXONOMY_PATH <- Sys.getenv("SPIRE_TAXONOMY_FILE", unset = NA)
if (is.na(SPIRE_TAXONOMY_PATH)) {
  local_tax <- file.path(SPIRE_LOCAL_DIR, "spire.per_cluster.taxonomy.tsv")
  if (!file.exists(local_tax)) {
    github_url <- paste0(GITHUB_REPO_BASE, "/SPIRE/spire.per_cluster.taxonomy.tsv")
    tryCatch({
      download_from_github(github_url, local_tax, "SPIRE taxonomy file")
    }, error = function(e) {
      message("Note: SPIRE taxonomy file not found. Please download manually.")
    })
  }
  SPIRE_TAXONOMY_PATH <- local_tax
}

SPIRE_MICRONTOLOGY_PATH <- Sys.getenv("SPIRE_MICRONTOLOGY_FILE", unset = NA)
if (is.na(SPIRE_MICRONTOLOGY_PATH)) {
  local_micro <- file.path(SPIRE_LOCAL_DIR, "spire_v1_microntology.tsv")
  if (!file.exists(local_micro)) {
    github_url <- paste0(GITHUB_REPO_BASE, "/SPIRE/spire_v1_microntology.tsv")
    tryCatch({
      download_from_github(github_url, local_micro, "SPIRE microntology file")
    }, error = function(e) {
      message("Note: SPIRE microntology file not found. Please download manually.")
    })
  }
  SPIRE_MICRONTOLOGY_PATH <- local_micro
}

# Genome directory (for checking existing downloads)
SPIRE_GENOME_DIR <- Sys.getenv("SPIRE_GENOME_DIR", unset = SPIRE_LOCAL_DIR)

SPIRE_CONFIG <- list(
  metadata_file = SPIRE_METADATA_PATH,
  taxonomy_file = SPIRE_TAXONOMY_PATH,
  microntology_file = SPIRE_MICRONTOLOGY_PATH,
  download_base_url = "http://spire.embl.de/download_file/",
  output_file = file.path(STRUO2_INPUT_DIR, "spire_MAGs_struo.tsv"),
  genome_pattern = "\\.fna\\.gz$",
  soil_habitat_keywords = c("soil", "rhizo", "forest", "cropland", "terrestrial", "litter", "agriculture"),
  genome_dir = SPIRE_GENOME_DIR
)

# GTDB 207-specific settings
GTDB_207_CONFIG <- list(
  output_file = file.path(STRUO2_INPUT_DIR, "gtdb_207_struo.tsv")
)

# Mycocosm-specific settings
MYCOCOSM_LOCAL_DIR <- get_source_dir("Mycocosm")
dir.create(MYCOCOSM_LOCAL_DIR, showWarnings = FALSE, recursive = TRUE)

# Genome directory - can be set via environment variable for custom locations
MYCOCOSM_GENOME_DIR <- Sys.getenv("MYCOCOSM_GENOME_DIR", unset = MYCOCOSM_LOCAL_DIR)

MYCOCOSM_CONFIG <- list(
  catalog_url = "https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/download-group?flt=&seq=all&pub=all&grp=fungi&srt=released&ord=desc",
  output_file = file.path(STRUO2_INPUT_DIR, "mycocosm_struo.tsv"),
  output_file_published = file.path(STRUO2_INPUT_DIR, "mycocosm_published_struo.tsv"),
  output_file_all = file.path(STRUO2_INPUT_DIR, "mycocosm_all_struo.tsv"),
  genome_pattern = ".fasta.gz",
  genome_dir = MYCOCOSM_GENOME_DIR,
  local_genome_dir = MYCOCOSM_LOCAL_DIR
)

# JGI GOLD-specific settings
JGI_GOLD_LOCAL_DIR <- get_source_dir("JGI_Gold")
JGI_GOLD_REF_DATA_DIR <- file.path(JGI_GOLD_LOCAL_DIR, "ref_data")
dir.create(JGI_GOLD_REF_DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# Metadata files - try local first, then download from GitHub if available
JGI_GOLD_DATA_FILE <- Sys.getenv("JGI_GOLD_DATA_FILE", unset = NA)
if (is.na(JGI_GOLD_DATA_FILE)) {
  local_data <- file.path(JGI_GOLD_REF_DATA_DIR, "goldData.xlsx")
  if (!file.exists(local_data)) {
    github_url <- paste0(GITHUB_REPO_BASE, "/JGI_Gold/ref_data/goldData.xlsx")
    tryCatch({
      download_from_github(github_url, local_data, "JGI GOLD data file")
    }, error = function(e) {
      message("Note: JGI GOLD data file not found. Please download from https://gold.jgi.doe.gov/")
    })
  }
  JGI_GOLD_DATA_FILE <- local_data
}

JGI_GOLD_ECOSYSTEM_FILE <- Sys.getenv("JGI_GOLD_ECOSYSTEM_FILE", unset = NA)
if (is.na(JGI_GOLD_ECOSYSTEM_FILE)) {
  local_ecosystem <- file.path(JGI_GOLD_REF_DATA_DIR, "GOLDs5levelEcosystemClassificationPaths.xlsx")
  if (!file.exists(local_ecosystem)) {
    github_url <- paste0(GITHUB_REPO_BASE, "/JGI_Gold/ref_data/GOLDs5levelEcosystemClassificationPaths.xlsx")
    tryCatch({
      download_from_github(github_url, local_ecosystem, "JGI GOLD ecosystem file")
    }, error = function(e) {
      message("Note: JGI GOLD ecosystem file not found. Please download from https://gold.jgi.doe.gov/")
    })
  }
  JGI_GOLD_ECOSYSTEM_FILE <- local_ecosystem
}

# NCBI genome directory (shared with RefSoil)
NCBI_GENOME_DIR <- Sys.getenv("NCBI_GENOME_DIR", unset = file.path(GENOME_DB_DIR, "ncbi_genomes"))
dir.create(NCBI_GENOME_DIR, showWarnings = FALSE, recursive = TRUE)

JGI_GOLD_CONFIG <- list(
  data_file = JGI_GOLD_DATA_FILE,
  ecosystem_file = JGI_GOLD_ECOSYSTEM_FILE,
  output_file = file.path(STRUO2_INPUT_DIR, "jgi_gold_struo.tsv"),
  genome_dir = NCBI_GENOME_DIR,
  local_genome_dir = NCBI_GENOME_DIR
)

# GEM/Nayfach-specific settings
GEM_LOCAL_DIR <- get_source_dir("GEM")
dir.create(GEM_LOCAL_DIR, showWarnings = FALSE, recursive = TRUE)

# Metadata files - try local first, then download from GitHub if available
GEM_SUPPLEMENTARY_FILE <- Sys.getenv("GEM_SUPPLEMENTARY_FILE", unset = NA)
if (is.na(GEM_SUPPLEMENTARY_FILE)) {
  local_supp <- file.path(GEM_LOCAL_DIR, "41587_2020_718_MOESM3_ESM.xlsx")
  if (!file.exists(local_supp)) {
    github_url <- paste0(GITHUB_REPO_BASE, "/GEM/41587_2020_718_MOESM3_ESM.xlsx")
    tryCatch({
      download_from_github(github_url, local_supp, "GEM supplementary file")
    }, error = function(e) {
      message("Note: GEM supplementary file not found. Please download from the publication.")
    })
  }
  GEM_SUPPLEMENTARY_FILE <- local_supp
}

GEM_CONFIG <- list(
  output_file = file.path(STRUO2_INPUT_DIR, "nayfach_MAGs_struo.tsv"),
  supplementary_file = GEM_SUPPLEMENTARY_FILE,
  download_base_url = "https://portal.nersc.gov/GEM/genomes/fna/",
  genome_pattern = "\\.fna\\.gz$",
  genome_dir = GEM_LOCAL_DIR,
  local_genome_dir = GEM_LOCAL_DIR
)

# RefSoil-specific settings
# RefSoil is a curated list of soil organisms from JGI GOLD
REFSOIL_LOCAL_DIR <- get_source_dir("RefSoil")
REFSOIL_REF_DATA_DIR <- file.path(REFSOIL_LOCAL_DIR, "ref_data")
dir.create(REFSOIL_REF_DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# Metadata files - try local first, then download from GitHub if available
REFSOIL_DATA_FILE <- Sys.getenv("REFSOIL_DATA_FILE", unset = NA)
if (is.na(REFSOIL_DATA_FILE)) {
  local_data <- file.path(REFSOIL_REF_DATA_DIR, "RefSoil_s1.xlsx")
  if (!file.exists(local_data)) {
    github_url <- paste0(GITHUB_REPO_BASE, "/RefSoil/ref_data/RefSoil_s1.xlsx")
    tryCatch({
      download_from_github(github_url, local_data, "RefSoil data file")
    }, error = function(e) {
      message("Note: RefSoil data file not found. Please download from the publication.")
    })
  }
  REFSOIL_DATA_FILE <- local_data
}

REFSOIL_CONFIG <- list(
  data_file = REFSOIL_DATA_FILE,
  output_file = file.path(STRUO2_INPUT_DIR, "refsoil_struo.tsv"),
  genome_dir = NCBI_GENOME_DIR,
  local_genome_dir = NCBI_GENOME_DIR
)

# Helper functions are now in helper_functions.r (loaded above)

