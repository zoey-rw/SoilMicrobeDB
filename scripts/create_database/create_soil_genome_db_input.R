# combine struo input files to create a fungi - prokaryote soil database
# This script combines all available struo input files from different genome sources
# into a single database file for use with Struo2

library(readr)
library(dplyr)
library(data.table)

# Load configuration
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg)
    dirname(normalizePath(script_path))
  } else {
    getwd()
  }
}, error = function(e) {
  getwd()
})
if (length(script_dir) == 0 || script_dir == "." || is.na(script_dir)) {
  script_dir <- getwd()
}
config_path <- file.path(script_dir, "config.R")
if (!file.exists(config_path)) {
  # Try alternative path from project root
  config_path <- file.path(script_dir, "..", "..", "scripts", "create_database", "config.R")
}
if (!file.exists(config_path)) {
  stop("Cannot find config.R. Please ensure you're running from the project directory.")
}
source(config_path)

# Load helper functions
if (!file.exists(HELPER_FUNCTIONS)) {
  stop(paste("Helper functions not found at:", HELPER_FUNCTIONS))
}
source(HELPER_FUNCTIONS)

# Read struo input files from all sources
# Function to safely read struo files with error handling
read_struo_file <- function(file_path, source_name) {
  if (!file.exists(file_path)) {
    warning(paste("Struo file not found for", source_name, ":", file_path, "\nSkipping this source."))
    return(NULL)
  }
  
  # Check file size - warn if suspiciously small (might be empty or header-only)
  file_size <- file.size(file_path)
  if (file_size < 1000) {
    warning(paste("Warning: Struo file for", source_name, "is very small (", file_size, "bytes).",
                  "It may be empty or contain only headers."))
  }
  
  tryCatch({
    df <- read_tsv(file_path, show_col_types = FALSE)
    
    # Check if file is empty or has no data rows
    if (nrow(df) == 0) {
      warning(paste("Struo file for", source_name, "is empty (no data rows). Skipping this source."))
      return(NULL)
    }
    
    # Check for required columns (fasta_file_path may be optional for some sources)
    required_cols <- c("ncbi_species_taxid", "accession")
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
      warning(paste("Struo file for", source_name, "is missing required columns:", 
                    paste(missing_cols, collapse = ", "), ". Skipping this source."))
      return(NULL)
    }
    
    # Add fasta_file_path if missing (some sources may not have it)
    if (!"fasta_file_path" %in% colnames(df)) {
      warning(paste("Struo file for", source_name, "is missing fasta_file_path column.",
                    "Adding placeholder (genome data may be available remotely)."))
      df$fasta_file_path <- NA_character_
    }
    
    cat("  ✓ Loaded", nrow(df), "genomes from", source_name, 
        paste0("(", format(file_size, big.mark = ","), " bytes)\n"))
    return(df)
  }, error = function(e) {
    warning(paste("Error reading struo file for", source_name, ":", e$message, "\nSkipping this source."))
    return(NULL)
  })
}

cat("Reading struo input files from all sources...\n")
cat("Note: Using existing struo files (genome data is available remotely if needed)\n")
cat("If files are small, regenerate complete files by running processing scripts without TEST_MODE\n\n")

# Read all available sources
# Note: If files are small (test files), regenerate complete files by running:
#   - GEM: Rscript scripts/create_database/process_genomes/GEM/download_nayfach_MAGs.r
#   - SPIRE: Rscript scripts/create_database/process_genomes/SPIRE/download_spire_MAGs.r  
#   - Mycocosm: Rscript scripts/create_database/process_genomes/Mycocosm/create_mycocosm_database_input.R
# (Run without TEST_MODE to process all genomes from remote directories)

nayfach_struo <- read_struo_file(GEM_CONFIG$output_file, "GEM") %>% 
  {if (!is.null(.)) mutate(., source = "GEM catalog") else .}

# Warn if GEM file is very small (likely test file)
if (!is.null(nayfach_struo) && nrow(nayfach_struo) < 100) {
  warning("GEM struo file contains only ", nrow(nayfach_struo), " genomes. This appears to be a test file.\n",
          "  To regenerate the complete file, run:\n",
          "  Rscript scripts/create_database/process_genomes/GEM/download_nayfach_MAGs.r\n",
          "  (without TEST_MODE to process all genomes from remote directories)")
}

SMAG_struo <- read_struo_file(SMAG_CONFIG$output_file, "SMAG")

JGI_struo <- read_struo_file(JGI_GOLD_CONFIG$output_file, "JGI GOLD")

# Use "all" version for Mycocosm if it's substantial, otherwise use published
myco_struo <- NULL
if (file.exists(MYCOCOSM_CONFIG$output_file_all) && file.size(MYCOCOSM_CONFIG$output_file_all) > 10000) {
  myco_struo <- read_struo_file(MYCOCOSM_CONFIG$output_file_all, "Mycocosm (all)")
} else {
  myco_struo <- read_struo_file(MYCOCOSM_CONFIG$output_file_published, "Mycocosm (published)")
}

# Warn if Mycocosm file is very small
if (!is.null(myco_struo) && nrow(myco_struo) < 100) {
  warning("Mycocosm struo file contains only ", nrow(myco_struo), " genomes. This appears to be a test file.\n",
          "  To regenerate the complete file, run:\n",
          "  Rscript scripts/create_database/process_genomes/Mycocosm/create_mycocosm_database_input.R\n",
          "  (without TEST_MODE to process all genomes from remote directories)")
}

spire_struo <- read_struo_file(SPIRE_CONFIG$output_file, "SPIRE")

# Warn if SPIRE file is very small
if (!is.null(spire_struo) && nrow(spire_struo) < 100) {
  warning("SPIRE struo file contains only ", nrow(spire_struo), " genomes. This appears to be a test file.\n",
          "  To regenerate the complete file, run:\n",
          "  Rscript scripts/create_database/process_genomes/SPIRE/download_spire_MAGs.r\n",
          "  (without TEST_MODE to process all genomes from remote directories)")
}

refsoil_struo <- read_struo_file(REFSOIL_CONFIG$output_file, "RefSoil")

# GTDB 207 - prefer filtered (full) version if available, otherwise use subsampled version
gtdb_207_filtered_file <- file.path(STRUO2_INPUT_DIR, "gtdb_207_filtered_struo.tsv")
if (file.exists(gtdb_207_filtered_file) && file.size(gtdb_207_filtered_file) > 1000000) {
  gtdb_207_struo <- read_struo_file(gtdb_207_filtered_file, "GTDB 207 (filtered)")
} else {
  gtdb_207_struo <- read_struo_file(GTDB_207_CONFIG$output_file, "GTDB 207")
}

# Combine all non-NULL data frames
struo_list <- list(nayfach_struo, SMAG_struo, JGI_struo, myco_struo, spire_struo, refsoil_struo, gtdb_207_struo)
struo_list <- struo_list[!sapply(struo_list, is.null)]

if (length(struo_list) == 0) {
  stop("No struo input files found! Please ensure at least one source has been processed.")
}

cat("\nCombining", length(struo_list), "sources...\n")

# Use local NCBI taxonomy directory
if (!file.exists(file.path(NCBI_TAX_DIR, "nodes.dmp"))) {
  stop(paste("NCBI taxonomy directory not found: ", NCBI_TAX_DIR,
             "\nPlease copy NCBI taxonomy files using: bash scripts/create_database/copy_ncbi_taxonomy_local.sh"))
}

cat("Using NCBI taxonomy directory:", NCBI_TAX_DIR, "\n")
nodes <- getnodes(taxdir = NCBI_TAX_DIR)
taxdir <- NCBI_TAX_DIR

main_df = rbindlist(struo_list, fill = TRUE) %>%
	filter(!is.na(ncbi_species_taxid))  %>%
	arrange(ncbi_species_taxid) %>%
	select(c("ncbi_species_taxid", "accession", "ncbi_organism_name",
					 "fasta_file_path", "source",
					 "is_published")) %>%
	distinct(accession, .keep_all = T)

cat("After combining and filtering:", nrow(main_df), "genomes\n")

main_df$rank = CHNOSZ::getrank(main_df$ncbi_species_taxid, taxdir, nodes = nodes)

# Optimized taxonomy reading: use data.table for fast filtering
unique_tax_ids <- unique(main_df$ncbi_species_taxid)
n_tax_ids <- length(unique_tax_ids)

cat("Reading NCBI taxonomy for", n_tax_ids, "unique tax IDs...\n")
cat("Using data.table for fast filtering (reading entire file but filtering in memory is faster)\n")

ranked_lineage_file <- file.path(taxdir, "rankedlineage.dmp")
if (!file.exists(ranked_lineage_file)) {
  stop(paste("rankedlineage.dmp not found at:", ranked_lineage_file))
}

# Read entire file with data.table (fast) and filter in R (much faster than awk for many tax IDs)
# rankedlineage.dmp format: tab-delimited with pipe separators, 20 columns total
# We only need columns 1, 3, 5, 7, 9, 11, 13, 15, 17, 19 (every other column starting from 1)
cat("Reading rankedlineage.dmp file (this may take 1-2 minutes for 315MB file)...\n")
ranked_lineage_all <- fread(ranked_lineage_file, sep = "\t", data.table = TRUE,
                            fill = TRUE, showProgress = TRUE)

# Select every other column (the actual data columns, skipping pipe separators)
# Columns are: 1=tax_id, 3=tax_name, 5=species, 7=genus, 9=family, 11=order, 13=class, 15=phylum, 17=kingdom, 19=superkingdom
col_indices <- seq(from = 1, to = min(20, ncol(ranked_lineage_all)), by = 2)
ranked_lineage_all <- ranked_lineage_all[, ..col_indices]
setnames(ranked_lineage_all, c("tax_id", "tax_name", "species", "genus", "family", 
                                "order", "class", "phylum", "kingdom", "superkingdom"))

cat("Filtering to", n_tax_ids, "unique tax IDs...\n")
# Filter to only needed tax IDs (data.table is very fast for this)
ncbi_taxonomy <- ranked_lineage_all[as.numeric(tax_id) %in% unique_tax_ids] %>%
  as_tibble() %>%
  mutate(across(everything(), ~trimws(.x))) %>%
  mutate(across(c(species, genus, family, order, class, phylum, kingdom, superkingdom), 
                ~ifelse(.x == "" | is.na(.x), NA_character_, .x))) %>%
  mutate(tax_id = as.numeric(tax_id))

cat("NCBI taxonomy filtering complete (", nrow(ncbi_taxonomy), "entries found)\n")

# Join taxonomy with main_df by tax_id (not using cbind which requires exact match)
main_df <- main_df %>% 
  arrange(ncbi_species_taxid) %>%
  mutate(ncbi_species_taxid = as.numeric(ncbi_species_taxid))

ncbi_taxonomy <- ncbi_taxonomy %>%
  arrange(tax_id) %>%
  mutate(kingdom = ifelse(kingdom == "", superkingdom, kingdom))

# Join by tax_id (left join to keep all genomes, even if taxonomy is missing)
main_df <- main_df %>%
  left_join(ncbi_taxonomy, by = c("ncbi_species_taxid" = "tax_id")) %>%
	mutate(species = ifelse(rank == "species", tax_name, species)) %>%
	mutate(phylum = ifelse(rank == "phylum", tax_name, phylum)) %>%
	mutate(class = ifelse(rank == "class", tax_name, class)) %>%
	mutate(genus = ifelse(rank == "genus", tax_name, genus)) %>%
	mutate(family = ifelse(rank == "family", tax_name, family)) %>%
	mutate(kingdom = ifelse(rank == "kingdom", tax_name, kingdom)) %>%
	mutate(order = ifelse(rank == "order", tax_name, order)) %>%
	mutate(ncbi_taxonomy =
				 	paste0("k__", kingdom,
				 				 ";p__", phylum,
				 				 ";c__",class,
				 				 ";o__",order,
				 				 ";f__",family,
				 				 ";g__", genus,
				 				 ";s__", species))
main_df$ncbi_organism_name <- main_df$tax_name

cat("Before Viridiplantae filter:", nrow(main_df), "genomes\n")
main_df <- main_df %>% filter(is.na(kingdom) | kingdom != "Viridiplantae")
cat("After Viridiplantae filter:", nrow(main_df), "genomes\n")
table(duplicated(main_df$accession))

table(main_df$kingdom)
table(main_df$source)
table(main_df$phylum)
table(is.na(main_df$fasta_file_path))
table(file.exists(main_df$fasta_file_path))
table(is.na(main_df$ncbi_species_taxid))

# Output file path
output_file <- file.path(STRUO2_INPUT_DIR, "soil_genome_db_struo.tsv")
cat("\nWriting combined database to:", output_file, "\n")
write_tsv(main_df, output_file)
cat("✓ Combined database written successfully!\n")
cat("  Total genomes:", nrow(main_df), "\n")
cat("  Unique accessions:", length(unique(main_df$accession)), "\n")



# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total genomes:", nrow(main_df), "\n")
cat("Unique accessions:", length(unique(main_df$accession)), "\n")
cat("\nGenomes by source:\n")
print(table(main_df$source))
cat("\nGenomes by superkingdom:\n")
print(table(main_df$superkingdom))

# Create summary plot
if (requireNamespace("ggplot2", quietly = TRUE)) {
  summary_plot <- main_df %>% 
    group_by(source, superkingdom) %>% 
    tally() %>%
    mutate(source = recode(source, 
                           "SPIRE_MAGs" = "SPIRE catalog",
                           "SPIRE" = "SPIRE catalog",
                           "SMAG" = "SMAG catalog",
                           "GTDB_207" = "GTDB r207",
                           "GTDB" = "GTDB r207",
                           "JGI_GOLD" = "JGI GOLD",
                           "RefSoil" = "RefSoil")) %>% 
    ungroup() %>% 
    ggplot(aes(x = reorder(source, -n), y = n, fill = superkingdom)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    coord_flip() +
    theme_bw(base_size = 20) +
    xlab("Genome source") + 
    scale_y_sqrt() +
    ggtitle("SoilMicrobeDB contents") +
    ylab("Genome count")
  
  # Save plot
  plot_file <- file.path(STRUO2_INPUT_DIR, "soil_genome_db_summary.png")
  ggsave(plot_file, summary_plot, width = 12, height = 8, dpi = 300)
  cat("\nSummary plot saved to:", plot_file, "\n")
} else {
  cat("\nggplot2 not available, skipping plot generation\n")
}



