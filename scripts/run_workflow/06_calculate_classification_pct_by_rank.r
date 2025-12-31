#!/usr/bin/env Rscript
# 06: Calculate classification percentage at each taxonomic rank
# Shows what % of total sequencing depth was classified to each specific rank
# Calculates before and after Architeuthis filtering
#
# Usage: Rscript scripts/run_workflow/06_calculate_classification_pct_by_rank.r
#
# Input:  *_kraken.kreport and *_filtered_kraken.kreport files from Steps 1-2
# Output: classification_pct_by_rank.csv (summary)
#         classification_pct_by_rank_per_sample.csv (detailed, preserved)
#
# Note: Calculates cumulative proportion (reads classified at this rank OR MORE SPECIFIC)
#       Uses cladeReads (cumulative, includes all children)
#       Denominator: total sequencing depth (seq_depth)
#       
#       For GTDB unfiltered: species and strain may have identical values because
#       all species-level classifications go to strain level. This is correct cumulative
#       behavior: species cumulative includes strain children, and if all reads are at
#       strain level, then species cumulative = strain cumulative.

library(tidyverse)
library(data.table)
library(future.apply)
library(pavian)
source("scripts/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")

# Set up parallel processing (will be adjusted later based on files to process)
n_workers <- min(18, availableCores() - 1)
cat("Available workers:", n_workers, "\n")

# Read sequencing depth for normalization
seq_depth_file <- "data/NEON_metagenome_classification/seq_depth_df.rds"
if(!file.exists(seq_depth_file)) {
    seq_depth_file <- "data/classification/analysis_files/seq_depth_df.rds"
}
if(!file.exists(seq_depth_file)) {
    stop("❌ MISSING FILE: seq_depth_df.rds not found!\n",
         "   Checked: data/NEON_metagenome_classification/seq_depth_df.rds\n",
         "   Checked: data/classification/analysis_files/seq_depth_df.rds\n",
         "   Please run 01_calculate_sequencing_depth.r first.")
}
seq_depth_df <- readRDS(seq_depth_file)

# Define fungal phyla for filtering
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")
fungal_pattern <- paste(fungal_phyla, collapse = "|")

# Build taxonomy_id -> is_fungi mapping from merged lineage files
# This is more reliable than parsing lineage from kreport files
cat("Building taxonomy_id -> is_fungi mapping from merged lineage files...\n")
lineage_dir_new <- "data/NEON_metagenome_classification/summary_files"
lineage_dir_old <- "data/classification/taxonomic_rank_summaries"
taxonomy_is_fungi <- data.frame(taxonomy_id = integer(), is_fungi = logical(), stringsAsFactors = FALSE)

databases_to_check <- c("soil_microbe_db", "pluspf", "gtdb_207", "gtdb_207_unfiltered")
for(db in databases_to_check) {
    # Check for species-level merged lineage file (most comprehensive)
    filename <- paste0(db, "_species_merged_lineage.csv")
    lineage_file <- NULL
    
    if(file.exists(file.path(lineage_dir_new, filename))) {
        lineage_file <- file.path(lineage_dir_new, filename)
    } else if(file.exists(file.path(lineage_dir_old, "species", filename))) {
        lineage_file <- file.path(lineage_dir_old, "species", filename)
    } else if(file.exists(file.path(lineage_dir_old, filename))) {
        lineage_file <- file.path(lineage_dir_old, filename)
    }
    
    if(!is.null(lineage_file) && file.exists(lineage_file)) {
        tryCatch({
            lineage_df <- read_csv(lineage_file, col_select = c("taxonomy_id", "lineage"), show_col_types = FALSE)
            if(nrow(lineage_df) > 0) {
                lineage_df <- lineage_df %>%
                    distinct(taxonomy_id, .keep_all = TRUE) %>%
                    mutate(is_fungi = grepl(fungal_pattern, lineage, ignore.case = TRUE)) %>%
                    select(taxonomy_id, is_fungi)
                taxonomy_is_fungi <- bind_rows(taxonomy_is_fungi, lineage_df)
                cat("  Loaded", nrow(lineage_df), "taxonomy IDs from", db, "\n")
            }
        }, error = function(e) {
            cat("  Warning: Could not read", lineage_file, ":", e$message, "\n")
        })
    }
}

# Also check genus and phylum files to catch any taxonomy_ids not in species files
for(rank in c("genus", "phylum")) {
    for(db in databases_to_check) {
        filename <- paste0(db, "_", rank, "_merged_lineage.csv")
        lineage_file <- NULL
        
        if(file.exists(file.path(lineage_dir_new, filename))) {
            lineage_file <- file.path(lineage_dir_new, filename)
        } else if(file.exists(file.path(lineage_dir_old, rank, filename))) {
            lineage_file <- file.path(lineage_dir_old, rank, filename)
        } else if(file.exists(file.path(lineage_dir_old, filename))) {
            lineage_file <- file.path(lineage_dir_old, filename)
        }
        
        if(!is.null(lineage_file) && file.exists(lineage_file)) {
            tryCatch({
                lineage_df <- read_csv(lineage_file, col_select = c("taxonomy_id", "lineage"), show_col_types = FALSE)
                if(nrow(lineage_df) > 0) {
                    lineage_df <- lineage_df %>%
                        distinct(taxonomy_id, .keep_all = TRUE) %>%
                        mutate(is_fungi = grepl(fungal_pattern, lineage, ignore.case = TRUE)) %>%
                        select(taxonomy_id, is_fungi)
                    # Only add taxonomy_ids we haven't seen yet
                    new_ids <- lineage_df %>%
                        filter(!taxonomy_id %in% taxonomy_is_fungi$taxonomy_id)
                    if(nrow(new_ids) > 0) {
                        taxonomy_is_fungi <- bind_rows(taxonomy_is_fungi, new_ids)
                    }
                }
            }, error = function(e) {
                # Silently skip if file doesn't exist or can't be read
            })
        }
    }
}

# Create a lookup table (data.table for fast joins)
if(nrow(taxonomy_is_fungi) > 0) {
    taxonomy_is_fungi <- taxonomy_is_fungi %>%
        group_by(taxonomy_id) %>%
        summarize(is_fungi = any(is_fungi), .groups = "drop") %>%
        arrange(taxonomy_id)
    setDT(taxonomy_is_fungi)
    cat("✓ Created taxonomy_id -> is_fungi mapping with", nrow(taxonomy_is_fungi), "unique taxonomy IDs\n")
    cat("  Fungal taxonomy IDs:", sum(taxonomy_is_fungi$is_fungi), "\n")
} else {
    cat("⚠ Warning: No taxonomy_id -> is_fungi mapping created. Will fall back to name-based identification.\n")
    taxonomy_is_fungi <- NULL
}

# Set up paths - check both new and old locations, plus HARDDRIVE
# Before filtering: original Kraken2 kreport files
kreport_before_local = "data/NEON_metagenome_classification/01_kraken_output"
kreport_before_local_old = "data/classification/01_kraken_output"
kreport_before_harddrive = "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/01_kraken_output"

# After filtering: filtered kreport files
kreport_after_local = "data/NEON_metagenome_classification/02_bracken_output"
kreport_after_local_old = "data/classification/02_bracken_output"
kreport_after_harddrive = "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/02_bracken_output"

# Find kreport files BEFORE filtering (original Kraken2 output)
kreport_before_files <- character(0)
if(dir.exists(kreport_before_local)) {
    kreport_before_files <- c(kreport_before_files, list.files(kreport_before_local, pattern = "_kraken.kreport", full.names = TRUE, recursive = TRUE))
}
if(dir.exists(kreport_before_local_old)) {
    kreport_before_files <- c(kreport_before_files, list.files(kreport_before_local_old, pattern = "_kraken.kreport", full.names = TRUE, recursive = TRUE))
    cat("✓ Found", length(list.files(kreport_before_local_old, pattern = "_kraken.kreport", full.names = TRUE, recursive = TRUE)), "kreport files (before filtering) in old location\n")
}
if(dir.exists(kreport_before_harddrive)) {
    kreport_before_harddrive_files <- list.files(kreport_before_harddrive, pattern = "_kraken.kreport", full.names = TRUE, recursive = TRUE)
    kreport_before_files <- c(kreport_before_files, kreport_before_harddrive_files)
    cat("✓ Found", length(kreport_before_harddrive_files), "kreport files (before filtering) on HARDDRIVE\n")
}

# Find kreport files AFTER filtering (Architeuthis filtered)
kreport_after_files <- character(0)
if(dir.exists(kreport_after_local)) {
    kreport_after_files <- c(kreport_after_files, list.files(kreport_after_local, pattern = "_filtered_kraken.kreport", full.names = TRUE, recursive = TRUE))
}
if(dir.exists(kreport_after_local_old)) {
    kreport_after_files <- c(kreport_after_files, list.files(kreport_after_local_old, pattern = "_filtered_kraken.kreport", full.names = TRUE, recursive = TRUE))
    cat("✓ Found", length(list.files(kreport_after_local_old, pattern = "_filtered_kraken.kreport", full.names = TRUE, recursive = TRUE)), "kreport files (after filtering) in old location\n")
}
if(dir.exists(kreport_after_harddrive)) {
    kreport_after_harddrive_files <- list.files(kreport_after_harddrive, pattern = "_filtered_kraken.kreport", full.names = TRUE, recursive = TRUE)
    kreport_after_files <- c(kreport_after_files, kreport_after_harddrive_files)
    cat("✓ Found", length(kreport_after_harddrive_files), "kreport files (after filtering) on HARDDRIVE\n")
}

cat("Found", length(kreport_before_files), "kreport files BEFORE filtering\n")
cat("Found", length(kreport_after_files), "kreport files AFTER filtering\n")

if(length(kreport_before_files) == 0 && length(kreport_after_files) == 0) {
    stop("No kreport files found! Checked:\n",
         "  Before (new): ", kreport_before_local, "\n",
         "  Before (old): ", kreport_before_local_old, "\n",
         "  Before (HARDDRIVE): ", kreport_before_harddrive, "\n",
         "  After (new): ", kreport_after_local, "\n",
         "  After (old): ", kreport_after_local_old, "\n",
         "  After (HARDDRIVE): ", kreport_after_harddrive, "\n")
}

# Function to parse sample ID from kreport filename
parse_sample_id <- function(filename, is_filtered = FALSE) {
    if(is_filtered) {
        samp_name <- gsub("_filtered_kraken.kreport", "", basename(filename))
    } else {
        samp_name <- gsub("_kraken.kreport", "", basename(filename))
    }
    
    # Parse sampleID and db_name
    parts <- strsplit(samp_name, "COMP_", fixed = TRUE)[[1]]
    if(length(parts) >= 2) {
        sampleID_temp <- parts[1]
        db_name <- parts[2]
        # Remove db_name suffix from sampleID
        sampleID <- sub(paste0("_", db_name), "", samp_name, fixed = TRUE)
        db_name <- gsub("_filtered|_kraken", "", db_name)
    } else {
        # Try alternative parsing for files without COMP_ separator
        if(grepl("soil_microbe_db", samp_name)) {
            sampleID <- gsub("_soil_microbe_db", "", samp_name)
            db_name <- "soil_microbe_db"
        } else {
            sampleID <- samp_name
            db_name <- "unknown"
        }
    }
    return(list(sampleID = sampleID, db_name = db_name, samp_name = samp_name))
}

# Optimized function: Read kreport file once and extract ALL ranks at once
read_kreport_all_ranks <- function(kreport_file, filter_status, taxonomy_is_fungi_map) {
    # Read kreport file - handle both 6-column and 8-column formats
    first_line <- readLines(kreport_file, n = 1)
    has_header <- grepl("^[a-zA-Z]", first_line)
    
    if(has_header) {
        report <- read_report3(kreport_file, keep_taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S", "S1"), has_header = TRUE)
    } else {
        # Check number of columns using fread for speed
        test_line <- fread(kreport_file, sep = "\t", nrows = 1, header = FALSE)
        n_cols <- ncol(test_line)
        
        if(n_cols == 8) {
            # 8-column format: use fread for speed
            report <- fread(kreport_file, sep = "\t", header = FALSE,
                          col.names = c("percentage", "cladeReads", "taxonReads", "nKmers", "n_unique_kmers", 
                                       "taxRank", "taxID", "name"),
                          quote = "", stringsAsFactors = FALSE)
            # Remove comment lines if any
            report <- report[!grepl("^#", report$name), ]
        } else if(n_cols == 6) {
            # 6-column format: use fread for speed
            report <- fread(kreport_file, sep = "\t", header = FALSE,
                          col.names = c("percentage", "cladeReads", "taxonReads", "taxRank", "taxID", "name"),
                          quote = "", stringsAsFactors = FALSE)
            # Remove comment lines if any
            report <- report[!grepl("^#", report$name), ]
        } else {
            # Try read_report3 as fallback
            report <- read_report3(kreport_file, keep_taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S", "S1"), has_header = FALSE)
        }
        
        # Process the report to match read_report3 format
        report$depth <- nchar(gsub("\\S.*", "", report$name)) / 2
        original_name <- gsub("^ *", "", report$name)
        report$name <- paste(tolower(report$taxRank), original_name, sep = "_")
        
        # Build taxLineage if it doesn't exist (needed for fungi identification at lower ranks)
        if(!"taxLineage" %in% names(report)) {
            report$taxLineage <- original_name
            rows_to_consider <- rep(FALSE, nrow(report))
            for(i in seq_len(nrow(report))) {
                if(i > 1 && report$depth[i] > 0) {
                    idx <- report$depth < report$depth[i] & rows_to_consider
                    if(any(idx)) {
                        my_row <- max(which(idx))
                        report$taxLineage[i] <- paste(report$taxLineage[my_row], 
                                                      report$taxLineage[i], sep = "|")
                    }
                }
                rows_to_consider[i] <- TRUE
            }
        }
    }
    
    if(is.null(report) || nrow(report) == 0) {
        return(NULL)
    }
    
    # Parse sample ID once
    is_filtered <- grepl("_filtered_kraken.kreport", kreport_file)
    sample_info <- parse_sample_id(kreport_file, is_filtered = is_filtered)
    
    # Rank mapping
    rank_map <- c("phylum" = "P", "class" = "C", "order" = "O", 
                  "family" = "F", "genus" = "G", "species" = "S", "strain" = "S1")
    
    # Extract all ranks at once using data.table for speed
    # Convert to data.table safely (don't modify in place if it's already a data.table)
    if(!inherits(report, "data.table")) {
        setDT(report)
    }
    
    # Identify fungi using taxonomy_id mapping (most reliable) or lineage/name (fallback)
    # Check if we have both the mapping and taxID column
    has_taxid_mapping <- !is.null(taxonomy_is_fungi_map) && nrow(taxonomy_is_fungi_map) > 0
    has_taxid_column <- "taxID" %in% names(report)
    
    if(has_taxid_mapping && has_taxid_column) {
        # Use taxonomy_id mapping from merged lineage files (most reliable)
        # Convert taxID to numeric if needed
        if(!is.numeric(report$taxID)) {
            report[, taxID := as.numeric(taxID)]
        }
        report <- merge(report, taxonomy_is_fungi_map, by.x = "taxID", by.y = "taxonomy_id", all.x = TRUE)
        report[, is_fungi := ifelse(is.na(is_fungi), FALSE, is_fungi)]
    } else {
        # Fallback: use taxLineage or name for fungi identification
        # Use global fungal_pattern defined at top of script
        if("taxLineage" %in% names(report)) {
            report[, is_fungi := grepl(fungal_pattern, taxLineage, ignore.case = TRUE)]
        } else {
            report[, is_fungi := grepl(fungal_pattern, name, ignore.case = TRUE)]
        }
    }
    
    # Process all ranks
    results_list <- list()
    
    for(rank_name in names(rank_map)) {
        rank_letter <- rank_map[rank_name]
        
        # Cumulative: sum cladeReads for this rank (includes all children/more specific ranks)
        # cladeReads = cumulative reads assigned at this rank OR MORE SPECIFIC
        # This gives us the proportion of reads classified at this rank or below
        all_domain_reads <- report[taxRank == rank_letter, sum(cladeReads, na.rm = TRUE)]
        fungi_reads <- report[taxRank == rank_letter & is_fungi == TRUE, sum(cladeReads, na.rm = TRUE)]
        
        # Store results
        if(all_domain_reads > 0) {
            results_list[[paste0(rank_name, "_all_domain")]] <- data.frame(
                sampleID = sample_info$sampleID,
                db_name = sample_info$db_name,
                samp_name = sample_info$samp_name,
                taxonomic_rank = rank_name,
                filter_status = filter_status,
                taxon_group = "all_domain",
                reads_classified = all_domain_reads,
                stringsAsFactors = FALSE
            )
        }
        
        if(fungi_reads > 0) {
            results_list[[paste0(rank_name, "_fungi")]] <- data.frame(
                sampleID = sample_info$sampleID,
                db_name = sample_info$db_name,
                samp_name = sample_info$samp_name,
                taxonomic_rank = rank_name,
                filter_status = filter_status,
                taxon_group = "fungi",
                reads_classified = fungi_reads,
                stringsAsFactors = FALSE
            )
        }
    }
    
    if(length(results_list) == 0) {
        return(NULL)
    }
    
    return(bind_rows(results_list))
}

# Set up output files - use old location if it exists, otherwise new location
if(dir.exists("data/classification/analysis_files")) {
    output_file <- "data/classification/analysis_files/classification_pct_by_rank_per_sample.csv"
    processed_log_file <- "data/classification/analysis_files/classification_pct_by_rank_processed_files.txt"
} else {
    output_file <- "data/NEON_metagenome_classification/summary_files/classification_pct_by_rank_per_sample.csv"
    processed_log_file <- "data/NEON_metagenome_classification/summary_files/classification_pct_by_rank_processed_files.txt"
}

# Create output directory if needed
output_dir <- dirname(output_file)
if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
}

# Load list of already processed files
processed_files <- character(0)
if(file.exists(processed_log_file)) {
    processed_files <- readLines(processed_log_file)
    cat("Found", length(processed_files), "already processed files\n")
}

# Combine all files with their filter status
all_files <- c(
    setNames(kreport_before_files, rep("before", length(kreport_before_files))),
    setNames(kreport_after_files, rep("after", length(kreport_after_files)))
)

# Identify new files to process
new_files <- all_files[!all_files %in% processed_files]

cat("Found", length(all_files), "total kreport files\n")
cat("Found", length(processed_files), "already processed files\n")
cat("Found", length(new_files), "new files to process\n")

# Load existing results if output file exists
existing_results <- NULL
if(file.exists(output_file)) {
    cat("Loading existing results from:", output_file, "\n")
    existing_results <- read_csv(output_file, show_col_types = FALSE)
    cat("Found", nrow(existing_results), "existing records\n")
    
    # Get list of files already in results (by sampleID, db_name, and filter_status)
    # Include db_name to distinguish between different databases for the same sample
    existing_combos <- existing_results %>%
        distinct(sampleID, db_name, filter_status) %>%
        mutate(file_key = paste0(sampleID, "_", db_name, "_", filter_status))
    
    # Create file keys for new files by parsing sampleID and db_name from filename
    file_keys_new <- sapply(seq_along(new_files), function(i) {
        kreport_file <- new_files[i]
        filter_status <- names(new_files)[i]
        sample_info <- parse_sample_id(kreport_file, is_filtered = (filter_status == "after"))
        paste0(sample_info$sampleID, "_", sample_info$db_name, "_", filter_status)
    })
    
    # Filter out files that are already in results
    new_files <- new_files[!file_keys_new %in% existing_combos$file_key]
    cat("After checking existing results:", length(new_files), "new files to process\n")
}

if(length(new_files) == 0) {
    cat("✓ All files have already been processed. Exiting.\n")
    plan(sequential)
} else {
    # Process new files in parallel
cat("Processing", length(new_files), "new kreport files in parallel...\n")

# Adjust number of workers based on files to process
n_workers_actual <- min(n_workers, length(new_files))
cat("Using", n_workers_actual, "parallel workers\n")

# Set up parallel plan
plan(multisession, workers = n_workers_actual)

# Process files in parallel
results <- future_lapply(seq_along(new_files), function(i) {
    kreport_file <- new_files[i]
    filter_status <- names(new_files)[i]
    
    tryCatch({
        read_kreport_all_ranks(kreport_file, filter_status, taxonomy_is_fungi)
    }, error = function(e) {
        cat("    Error processing", basename(kreport_file), ":", e$message, "\n")
        return(NULL)
    })
}, future.seed = TRUE)

# Clean up parallel workers
plan(sequential)

# Combine new results
new_results <- bind_rows(Filter(Negate(is.null), results))

# Append to existing results if present
if(!is.null(existing_results) && nrow(existing_results) > 0) {
    all_results <- bind_rows(existing_results, new_results)
    cat("Appended", nrow(new_results), "new records to", nrow(existing_results), "existing records\n")
} else {
    all_results <- new_results
    cat("Created", nrow(all_results), "new records\n")
}

# Merge with sequencing depth to calculate percentage
# Sequencing depth is sample-specific, not database-specific
if(nrow(all_results) > 0) {
    # Remove any existing seq_depth columns (from previous runs or conflicts)
    all_results <- all_results %>%
        select(-any_of(c("seq_depth", "seq_depth.x", "seq_depth.y", "pct_classified")))
    
    # Join with sequencing depth
    all_results <- all_results %>%
        left_join(seq_depth_df %>% select(sampleID, seq_depth), by = "sampleID") %>%
        mutate(
            reads_classified = replace_na(reads_classified, 0),
            pct_classified = (reads_classified / seq_depth) * 100
        )
    
    # Report missing seq_depth
    missing_seq_depth <- sum(is.na(all_results$seq_depth))
    if(missing_seq_depth > 0) {
        cat("⚠ WARNING:", missing_seq_depth, "records missing seq_depth after join\n")
        cat("   Missing sampleIDs:", paste(unique(all_results$sampleID[is.na(all_results$seq_depth)]), collapse = ", "), "\n")
    }
} else {
    stop("No results generated! Check that kreport files are readable and contain data.")
}

# ============================================================================
# Aggregate results
# ============================================================================

all_domain_results <- all_results %>% filter(taxon_group == "all_domain")
fungi_results <- all_results %>% filter(taxon_group == "fungi")

all_domain_summary <- all_domain_results %>%
    group_by(db_name, taxonomic_rank, filter_status) %>%
    summarize(
        mean_pct_classified = mean(pct_classified, na.rm = TRUE),
        median_pct_classified = median(pct_classified, na.rm = TRUE),
        sd_pct_classified = sd(pct_classified, na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
    ) %>%
    mutate(taxon_group = "all_domain") %>%
    select(taxon_group, db_name, taxonomic_rank, filter_status, mean_pct_classified, median_pct_classified, sd_pct_classified, n_samples)

fungi_summary <- fungi_results %>%
    group_by(db_name, taxonomic_rank, filter_status) %>%
    summarize(
        mean_pct_classified = mean(pct_classified, na.rm = TRUE),
        median_pct_classified = median(pct_classified, na.rm = TRUE),
        sd_pct_classified = sd(pct_classified, na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
    ) %>%
    mutate(taxon_group = "fungi") %>%
    select(taxon_group, db_name, taxonomic_rank, filter_status, mean_pct_classified, median_pct_classified, sd_pct_classified, n_samples)

# ============================================================================
# Combine and save results
# ============================================================================

classification_by_rank <- bind_rows(all_domain_summary, fungi_summary) %>%
    arrange(taxon_group, db_name, taxonomic_rank, filter_status)

# Save summary results - use old location if it exists, otherwise new location
if(dir.exists("data/classification/analysis_files")) {
    summary_output <- "data/classification/analysis_files/classification_pct_by_rank.csv"
} else {
    summary_output <- "data/NEON_metagenome_classification/summary_files/classification_pct_by_rank.csv"
    if(!dir.exists(dirname(summary_output))) {
        dir.create(dirname(summary_output), recursive = TRUE)
    }
}
write_csv(classification_by_rank, summary_output)
cat("✅ Saved classification percentages by rank to:", summary_output, "\n")
cat("   Note: Percentages are % of total sequencing depth (seq_depth)\n")
cat("   Values are CUMULATIVE (reads classified at this rank OR HIGHER/more specific)\n")
cat("   filter_status indicates 'before' or 'after' Architeuthis filtering\n")

# Also save per-sample results
output_dir <- dirname(output_file)
if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
}
write_csv(all_results, output_file)
cat("✅ Saved per-sample classification percentages by rank to:", output_file, "\n")
cat("   Total records:", nrow(all_results), "\n")

# Update processed files log
processed_files <- c(processed_files, new_files)
writeLines(processed_files, processed_log_file)
cat("✅ Updated processed files log:", processed_log_file, "\n")
cat("   Total processed files:", length(processed_files), "\n")

# Clean up parallel workers
plan(sequential)
cat("✅ Script completed successfully\n")
}
