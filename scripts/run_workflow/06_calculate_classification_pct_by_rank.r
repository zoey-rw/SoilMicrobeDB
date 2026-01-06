#!/usr/bin/env Rscript
# 06: Calculate classification percentage at each taxonomic rank
# Shows what % of total sequencing depth was classified to each specific rank
# Calculates before and after Architeuthis filtering
#
# Usage: Rscript scripts/run_workflow/06_calculate_classification_pct_by_rank.r
#
# Input:  For "before": *_kraken.kreport files from Step 1
#         For "after": merged lineage files from Step 3 (*_*_merged_lineage.csv)
# Output: classification_pct_by_rank.csv (summary)
#         classification_pct_by_rank_per_sample.csv (detailed, preserved)
#
# Note: Calculates cumulative proportion (reads classified at this rank OR MORE SPECIFIC)
#       For "after" filtering: Uses merged lineage files (simplified approach)
#       For "before" filtering: Uses kreport files (no merged files exist)
#       Denominator: total sequencing depth (seq_depth)

library(tidyverse)
library(data.table)
library(future.apply)
library(pavian)
source("scripts/helper_functions.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")

# TEST MODE: Set to a sample ID or vector of sample IDs to test (e.g., "ABBY_001-M-20170607-COMP")
# Set to NULL to process all samples
TEST_SAMPLE_ID <- NULL

# Set up parallel processing
n_workers <- min(18, availableCores() - 1)
cat("Available workers:", n_workers, "\n")

# Read sequencing depth for normalization
# Check both locations and use the one with more samples (more complete)
seq_depth_file_new <- "data/NEON_metagenome_classification/seq_depth_df.rds"
seq_depth_file_old <- "data/classification/analysis_files/seq_depth_df.rds"

seq_depth_file <- NULL
if(file.exists(seq_depth_file_new) && file.exists(seq_depth_file_old)) {
    # Both exist - use the one with more samples
    df_new <- readRDS(seq_depth_file_new)
    df_old <- readRDS(seq_depth_file_old)
    if(nrow(df_old) > nrow(df_new)) {
        seq_depth_file <- seq_depth_file_old
        cat("Using", seq_depth_file_old, "(", nrow(df_old), "samples)\n")
    } else {
        seq_depth_file <- seq_depth_file_new
        cat("Using", seq_depth_file_new, "(", nrow(df_new), "samples)\n")
    }
} else if(file.exists(seq_depth_file_old)) {
    seq_depth_file <- seq_depth_file_old
} else if(file.exists(seq_depth_file_new)) {
    seq_depth_file <- seq_depth_file_new
} else {
    stop("❌ MISSING FILE: seq_depth_df.rds not found!\n",
         "   Checked: data/NEON_metagenome_classification/seq_depth_df.rds\n",
         "   Checked: data/classification/analysis_files/seq_depth_df.rds\n",
         "   Please run 01_calculate_sequencing_depth.r first.")
}
seq_depth_df <- readRDS(seq_depth_file)
# Ensure sampleID is character and trim whitespace
if(!"sampleID" %in% names(seq_depth_df)) {
    stop("seq_depth_df does not have 'sampleID' column. Columns:", paste(names(seq_depth_df), collapse = ", "))
}
seq_depth_df$sampleID <- as.character(trimws(seq_depth_df$sampleID))

# Debug in test mode
if(!is.null(TEST_SAMPLE_ID)) {
    cat("seq_depth_df loaded:", nrow(seq_depth_df), "samples\n")
    if(length(TEST_SAMPLE_ID) == 1) {
        cat("Test sample in seq_depth_df:", TEST_SAMPLE_ID %in% seq_depth_df$sampleID, "\n")
        if(TEST_SAMPLE_ID %in% seq_depth_df$sampleID) {
            cat("Test sample seq_depth:", seq_depth_df$seq_depth[seq_depth_df$sampleID == TEST_SAMPLE_ID], "\n")
        }
    } else {
        cat("Test samples in seq_depth_df:", sum(TEST_SAMPLE_ID %in% seq_depth_df$sampleID), "out of", length(TEST_SAMPLE_ID), "\n")
    }
}

# Define fungal phyla for filtering
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")
fungal_pattern <- paste(fungal_phyla, collapse = "|")

# Rank hierarchy for cumulative calculations (most specific to least specific)
rank_hierarchy <- c("strain" = "S1", "species" = "S", "genus" = "G", "family" = "F", 
                    "order" = "O", "class" = "C", "phylum" = "P")
rank_order <- c("strain", "species", "genus", "family", "order", "class", "phylum")

# ============================================================================
# Process "AFTER" filtering using merged lineage files
# ============================================================================

cat("Processing 'after' filtering data from merged lineage files...\n")
if(!is.null(TEST_SAMPLE_ID)) {
    cat("⚠ TEST MODE: Processing only sample:", TEST_SAMPLE_ID, "\n")
}

# Define paths once
lineage_dir_new <- "data/NEON_metagenome_classification/summary_files"
lineage_dir_old <- "data/classification/taxonomic_rank_summaries"
# IMPORTANT: Order matters! Put longer/more specific patterns first
# "gtdb_207_unfiltered" must come before "gtdb_207" to match correctly
databases_to_check <- c("soil_microbe_db", "pluspf", "gtdb_207_unfiltered", "gtdb_207")

# Helper function to find file path
find_lineage_file <- function(db, rank) {
    filename <- paste0(db, "_", rank, "_merged_lineage.csv")
    paths <- c(
        file.path(lineage_dir_new, filename),
        file.path(lineage_dir_old, rank, filename),
        file.path(lineage_dir_old, filename)
    )
    existing <- paths[file.exists(paths)]
    if(length(existing) > 0) return(existing[1]) else return(NULL)
}

# Vectorized file discovery
cat("Discovering merged lineage files...\n")
file_grid <- expand.grid(db = databases_to_check, rank = rank_order, stringsAsFactors = FALSE)
file_grid$path <- mapply(find_lineage_file, file_grid$db, file_grid$rank, SIMPLIFY = FALSE)
# Filter out NULL paths
file_grid <- file_grid[!sapply(file_grid$path, is.null), ]
file_grid$path <- unlist(file_grid$path)

cat("Found", nrow(file_grid), "merged lineage files to process\n")

# Read all files efficiently
cat("Reading merged lineage files...\n")
all_rank_data_list <- lapply(seq_len(nrow(file_grid)), function(i) {
    db <- file_grid$db[i]
    rank <- file_grid$rank[i]
    path <- file_grid$path[i]
    
    tryCatch({
        rank_data <- fread(path, 
                          select = c("sample_id", "name", "taxonomy_id", "taxonomy_lvl", 
                                    "new_est_reads", "lineage"),
                          colClasses = list(integer = "taxonomy_id", 
                                           character = c("sample_id", "name", "taxonomy_lvl", "lineage"),
                                           numeric = "new_est_reads"))
        
            if(nrow(rank_data) > 0) {
                rank_letter <- rank_hierarchy[rank]
                rank_data <- rank_data[taxonomy_lvl == rank_letter]
                
                # Filter to test sample if in test mode (before adding columns for efficiency)
                if(!is.null(TEST_SAMPLE_ID) && nrow(rank_data) > 0) {
                    rank_data <- rank_data[grepl(TEST_SAMPLE_ID, sample_id, fixed = TRUE)]
                }
                
                if(nrow(rank_data) > 0) {
                    rank_data[, `:=`(rank_name = rank, db_source = db)]
                    cat("  Loaded", nrow(rank_data), "records from", db, rank, "rank\n")
                    return(rank_data)
                }
            }
        return(NULL)
    }, error = function(e) {
        path_display <- if(!is.null(path) && is.character(path)) basename(path) else "unknown file"
        cat("  Warning: Could not read", path_display, ":", conditionMessage(e), "\n")
        return(NULL)
    })
})

# Remove NULL entries and combine
all_rank_data_list <- Filter(Negate(is.null), all_rank_data_list)
if(length(all_rank_data_list) == 0) {
    stop("No merged lineage files could be read!")
}

cat("Combining", length(all_rank_data_list), "rank files...\n")
combined_data <- rbindlist(all_rank_data_list, fill = TRUE)
cat("Total records:", nrow(combined_data), "\n")

# Filter to test sample if in test mode
if(!is.null(TEST_SAMPLE_ID)) {
    # Filter sample_id that contains the test sample ID
    combined_data <- combined_data[grepl(TEST_SAMPLE_ID, sample_id, fixed = TRUE)]
    cat("After filtering to test sample:", nrow(combined_data), "records\n")
}

# Optimized regex: single pass for sampleID parsing
cat("Parsing sampleIDs and identifying fungi...\n")

# Remove _filtered suffix
combined_data[, sampleID_temp := sub("_{1,2}filtered$", "", sample_id, fixed = FALSE)]

# Simplified parsing: remove database and rank suffixes
db_pattern <- paste(databases_to_check, collapse = "|")
rank_pattern <- paste(rank_order, collapse = "|")

# Initialize
combined_data[, sampleID := sampleID_temp]

# For -COMP patterns: remove everything after -COMP_ except keep -COMP
# Pattern: -COMP_db_rank or -COMP_db -> keep -COMP
# Check if sampleID contains -COMP_ (not just starts with it)
comp_rows <- grepl("-COMP_", combined_data$sampleID_temp, fixed = TRUE)
if(any(comp_rows)) {
    combined_data[comp_rows, sampleID := sub("-COMP_.*$", "-COMP", sampleID_temp)]
}

# For non-COMP patterns: remove _db_rank or _db
# Only apply to rows that don't have -COMP in the original (sampleID_temp)
# and haven't been processed yet (sampleID still equals sampleID_temp)
non_comp_rows <- !grepl("-COMP", combined_data$sampleID_temp, fixed = TRUE) & 
                 combined_data$sampleID == combined_data$sampleID_temp
if(any(non_comp_rows)) {
    # Pattern: _db_rank or __db_rank -> remove
    combined_data[non_comp_rows, sampleID := sub(paste0("_{1,2}(", db_pattern, ")_{1,2}(", rank_pattern, ")$"), "", sampleID_temp)]
    # Pattern: _db or __db -> remove (only if still needed)
    still_needs_processing <- non_comp_rows & combined_data$sampleID == combined_data$sampleID_temp
    if(any(still_needs_processing)) {
        combined_data[still_needs_processing, sampleID := sub(paste0("_{1,2}(", db_pattern, ")$"), "", sampleID_temp)]
    }
}

# Clean trailing underscores
combined_data[, sampleID := sub("_+$", "", sampleID, fixed = FALSE)]
combined_data[, sampleID_temp := NULL]

# Debug output in test mode
if(!is.null(TEST_SAMPLE_ID)) {
    cat("\nSampleID parsing test (all unique combinations):\n")
    test_df <- combined_data[, .(sample_id, sampleID)] %>% unique()
    print(test_df)
    cat("\nChecking if sampleIDs match seq_depth_df:\n")
    seq_df <- readRDS(seq_depth_file)
    matches <- sum(unique(combined_data$sampleID) %in% seq_df$sampleID)
    cat("Matching sampleIDs:", matches, "out of", length(unique(combined_data$sampleID)), "\n")
}

# Extract db_name efficiently (single regex pass)
db_pattern_regex <- paste0("(", db_pattern, ")")
combined_data[, db_name := stringr::str_extract(sample_id, db_pattern_regex)]
combined_data[is.na(db_name), db_name := db_source]

# Identify fungi (vectorized)
combined_data[, is_fungi := grepl(fungal_pattern, lineage, ignore.case = TRUE, perl = TRUE)]

# Optimized cumulative calculation using data.table grouping
cat("Calculating cumulative reads by rank...\n")

# Convert rank_name to factor with hierarchy order for proper sorting
# rank_order is from most specific to least specific
combined_data[, rank_name := factor(rank_name, levels = rank_order)]

# Calculate reads per rank per sample
# IMPORTANT: Each rank file (species, genus, phylum, etc.) already contains
# reads classified at that level OR MORE SPECIFIC. Therefore, we do NOT need
# to calculate cumulative sums - we use the reads directly from each rank file.
# Group by sampleID, db_name, and rank_name only (not sample_id, which varies by file)
reads_by_rank <- combined_data[, .(
    reads_classified = sum(new_est_reads, na.rm = TRUE),
    fungi_reads = sum(new_est_reads[is_fungi == TRUE], na.rm = TRUE)
), by = .(sampleID, db_name, rank_name)]

# NO CUMULATIVE SUM NEEDED - each rank file already represents reads at that level or below
# The merged_lineage files are already aggregated, so we use reads_classified directly

# Reshape to long format for output
after_results_list <- list()

# All domain results
# Use parsed sampleID (not sample_id) for samp_name to ensure consistency
all_domain_results <- reads_by_rank[reads_classified > 0, .(
    sampleID, db_name, samp_name = sampleID,
    taxonomic_rank = as.character(rank_name), filter_status = "after",
    taxon_group = "all_domain", reads_classified = reads_classified
)]
if(nrow(all_domain_results) > 0) {
    after_results_list[["all_domain"]] <- all_domain_results
}

# Fungi results
fungi_results <- reads_by_rank[fungi_reads > 0, .(
    sampleID, db_name, samp_name = sampleID,
    taxonomic_rank = as.character(rank_name), filter_status = "after",
    taxon_group = "fungi", reads_classified = fungi_reads
)]
if(nrow(fungi_results) > 0) {
    after_results_list[["fungi"]] <- fungi_results
}

after_results <- if(length(after_results_list) > 0) {
    result <- rbindlist(after_results_list, fill = TRUE)
    # Convert to data.frame for compatibility with dplyr joins
    result <- as.data.frame(result)
    # Ensure sampleID is character (not factor) and trim whitespace
    result$sampleID <- as.character(trimws(result$sampleID))
    cat("✓ Processed", nrow(result), "records from merged lineage files (after filtering)\n")
    
    # Debug in test mode
    if(!is.null(TEST_SAMPLE_ID)) {
        cat("  SampleIDs in after_results:", paste(unique(result$sampleID), collapse = ", "), "\n")
        cat("  SampleID types:", class(result$sampleID), "\n")
    }
    
    result
} else {
    cat("⚠ No 'after' filtering results generated\n")
    data.frame(sampleID = character(), db_name = character(), samp_name = character(),
               taxonomic_rank = character(), filter_status = character(),
               taxon_group = character(), reads_classified = numeric(),
               stringsAsFactors = FALSE)
}

# ============================================================================
# Process "BEFORE" filtering using kreport files
# ============================================================================

cat("\nProcessing 'before' filtering data from kreport files...\n")

# Build taxonomy_id -> is_fungi mapping from merged lineage files (reuse from above)
taxonomy_is_fungi <- NULL
if(nrow(after_results) > 0) {
    # Extract unique taxonomy_id -> is_fungi mapping from the combined data we already loaded
    # We'll rebuild it from the merged files
    taxonomy_is_fungi_list <- list()
    
    for(db in databases_to_check) {
        for(rank in c("species", "genus", "phylum")) {
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
                    lineage_df <- read_csv(lineage_file, col_select = c("taxonomy_id", "lineage"), 
                                          col_types = cols(taxonomy_id = col_integer(), lineage = col_character()),
                                          show_col_types = FALSE)
                    if(nrow(lineage_df) > 0) {
                        lineage_df <- lineage_df %>%
                            mutate(taxonomy_id = as.integer(round(taxonomy_id))) %>%
                            filter(!is.na(taxonomy_id)) %>%
                            distinct(taxonomy_id, .keep_all = TRUE) %>%
                            mutate(is_fungi = grepl(fungal_pattern, lineage, ignore.case = TRUE)) %>%
                            select(taxonomy_id, is_fungi)
                        taxonomy_is_fungi_list[[paste0(db, "_", rank)]] <- lineage_df
                    }
                }, error = function(e) {
                    # Silently skip
                })
            }
        }
    }
    
    if(length(taxonomy_is_fungi_list) > 0) {
        taxonomy_is_fungi <- bind_rows(taxonomy_is_fungi_list) %>%
            mutate(taxonomy_id = as.integer(round(taxonomy_id))) %>%
            filter(!is.na(taxonomy_id)) %>%
            group_by(taxonomy_id) %>%
            summarize(is_fungi = any(is_fungi), .groups = "drop") %>%
            arrange(taxonomy_id)
        setDT(taxonomy_is_fungi)
        taxonomy_is_fungi[, taxonomy_id := as.integer(taxonomy_id)]
        cat("✓ Created taxonomy_id -> is_fungi mapping with", nrow(taxonomy_is_fungi), "unique taxonomy IDs\n")
    }
}

# Function to parse sample ID from kreport filename
parse_sample_id <- function(filename, is_filtered = FALSE) {
    if(is_filtered) {
        samp_name <- gsub("_{1,2}filtered_kraken\\.kreport$", "", basename(filename))
    } else {
        samp_name <- gsub("_kraken\\.kreport$", "", basename(filename))
    }
    
    parts <- strsplit(samp_name, "COMP_", fixed = TRUE)[[1]]
    if(length(parts) >= 2) {
        sampleID_temp <- parts[1]
        db_name <- parts[2]
        sampleID <- sub(paste0("__", db_name), "", samp_name, fixed = TRUE)
        if(sampleID == samp_name) {
            sampleID <- sub(paste0("_", db_name), "", samp_name, fixed = TRUE)
        }
        db_name <- gsub("_{1,2}filtered|_kraken", "", db_name)
    } else {
        if(grepl("soil_microbe_db", samp_name)) {
            sampleID <- gsub("_{1,2}soil_microbe_db", "", samp_name)
            db_name <- "soil_microbe_db"
        } else {
            sampleID <- samp_name
            db_name <- "unknown"
        }
    }
    # Clean up sampleID: remove _filtered, database, and rank suffixes
    sampleID <- sub("_{1,2}filtered$", "", sampleID, fixed = FALSE)
    db_pattern <- paste(c("gtdb_207", "gtdb_207_unfiltered", "pluspf", "soil_microbe_db"), collapse = "|")
    rank_pattern <- paste(c("species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"), collapse = "|")
    # Remove _db_rank or _db patterns
    sampleID <- sub(paste0("_{1,2}(", db_pattern, ")_{1,2}(", rank_pattern, ")$"), "", sampleID)
    sampleID <- sub(paste0("_{1,2}(", db_pattern, ")$"), "", sampleID)
    # Clean trailing underscores
    sampleID <- sub("_+$", "", sampleID, fixed = FALSE)
    # Use parsed sampleID as samp_name for consistency
    return(list(sampleID = sampleID, db_name = db_name, samp_name = sampleID))
}

# Function to read kreport file and extract all ranks (optimized)
read_kreport_all_ranks <- function(kreport_file, filter_status, taxonomy_is_fungi_map) {
    first_line <- readLines(kreport_file, n = 1)
    has_header <- grepl("^[a-zA-Z]", first_line)
    
    if(has_header) {
        report <- read_report3(kreport_file, keep_taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S", "S1"), has_header = TRUE)
    } else {
        test_line <- fread(kreport_file, sep = "\t", nrows = 1, header = FALSE)
        n_cols <- ncol(test_line)
        
        if(n_cols == 8) {
            report <- fread(kreport_file, sep = "\t", header = FALSE,
                          col.names = c("percentage", "cladeReads", "taxonReads", "nKmers", "n_unique_kmers", 
                                       "taxRank", "taxID", "name"),
                          quote = "", stringsAsFactors = FALSE)
            report <- report[!grepl("^#", report$name), ]
        } else if(n_cols == 6) {
            report <- fread(kreport_file, sep = "\t", header = FALSE,
                          col.names = c("percentage", "cladeReads", "taxonReads", "taxRank", "taxID", "name"),
                          quote = "", stringsAsFactors = FALSE)
            report <- report[!grepl("^#", report$name), ]
        } else {
            report <- read_report3(kreport_file, keep_taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S", "S1"), has_header = FALSE)
        }
        
        report$depth <- nchar(gsub("\\S.*", "", report$name)) / 2
        original_name <- gsub("^ *", "", report$name)
        report$name <- paste(tolower(report$taxRank), original_name, sep = "_")
        
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
    
    is_filtered <- grepl("_filtered_kraken.kreport", kreport_file)
    sample_info <- parse_sample_id(kreport_file, is_filtered = is_filtered)
    
    rank_map <- c("phylum" = "P", "class" = "C", "order" = "O", 
                  "family" = "F", "genus" = "G", "species" = "S", "strain" = "S1")
    
    if(!inherits(report, "data.table")) {
        setDT(report)
    }
    
    # Optimization: Filter to relevant ranks BEFORE merge (memory efficient)
    target_ranks <- rank_map
    report <- report[taxRank %in% target_ranks]
    
    if(nrow(report) == 0) {
        return(NULL)
    }
    
    has_taxid_mapping <- !is.null(taxonomy_is_fungi_map) && nrow(taxonomy_is_fungi_map) > 0
    has_taxid_column <- "taxID" %in% names(report)
    
    if(has_taxid_mapping && has_taxid_column) {
        # Optimization: Use setkey for fast merge
        if(!is.numeric(report$taxID)) {
            report[, taxID := as.numeric(taxID)]
        }
        report[, taxID := as.integer(taxID)]
        taxonomy_is_fungi_map[, taxonomy_id := as.integer(taxonomy_id)]
        setkey(taxonomy_is_fungi_map, taxonomy_id)
        setkey(report, taxID)
        report <- taxonomy_is_fungi_map[report]
        report[, is_fungi := fifelse(is.na(is_fungi), FALSE, is_fungi)]
    } else {
        if("taxLineage" %in% names(report)) {
            report[, is_fungi := grepl(fungal_pattern, taxLineage, ignore.case = TRUE, perl = TRUE)]
        } else {
            report[, is_fungi := grepl(fungal_pattern, name, ignore.case = TRUE, perl = TRUE)]
        }
    }
    
    # Optimization: Calculate all ranks in one pass using data.table grouping
    rank_map_dt <- data.table(rank_name = names(rank_map), taxRank = rank_map)
    setkey(rank_map_dt, taxRank)
    setkey(report, taxRank)
    report <- rank_map_dt[report]
    
    # Aggregate by rank in one pass
    results_dt <- report[, .(
        all_domain_reads = sum(cladeReads, na.rm = TRUE),
        fungi_reads = sum(cladeReads[is_fungi == TRUE], na.rm = TRUE)
    ), by = .(rank_name)]
    
    # Build results list
    results_list <- list()
    for(i in seq_len(nrow(results_dt))) {
        rank_name <- results_dt$rank_name[i]
        all_domain_reads <- results_dt$all_domain_reads[i]
        fungi_reads <- results_dt$fungi_reads[i]
        
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

# Find kreport files BEFORE filtering
kreport_before_local = "data/NEON_metagenome_classification/01_kraken_output"
kreport_before_local_old = "data/classification/01_kraken_output"
kreport_before_harddrive = "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/01_kraken_output"

kreport_before_files <- character(0)
if(dir.exists(kreport_before_local)) {
    kreport_before_files <- c(kreport_before_files, list.files(kreport_before_local, pattern = "_kraken.kreport", full.names = TRUE, recursive = TRUE))
}
if(dir.exists(kreport_before_local_old)) {
    kreport_before_files <- c(kreport_before_files, list.files(kreport_before_local_old, pattern = "_kraken.kreport", full.names = TRUE, recursive = TRUE))
}
if(dir.exists(kreport_before_harddrive)) {
    kreport_before_files <- c(kreport_before_files, list.files(kreport_before_harddrive, pattern = "_kraken.kreport", full.names = TRUE, recursive = TRUE))
}

cat("Found", length(kreport_before_files), "kreport files BEFORE filtering\n")

# Set up output files (consolidate path management)
# Use NEON_metagenome_classification directory structure
base_dir <- if(dir.exists("data/NEON_metagenome_classification/analysis_files")) {
    "data/NEON_metagenome_classification/analysis_files"
} else if(dir.exists("data/NEON_metagenome_classification")) {
    "data/NEON_metagenome_classification"
} else if(dir.exists("data/classification/analysis_files")) {
    "data/classification/analysis_files"
} else {
    "data/classification/analysis_files"
}
output_file <- file.path(base_dir, "classification_pct_by_rank_per_sample.csv")
processed_log_file <- file.path(base_dir, "classification_pct_by_rank_processed_files.txt")
summary_output <- file.path(base_dir, "classification_pct_by_rank.csv")

output_dir <- dirname(output_file)
if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Load existing results
existing_results <- NULL
if(file.exists(output_file)) {
    existing_results <- read_csv(output_file, show_col_types = FALSE)
    cat("Found", nrow(existing_results), "existing records\n")
}

# Process "before" filtering files
before_results <- NULL
if(length(kreport_before_files) > 0) {
    # Check which files need processing
    processed_files <- character(0)
    if(file.exists(processed_log_file)) {
        processed_files <- readLines(processed_log_file)
    }
    
    # Get existing combos from results
    existing_combos <- NULL
    if(!is.null(existing_results) && nrow(existing_results) > 0) {
        existing_combos <- existing_results %>%
            filter(filter_status == "before") %>%
            distinct(sampleID, db_name) %>%
            mutate(file_key = paste0(sampleID, "_", db_name))
    }
    
    # Identify files to process
    files_to_process <- kreport_before_files
    if(!is.null(existing_combos)) {
        file_keys_new <- sapply(files_to_process, function(f) {
            sample_info <- parse_sample_id(f, is_filtered = FALSE)
            paste0(sample_info$sampleID, "_", sample_info$db_name)
        })
        files_to_process <- files_to_process[!file_keys_new %in% existing_combos$file_key]
    }
    
    if(length(files_to_process) > 0) {
        cat("Processing", length(files_to_process), "kreport files (before filtering)...\n")
        n_workers_actual <- min(n_workers, length(files_to_process))
        plan(multisession, workers = n_workers_actual)
        
        before_results_list <- future_lapply(files_to_process, function(kreport_file) {
            tryCatch({
                read_kreport_all_ranks(kreport_file, "before", taxonomy_is_fungi)
            }, error = function(e) {
                cat("    Error processing", basename(kreport_file), ":", e$message, "\n")
                return(NULL)
            })
        }, future.seed = TRUE)
        
        plan(sequential)
        before_results <- bind_rows(Filter(Negate(is.null), before_results_list))
        cat("✓ Processed", nrow(before_results), "records from kreport files (before filtering)\n")
    } else {
        cat("✓ All 'before' filtering files already processed\n")
    }
}

# Combine results
# Always replace "after" filtering results with new ones (from merged files)
# Keep existing "before" filtering results unless new ones are available
if(!is.null(existing_results) && nrow(existing_results) > 0) {
    # Keep only "before" filtering results from existing file (exclude old "after" results)
    existing_before <- existing_results %>% filter(filter_status == "before")
    
    # Combine: existing before + new before (if any) + new after (replace old after)
    result_list <- list()
    if(nrow(existing_before) > 0) {
        result_list[["existing_before"]] <- existing_before
    }
    if(!is.null(before_results) && nrow(before_results) > 0) {
        result_list[["new_before"]] <- before_results
    }
    if(nrow(after_results) > 0) {
        result_list[["new_after"]] <- after_results
    }
    
    all_results <- if(length(result_list) > 0) {
        bind_rows(result_list)
    } else {
        data.frame()
    }
} else {
    result_list <- list()
    if(!is.null(before_results) && nrow(before_results) > 0) {
        result_list[["before"]] <- before_results
    }
    if(nrow(after_results) > 0) {
        result_list[["after"]] <- after_results
    }
    all_results <- if(length(result_list) > 0) {
        bind_rows(result_list)
    } else {
        data.frame()
    }
}

# Merge with sequencing depth
if(nrow(all_results) > 0) {
    # Debug in test mode
    if(!is.null(TEST_SAMPLE_ID)) {
        cat("\nBefore join:\n")
        after_in_results <- all_results %>% filter(filter_status == "after")
        cat("  After filtering sampleIDs:", paste(unique(after_in_results$sampleID), collapse = ", "), "\n")
        cat("  After filtering sampleID type:", class(after_in_results$sampleID), "\n")
        cat("  seq_depth_df sampleIDs (first 5):", paste(head(seq_depth_df$sampleID, 5), collapse = ", "), "\n")
        cat("  seq_depth_df sampleID type:", class(seq_depth_df$sampleID), "\n")
        test_match <- "ABBY_001-M-20170607-COMP" %in% seq_depth_df$sampleID
        cat("  Test ID in seq_depth_df:", test_match, "\n")
        if(test_match) {
            cat("  Test ID seq_depth:", seq_depth_df$seq_depth[seq_depth_df$sampleID == "ABBY_001-M-20170607-COMP"], "\n")
        }
    }
    
    # Ensure both are character for join
    all_results$sampleID <- as.character(trimws(all_results$sampleID))
    seq_depth_df$sampleID <- as.character(trimws(seq_depth_df$sampleID))
    
    all_results <- all_results %>%
        select(-any_of(c("seq_depth", "seq_depth.x", "seq_depth.y", "pct_classified"))) %>%
        left_join(seq_depth_df %>% select(sampleID, seq_depth), by = "sampleID") %>%
        mutate(
            reads_classified = replace_na(reads_classified, 0),
            pct_classified = (reads_classified / seq_depth) * 100
        )
    
    # Debug in test mode
    if(!is.null(TEST_SAMPLE_ID)) {
        after_after_join <- all_results %>% filter(filter_status == "after")
        cat("\nAfter join:\n")
        cat("  Records with seq_depth:", sum(!is.na(after_after_join$seq_depth)), "out of", nrow(after_after_join), "\n")
        if(sum(!is.na(after_after_join$seq_depth)) > 0) {
            cat("  Sample with seq_depth:\n")
            print(after_after_join %>% filter(!is.na(seq_depth)) %>% head(2) %>% select(sampleID, seq_depth, pct_classified))
        }
    }
    
    missing_seq_depth <- sum(is.na(all_results$seq_depth))
    if(missing_seq_depth > 0) {
        cat("⚠ WARNING:", missing_seq_depth, "records missing seq_depth\n")
        if(!is.null(TEST_SAMPLE_ID)) {
            missing_after <- all_results %>% filter(filter_status == "after" & is.na(seq_depth))
            cat("  Missing sampleIDs (after):", paste(unique(missing_after$sampleID), collapse = ", "), "\n")
        }
    }
} else {
    stop("No results generated!")
}

# Aggregate results
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

classification_by_rank <- bind_rows(all_domain_summary, fungi_summary) %>%
    arrange(taxon_group, db_name, taxonomic_rank, filter_status)

# Save results
if(!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
}
write_csv(classification_by_rank, summary_output)
cat("✅ Saved classification percentages by rank to:", summary_output, "\n")

write_csv(all_results, output_file)
cat("✅ Saved per-sample classification percentages by rank to:", output_file, "\n")
cat("   Total records:", nrow(all_results), "\n")

cat("✅ Script completed successfully\n")
