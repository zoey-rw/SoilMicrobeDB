#!/usr/bin/env Rscript
# 04: Extract and reshape scoring information from Architeuthis _scores.output files
# Processes scoring files incrementally and appends results to a common file
# Allows deletion of _scores.output files after information is extracted
#
# Usage: Rscript scripts/run_workflow/04_reshape_score_reads.r
#
# Input:  *_scores.output files from Step 2
# Output: filter_results_summary.csv (preserved, used by visualization scripts)

library(tidyverse)
library(data.table)
library(future.apply)

# Helper function to normalize samp_name (remove trailing underscores and normalize format)
normalize_samp_name = function(samp_name) {
    # Remove trailing underscores
    samp_name = sub("_+$", "", samp_name)
    # Remove trailing underscores after database names
    samp_name = sub("_(gtdb_207_unfiltered|gtdb_207|soil_microbe_db|pluspf)_+$", "_\\1", samp_name)
    return(samp_name)
}

# Function to read and summarize large files in chunks
# Returns summary dataframe or NULL if error
read_and_summarize_chunked = function(file_path, samp_name, seq_depth_df, 
                                      chunk_size = 5e6, max_entropy = 0.1,
                                      max_multiplicity = 2, min_consistency = 0.9) {
    # Normalize samp_name (remove trailing underscores)
    samp_name = normalize_samp_name(samp_name)
    
    tryCatch({
        # Initialize aggregation variables
        total_rows = 0
        sum_consistency = 0
        sum_multiplicity = 0
        sum_entropy = 0
        sum_confidence = 0
        n_valid_rows = 0
        n_passing = 0
        
        # Read header to get column names
        header_chunk = NULL
        tryCatch({
            header_chunk = suppressWarnings(data.table::fread(
                file_path, 
                nrows = 1,
                showProgress = FALSE
            ))
        }, error = function(e) {
            if(grepl("character strings are limited to 2\\^31-1 bytes", e$message)) {
                cat("    ERROR: File contains lines too long for R to process\n")
            } else {
                cat("    ERROR reading header:", e$message, "\n")
            }
        })
        
        if(is.null(header_chunk) || nrow(header_chunk) == 0) {
            return(NULL)
        }
        
        # Check required columns exist
        required_cols = c("n_kmers", "consistency", "multiplicity", "entropy", "confidence")
        missing_cols = setdiff(required_cols, names(header_chunk))
        if(length(missing_cols) > 0) {
            cat("    ERROR: Missing required columns:", paste(missing_cols, collapse=", "), "\n")
            return(NULL)
        }
        
        # Read first data chunk (skip header row)
        first_chunk = NULL
        tryCatch({
            first_chunk = suppressWarnings(data.table::fread(
                file_path, 
                skip = 1,
                nrows = chunk_size,
                showProgress = FALSE
            ))
        }, error = function(e) {
            if(grepl("character strings are limited to 2\\^31-1 bytes", e$message)) {
                cat("    ERROR: File contains lines too long for R to process\n")
            } else {
                cat("    ERROR reading first chunk:", e$message, "\n")
            }
        })
        
        if(is.null(first_chunk) || nrow(first_chunk) == 0) {
            return(NULL)
        }
        
        # Set column names from header
        if(!identical(names(first_chunk), names(header_chunk))) {
            names(first_chunk) = names(header_chunk)
        }
        
        # Process first chunk
        valid_chunk = first_chunk %>% filter(n_kmers > 0)
        if(nrow(valid_chunk) > 0) {
            sum_consistency = sum_consistency + sum(valid_chunk$consistency, na.rm = TRUE)
            sum_multiplicity = sum_multiplicity + sum(valid_chunk$multiplicity, na.rm = TRUE)
            sum_entropy = sum_entropy + sum(valid_chunk$entropy, na.rm = TRUE)
            sum_confidence = sum_confidence + sum(valid_chunk$confidence, na.rm = TRUE)
            n_valid_rows = n_valid_rows + nrow(valid_chunk)
            
            # Count passing reads
            passing = valid_chunk %>%
                filter(multiplicity <= max_multiplicity & 
                       consistency >= min_consistency & 
                       entropy <= max_entropy)
            n_passing = n_passing + nrow(passing)
        }
        
        total_rows = total_rows + nrow(first_chunk)
        
        # Check if we've already read the entire file (first chunk was smaller than chunk_size)
        # If so, we don't need to read more chunks
        if(nrow(first_chunk) < chunk_size) {
            # Already read entire file in first chunk
        } else {
            # Read remaining chunks (skip header + rows already read)
            skip_rows = 1 + nrow(first_chunk)  # 1 for header + rows already read
            chunk_num = 1
            
            while(TRUE) {
                chunk_num = chunk_num + 1
                chunk_df = NULL
                tryCatch({
                    chunk_df = suppressWarnings(data.table::fread(
                        file_path, 
                        skip = skip_rows, 
                        nrows = chunk_size,
                        showProgress = FALSE
                    ))
                }, error = function(e) {
                    if(grepl("character strings are limited to 2\\^31-1 bytes", e$message)) {
                        cat("    ERROR: File contains lines too long for R to process at chunk", chunk_num, "\n")
                    } else if(grepl("skip=.*but the input only has", e$message)) {
                        # End of file reached - this is expected, not an error
                        # Don't print error, just break
                    } else {
                        cat("    ERROR reading chunk", chunk_num, ":", e$message, "\n")
                    }
                    # Set chunk_df to NULL to signal end of file
                    chunk_df <<- NULL
                })
                
                if(is.null(chunk_df) || nrow(chunk_df) == 0) {
                    break  # End of file
                }
                
                # Set column names from header (fread with skip might not preserve them)
                if(!identical(names(chunk_df), names(header_chunk))) {
                    if(length(names(chunk_df)) == length(names(header_chunk))) {
                        names(chunk_df) = names(header_chunk)
                    } else {
                        cat("    WARNING: Chunk", chunk_num, "column count mismatch, skipping remainder\n")
                        break
                    }
                }
                
                # Check if required columns exist
                if(!all(required_cols %in% names(chunk_df))) {
                    cat("    WARNING: Chunk", chunk_num, "missing required columns, skipping remainder\n")
                    break
                }
                
                # Aggregate statistics from this chunk
                valid_chunk = chunk_df %>% filter(n_kmers > 0)
                if(nrow(valid_chunk) > 0) {
                    sum_consistency = sum_consistency + sum(valid_chunk$consistency, na.rm = TRUE)
                    sum_multiplicity = sum_multiplicity + sum(valid_chunk$multiplicity, na.rm = TRUE)
                    sum_entropy = sum_entropy + sum(valid_chunk$entropy, na.rm = TRUE)
                    sum_confidence = sum_confidence + sum(valid_chunk$confidence, na.rm = TRUE)
                    n_valid_rows = n_valid_rows + nrow(valid_chunk)
                    
                    # Count passing reads
                    passing = valid_chunk %>%
                        filter(multiplicity <= max_multiplicity & 
                               consistency >= min_consistency & 
                               entropy <= max_entropy)
                    n_passing = n_passing + nrow(passing)
                }
                
                total_rows = total_rows + nrow(chunk_df)
                skip_rows = skip_rows + nrow(chunk_df)
                
                if(nrow(chunk_df) < chunk_size) {
                    break  # Last chunk
                }
                
                if(chunk_num %% 20 == 0) {
                    cat("    Processed", chunk_num, "chunks (", format(total_rows, big.mark=","), "rows)...\n")
                }
            }  # end while
        }  # end else
        
        if(n_valid_rows == 0) {
            return(NULL)
        }
        
        # Create summary from aggregated statistics
        db_name_val = case_when(
            grepl("gtdb_207_unfiltered", samp_name) ~ "gtdb_207_unfiltered",
            grepl("gtdb_207", samp_name) & !grepl("gtdb_207_unfiltered", samp_name) ~ "gtdb_207",
            grepl("soil_microbe_db", samp_name) ~ "soil_microbe_db",
            grepl("pluspf", samp_name) ~ "pluspf",
            TRUE ~ NA_character_
        )
        
        sampleID = sub("_(gtdb_207_unfiltered|gtdb_207|soil_microbe_db|pluspf)$", "", samp_name)
        
        # Calculate means
        mean_consistency = sum_consistency / n_valid_rows
        mean_multiplicity = sum_multiplicity / n_valid_rows
        mean_entropy = sum_entropy / n_valid_rows
        mean_confidence = sum_confidence / n_valid_rows
        pct_passing = n_passing / n_valid_rows
        
        # Get sequencing depth
        seq_depth_row = seq_depth_df %>% filter(sampleID == !!sampleID)
        if(nrow(seq_depth_row) > 0) {
            seq_depth = seq_depth_row$seq_depth[1]
            percent_classified = n_valid_rows / seq_depth
            percent_passing_overall = n_passing / seq_depth
        } else {
            seq_depth = NA
            percent_classified = NA
            percent_passing_overall = NA
        }
        
        # Create summary dataframe in the expected format
        file_summary = tibble(
            db_name = db_name_val,
            sampleID = sampleID,
            samp_name = samp_name,
            metric = c("mean_consistency", "mean_multiplicity", "mean_entropy", "mean_confidence", 
                      "percent_classified", "percent_passing"),
            value = c(mean_consistency, mean_multiplicity, mean_entropy, mean_confidence,
                     percent_classified, percent_passing_overall)
        ) %>% filter(!is.na(value))
        
        return(file_summary)
    }, error = function(e) {
        cat("    ERROR in chunked reading:", e$message, "\n")
        return(NULL)
    })
}

# Function to summarize the "scores" output from Architeuthis
# Format is one row per classified read
summarize_filter_scores = function(scores_df, seq_depth_df, max_entropy = .1,
                                   max_multiplicity = 2,
                                   min_consistency = .9) {
    
    scores_df = scores_df %>% mutate(db_name = ifelse(grepl("soil_microbe_db", samp_name), "soil_microbe_db",
                                                      ifelse(grepl("pluspf", samp_name), "pluspf",
                                                             ifelse(grepl("gtdb_207_unfiltered", samp_name), "gtdb_207_unfiltered",
                                                                    ifelse(grepl("gtdb", samp_name), "gtdb_207", NA
                                                                    )))))
    
    # Add column for passing or not
    scores_df = scores_df %>% mutate(pass_filter =
                                         ifelse(multiplicity <= max_multiplicity &
                                                    consistency >= min_consistency &
                                                    entropy <= max_entropy, 1, 0))
    
    # Summarize per-read scores
    score_summary_df = scores_df  %>%
        filter(n_kmers > 0) %>%
        group_by(samp_name) %>%
        reframe(mean_consistency = mean(consistency, na.rm=T),
                mean_multiplicity = mean(multiplicity, na.rm=T),
                mean_entropy = mean(entropy, na.rm=T),
                mean_confidence = mean(confidence,  na.rm=T)) %>% distinct
    
    # Percent passing or not
    filter_summary_df = scores_df %>%
        group_by(db_name, samp_name) %>%
        add_tally(name = "n_scored_reads") %>%
        group_by(db_name, samp_name, pass_filter) %>%
        add_tally(name = "n_reads") %>%
        group_by(db_name, samp_name, pass_filter, n_scored_reads, n_reads) %>%
        reframe(pct_of_classified_passing = n_reads/n_scored_reads) %>% 
        filter(pass_filter == 1) %>% ungroup %>%
        select(-pass_filter) %>% distinct() %>%
        mutate(
            # Remove database suffix from samp_name to get sampleID
            # Database names are: gtdb_207, gtdb_207_unfiltered, soil_microbe_db, pluspf
            # _filtered in filenames (like _filtered_kraken.kreport) is about file type, not database name
            sampleID = sub("_(gtdb_207_unfiltered|gtdb_207|soil_microbe_db|pluspf)$", "", samp_name)
        ) %>%
        mutate(n_total_classified_reads = n_scored_reads) %>%
        select(-n_scored_reads)
    
    # Include sequencing depth to get overall percent
    filter_summary_df = left_join(filter_summary_df, seq_depth_df, by=join_by(sampleID)) %>%
        mutate(percent_classified = n_total_classified_reads / seq_depth,
               percent_passing = n_reads / seq_depth)
    
    # Combine and reshape for output
    # Check for duplicates before joining to avoid cartesian products
    score_summary_df_unique = score_summary_df %>% distinct(samp_name, .keep_all = TRUE)
    filter_summary_df_unique = filter_summary_df %>% distinct(samp_name, .keep_all = TRUE)
    
    # Use inner_join to avoid cartesian products that can occur with full_join
    # This ensures we only keep samples present in both dataframes
    out_df = inner_join(score_summary_df_unique, filter_summary_df_unique, by=join_by(samp_name)) %>% 
        pivot_longer(cols=-c(db_name, sampleID, samp_name), names_to = "metric", values_drop_na = TRUE) %>% 
        distinct()
    
    out_df
}

# Set up paths - check both new and old locations, plus HARDDRIVE
filter_reads_dir_local = "data/NEON_metagenome_classification/02_bracken_output"
filter_reads_dir_local_old = "data/classification/02_bracken_output"
filter_reads_dir_harddrive = "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/02_bracken_output"
filter_reads_dir_cluster = "/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/02_bracken_output"

# Output files - check if old location exists, otherwise use new location
if(dir.exists("data/classification/analysis_files")) {
    output_file = "data/classification/analysis_files/filter_results_summary.csv"
    processed_log_file = "data/classification/analysis_files/filter_results_processed_files.txt"
} else {
    output_file = "data/summary_files/filter_results_summary.csv"
    processed_log_file = "data/summary_files/filter_results_processed_files.txt"
}

# Create output directory if it doesn't exist
output_dir = dirname(output_file)
if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
}

# Check which directories exist
dirs_to_search = character(0)
if(dir.exists(filter_reads_dir_local)) {
    dirs_to_search = c(dirs_to_search, filter_reads_dir_local)
    cat("✓ Found new location:", filter_reads_dir_local, "\n")
}
if(dir.exists(filter_reads_dir_local_old)) {
    dirs_to_search = c(dirs_to_search, filter_reads_dir_local_old)
    cat("✓ Found old location:", filter_reads_dir_local_old, "\n")
}
if(dir.exists(filter_reads_dir_harddrive)) {
    dirs_to_search = c(dirs_to_search, filter_reads_dir_harddrive)
    cat("✓ Found HARDDRIVE directory:", filter_reads_dir_harddrive, "\n")
}
if(dir.exists(filter_reads_dir_cluster)) {
    dirs_to_search = c(dirs_to_search, filter_reads_dir_cluster)
    cat("✓ Found cluster directory:", filter_reads_dir_cluster, "\n")
}

if(length(dirs_to_search) == 0) {
    stop("❌ MISSING DIRECTORY: Filter scores directory not found!\n",
         "   Checked: ", filter_reads_dir_local, "\n",
         "   Checked: ", filter_reads_dir_local_old, "\n",
         "   Checked: ", filter_reads_dir_harddrive, "\n",
         "   This directory should contain Architeuthis filter score files (*_scores.output).\n",
         "   Please add this directory and the score files to continue.")
}

# Load sequencing depth data - check both new and old locations
seq_depth_file = "data/NEON_metagenome_classification/seq_depth_df.rds"
if(!file.exists(seq_depth_file)) {
    seq_depth_file = "data/classification/analysis_files/seq_depth_df.rds"
}
seq_depth_df <- readRDS(seq_depth_file)
# Remove columns if they exist
if("db_name" %in% names(seq_depth_df)) {
    seq_depth_df <- seq_depth_df %>% select(-db_name)
}
if("identified_reads" %in% names(seq_depth_df)) {
    seq_depth_df <- seq_depth_df %>% select(-identified_reads)
}

# Find all scoring files in all directories
cat("Searching for scoring files...\n")
filter_scores_list = character(0)
for(dir_path in dirs_to_search) {
    files_in_dir = list.files(dir_path, 
                              pattern = "_scores.output", 
                              recursive = T, 
                              full.names = T)
    filter_scores_list = c(filter_scores_list, files_in_dir)
    cat("  Found", length(files_in_dir), "files in", dir_path, "\n")
}

if(length(filter_scores_list) == 0) {
    stop("❌ MISSING FILES: No filter score files found!\n",
         "   Expected files matching pattern '*_scores.output' in: ", filter_reads_dir, "\n",
         "   These are output files from Architeuthis filter analysis.\n",
         "   Please add the score files to continue.")
}


cat("Found", length(filter_scores_list), "scoring files\n")

# Diagnostic: Show file distribution by directory
if(length(filter_scores_list) > 0) {
    cat("\nFile distribution by directory:\n")
    for(dir_path in dirs_to_search) {
        files_in_this_dir = sum(grepl(paste0("^", dir_path), filter_scores_list))
        if(files_in_this_dir > 0) {
            cat("  ", dir_path, ":", files_in_this_dir, "files\n")
        }
    }
    
    # Check if we're finding enough files (expected: 500+ per database = 2000+ total)
    expected_min_files = 1500  # Conservative estimate: 4 databases * ~400 samples each
    if(length(filter_scores_list) < expected_min_files) {
        cat("\n⚠️  WARNING: Found only", length(filter_scores_list), "files, but expected", expected_min_files, "+ files\n")
        cat("  This suggests files may be in a location not being searched.\n")
        cat("  If files are on the cluster, this script may need to be run on the cluster.\n")
        cat("  Or files may be in a different directory structure.\n")
    }
}

# Load list of already processed files
processed_files = character(0)
processed_basenames = character(0)
if(file.exists(processed_log_file)) {
    processed_files = readLines(processed_log_file)
    # Normalize processed file paths (convert relative to absolute if they exist)
    processed_files_normalized = sapply(processed_files, function(x) {
        if(file.exists(x)) {
            normalizePath(x)
        } else {
            # If file doesn't exist, try to normalize just the basename match
            x
        }
    })
    processed_basenames = basename(processed_files_normalized)
    cat("Found", length(processed_files), "already processed files\n")
}

# Normalize all found files to absolute paths for comparison
filter_scores_list_normalized = normalizePath(filter_scores_list)
filter_scores_basenames = basename(filter_scores_list_normalized)

# Load existing results if file exists
existing_results = NULL
existing_samp_names = character(0)

if(file.exists(output_file)) {
    cat("Loading existing results from:", output_file, "\n")
    existing_results = data.table::fread(output_file, showProgress = FALSE)
} else {
    cat("No existing results file found. Starting fresh.\n")
}

if(!is.null(existing_results)) {
    # Normalize samp_names in existing results (remove trailing underscores)
    if("samp_name" %in% names(existing_results)) {
        existing_results = existing_results %>%
            mutate(samp_name = normalize_samp_name(samp_name))
    }
    
    # Fix incorrect sampleID values by extracting from samp_name
    if("samp_name" %in% names(existing_results) && "sampleID" %in% names(existing_results)) {
        cat("  Fixing sampleID values extracted from samp_name...\n")
        existing_results = existing_results %>%
            mutate(
                # Extract sampleID from samp_name (remove database suffix)
                # Database names are: gtdb_207, gtdb_207_unfiltered, soil_microbe_db, pluspf
                sampleID = sub("_(gtdb_207_unfiltered|gtdb_207|soil_microbe_db|pluspf)$", "", samp_name)
            )
        
        # Count how many were fixed
        invalid_sampleIDs = existing_results %>%
            filter(sampleID %in% c("gtdb_207", "gtdb_207_unfiltered", "soil_microbe_db", "pluspf") | 
                   grepl("^(ABBY|BARR|BART|BLAN|CLBJ|CPER|DSNY|GRSM|GUAN|HARV|ORNL|TALL)", sampleID) == FALSE)
        
        if(nrow(invalid_sampleIDs) > 0) {
            cat("  Fixed", nrow(invalid_sampleIDs), "rows with incorrect sampleID values\n")
        }
    }
    
    if("samp_name" %in% names(existing_results)) {
        existing_samp_names = unique(existing_results$samp_name)
        cat("  Found", length(existing_samp_names), "unique samp_names in existing results\n")
    }
}

# ============================================================================
# STEP 1: Build comprehensive inventory of expected vs actual state
# ============================================================================
cat("\n=== Building comprehensive file inventory ===\n")

# Load lineage files to determine what samples SHOULD exist
lineage_dir_new = "data/NEON_metagenome_classification/summary_files"
lineage_dir_old = "data/classification/taxonomic_rank_summaries/species"
databases_to_check = c("soil_microbe_db", "pluspf", "gtdb_207", "gtdb_207_unfiltered")

expected_samples = tibble()
for(db in databases_to_check) {
    filename = paste0(db, "_species_merged_lineage.csv")
    lineage_file = NULL
    
    if(file.exists(file.path(lineage_dir_new, filename))) {
        lineage_file = file.path(lineage_dir_new, filename)
    } else if(file.exists(file.path(lineage_dir_old, filename))) {
        lineage_file = file.path(lineage_dir_old, filename)
    }
    
    if(!is.null(lineage_file) && file.exists(lineage_file)) {
        tryCatch({
            lineage_df = read_csv(lineage_file, col_select = "sample_id", show_col_types = FALSE)
            if("sample_id" %in% names(lineage_df) && nrow(lineage_df) > 0) {
                lineage_samples = unique(lineage_df$sample_id)
                
                # Extract sampleID and determine expected samp_name format
                # Lineage files have format: {sampleID}-COMP_{db} or {sampleID}_{db}
                lineage_sampleIDs = sapply(lineage_samples, function(samp_name) {
                    samp_name = sub("_filtered$", "", samp_name)  # Remove _filtered if present
                    # Remove database suffix to get base sampleID
                    base_id = sub("_(soil_microbe_db|pluspf|gtdb_207_unfiltered|gtdb_207)$", "", samp_name)
                    # Remove -COMP if present (some lineage files already have it)
                    base_id = sub("-COMP$", "", base_id)
                    return(base_id)
                })
                
                # Build expected samp_names (format: {sampleID}-COMP_{db})
                # Note: Actual files use this format, so we match it
                expected_samp_names = sapply(lineage_sampleIDs, function(samp_id) {
                    if(db == "gtdb_207_unfiltered") {
                        paste0(samp_id, "-COMP_gtdb_207_unfiltered")
                    } else if(db == "gtdb_207") {
                        paste0(samp_id, "-COMP_gtdb_207")
                    } else {
                        paste0(samp_id, "-COMP_", db)
                    }
                })
                # Normalize to remove any trailing underscores
                expected_samp_names = normalize_samp_name(expected_samp_names)
                
                expected_samples = bind_rows(
                    expected_samples,
                    tibble(
                        sampleID = lineage_sampleIDs,
                        samp_name = expected_samp_names,
                        db_name = db,
                        expected_scores_file = paste0(expected_samp_names, "_scores.output")
                    )
                )
            }
        }, error = function(e) {
            cat("  Warning: Could not read lineage file for", db, ":", e$message, "\n")
        })
    }
}

cat("  Expected samples from lineage files:", nrow(expected_samples), "\n")
if(nrow(expected_samples) > 0) {
    cat("  Distribution by database:\n")
    db_counts = expected_samples %>% count(db_name)
    print(db_counts)
}

# Build inventory of found files
found_files_inventory = tibble(
    file_path = filter_scores_list,
    basename = filter_scores_basenames,
    samp_name = normalize_samp_name(sub("_scores.output$", "", filter_scores_basenames))
) %>%
    mutate(
        db_name = case_when(
            grepl("soil_microbe_db", samp_name) ~ "soil_microbe_db",
            grepl("pluspf", samp_name) ~ "pluspf",
            grepl("gtdb_207_unfiltered", samp_name) ~ "gtdb_207_unfiltered",
            grepl("gtdb_207", samp_name) ~ "gtdb_207",  # More specific: gtdb_207 (not unfiltered)
            TRUE ~ NA_character_
        )
    )

cat("  Found files:", nrow(found_files_inventory), "\n")
if(nrow(found_files_inventory) > 0) {
    cat("  Distribution by database:\n")
    found_db_counts = found_files_inventory %>% 
        filter(!is.na(db_name)) %>%
        count(db_name)
    print(found_db_counts)
}

# Build inventory of processed log entries
processed_inventory = tibble(
    basename = processed_basenames,
    samp_name = normalize_samp_name(sub("_scores.output$", "", processed_basenames))
) %>%
    mutate(
        db_name = case_when(
            grepl("soil_microbe_db", samp_name) ~ "soil_microbe_db",
            grepl("pluspf", samp_name) ~ "pluspf",
            grepl("gtdb_207_unfiltered", samp_name) ~ "gtdb_207_unfiltered",
            grepl("gtdb_207", samp_name) ~ "gtdb_207",  # More specific: gtdb_207 (not unfiltered)
            TRUE ~ NA_character_
        )
    )

cat("  Files in processed log:", nrow(processed_inventory), "\n")

# Build inventory of existing results
results_inventory = NULL
if(!is.null(existing_results) && "samp_name" %in% names(existing_results) && "db_name" %in% names(existing_results)) {
    results_inventory = existing_results %>%
        filter(db_name %in% databases_to_check) %>%
        distinct(samp_name, db_name) %>%
        mutate(in_results = TRUE)
    
    cat("  Unique samp_names in results:", nrow(results_inventory), "\n")
    if(nrow(results_inventory) > 0) {
        cat("  Distribution by database:\n")
        results_db_counts = results_inventory %>% count(db_name)
        print(results_db_counts)
    }
}

# ============================================================================
# STEP 2: Make intelligent processing decisions
# ============================================================================
cat("\n=== Determining which files need processing ===\n")

# Start with all found files
files_to_process = found_files_inventory %>%
    mutate(
        in_processed_log = basename %in% processed_inventory$basename,
        in_results = FALSE,
        results_db_name = NA_character_,
        needs_processing = FALSE,
        reason = NA_character_
    )

# Check which files are in results
# Note: results_db_name is already initialized above, so we just need to populate it
if(!is.null(results_inventory) && nrow(results_inventory) > 0) {
    # Join to get db_name from results (will create db_name column from results)
    results_lookup = results_inventory %>% 
        select(samp_name, db_name)
    
    files_to_process = files_to_process %>%
        left_join(results_lookup, by = "samp_name", suffix = c("", "_from_results")) %>%
        mutate(
            results_db_name = if_else(!is.na(db_name_from_results), db_name_from_results, results_db_name),
            in_results = !is.na(db_name_from_results)
        ) %>%
        select(-db_name_from_results)  # Remove temporary column
}
# If results_inventory is NULL or empty, results_db_name stays as NA_character_ (already initialized)

# Decision logic: Process all files that are NOT in results
# Note: Files with wrong db_name in results will be fixed in existing results, not reprocessed
# We process ALL files that exist and aren't in results (regardless of processed log status)
files_to_process = files_to_process %>%
    mutate(
        # Process files that are NOT in results
        # Files in results with wrong db_name will be fixed separately, not reprocessed
        needs_processing = !in_results,
        
        reason = case_when(
            !in_results ~ "Not in results",
            in_results & !is.na(db_name) & !is.na(results_db_name) & db_name != results_db_name ~ "Wrong db_name in results (will fix in results, not reprocess)",
            TRUE ~ NA_character_
        )
    )

# Count decisions
new_files_count = sum(!files_to_process$in_results & files_to_process$needs_processing, na.rm = TRUE)
wrong_db_count = sum(files_to_process$in_results & !is.na(files_to_process$db_name) & 
                     !is.na(files_to_process$results_db_name) & 
                     files_to_process$db_name != files_to_process$results_db_name, na.rm = TRUE)
total_needs_processing = sum(files_to_process$needs_processing, na.rm = TRUE)

cat("  New files to process (not in results):", new_files_count, "\n")
cat("  Files with wrong db_name in results:", wrong_db_count, "\n")
cat("    (db_name will be fixed in existing results, files won't be reprocessed)\n")
cat("  Total files to process:", new_files_count, "\n")

# Show examples
if(total_needs_processing > 0) {
    examples = files_to_process %>%
        filter(needs_processing) %>%
        head(5) %>%
        select(samp_name, reason)
    cat("  Examples:\n")
    for(i in 1:nrow(examples)) {
        cat("    -", examples$samp_name[i], ":", examples$reason[i], "\n")
    }
}

# Extract files to process - only files that are NOT in results
# Files with wrong db_name will be fixed in existing results, not reprocessed
new_files = files_to_process %>%
    filter(needs_processing & !in_results) %>%
    pull(file_path)

cat("\n=== Final processing queue ===\n")
cat("  Files to process:", length(new_files), "\n")

# Diagnostic: Show breakdown of why files are/aren't being processed
cat("\n=== Diagnostic: File processing breakdown ===\n")
cat("  Total files found:", nrow(files_to_process), "\n")
cat("  Files in processed log:", sum(files_to_process$in_processed_log, na.rm = TRUE), "\n")
cat("  Files in results:", sum(files_to_process$in_results, na.rm = TRUE), "\n")
cat("  Files NOT in results:", sum(!files_to_process$in_results, na.rm = TRUE), "\n")
cat("  Files marked as needs_processing:", sum(files_to_process$needs_processing, na.rm = TRUE), "\n")
cat("  Files to actually process (not in results):", length(new_files), "\n")

# Show sample of files in results vs not in results
if(sum(files_to_process$in_results, na.rm = TRUE) > 0) {
    cat("\n  Sample of files IN results:\n")
    in_results_sample = files_to_process %>% 
        filter(in_results) %>% 
        head(3) %>% 
        select(samp_name, db_name, results_db_name)
    print(in_results_sample)
}

if(sum(!files_to_process$in_results, na.rm = TRUE) > 0) {
    cat("\n  Sample of files NOT in results:\n")
    not_in_results_sample = files_to_process %>% 
        filter(!in_results) %>% 
        head(10) %>% 
        select(samp_name, db_name, in_processed_log)
    print(not_in_results_sample)
}

# Check if expected samples from lineage are being found
if(nrow(expected_samples) > 0) {
    cat("\n  Expected samples from lineage:", nrow(expected_samples), "\n")
    expected_samp_names = expected_samples$samp_name
    found_samp_names = files_to_process$samp_name
    
    missing_expected = setdiff(expected_samp_names, found_samp_names)
    if(length(missing_expected) > 0) {
        cat("  ⚠️  Found", length(missing_expected), "expected samples from lineage that don't have files\n")
        cat("  Examples:", paste(head(missing_expected, 5), collapse=", "), "\n")
        
        # Check if these expected samples have files that exist but weren't found
        cat("  Checking if files exist for missing expected samples...\n")
        missing_with_files = character(0)
        for(samp_name in head(missing_expected, 50)) {  # Check first 50
            expected_file = paste0(samp_name, "_scores.output")
            for(dir_path in dirs_to_search) {
                if(dir.exists(dir_path)) {
                    file_path = file.path(dir_path, expected_file)
                    if(file.exists(file_path)) {
                        missing_with_files = c(missing_with_files, file_path)
                        break
                    }
                }
            }
        }
        
        if(length(missing_with_files) > 0) {
            cat("  ⚠️  Found", length(missing_with_files), "files for expected samples that weren't in initial file search!\n")
            cat("  Adding these files to processing queue...\n")
            new_files = c(new_files, missing_with_files)
            new_files = unique(new_files)
            cat("  ➕ Added", length(missing_with_files), "missing files. Total files to process now:", length(new_files), "\n")
        }
    }
    
    found_but_not_expected = setdiff(found_samp_names, expected_samp_names)
    if(length(found_but_not_expected) > 0) {
        cat("  Found", length(found_but_not_expected), "files that aren't in expected samples from lineage\n")
    }
    
    # Check if expected samples that ARE found are incorrectly marked as "in results"
    expected_found = intersect(expected_samp_names, found_samp_names)
    if(length(expected_found) > 0) {
        expected_found_in_results = files_to_process %>%
            filter(samp_name %in% expected_found & in_results) %>%
            nrow()
        expected_found_not_in_results = files_to_process %>%
            filter(samp_name %in% expected_found & !in_results) %>%
            nrow()
        
        cat("  Expected samples that were found:", length(expected_found), "\n")
        cat("    In results:", expected_found_in_results, "\n")
        cat("    NOT in results:", expected_found_not_in_results, "\n")
        
        if(expected_found_not_in_results < length(missing_expected)) {
            cat("  ⚠️  WARNING: Many expected samples are marked as 'in results' but diagnostic shows they're missing!\n")
            cat("  This suggests samp_name matching between files and results might be incorrect.\n")
        }
    }
}

# ============================================================================
# STEP 3: Fix db_name in existing results (for records with wrong db_name)
# ============================================================================
if(!is.null(existing_results) && "db_name" %in% names(existing_results) && "samp_name" %in% names(existing_results)) {
    cat("\n=== Fixing db_name in existing results ===\n")
    
    # Find records with wrong db_name
    wrong_db_records = files_to_process %>%
        filter(in_results & !is.na(db_name) & !is.na(results_db_name) & db_name != results_db_name) %>%
        select(samp_name, db_name) %>%
        distinct()
    
    if(nrow(wrong_db_records) > 0) {
        cat("  Found", nrow(wrong_db_records), "samp_names with wrong db_name in results\n")
        cat("  Fixing db_name assignments...\n")
        
        # Fix db_name in existing results
        for(i in 1:nrow(wrong_db_records)) {
            samp = wrong_db_records$samp_name[i]
            correct_db = wrong_db_records$db_name[i]
            
            # Update db_name for this samp_name
            existing_results = existing_results %>%
                mutate(db_name = if_else(samp_name == samp, correct_db, db_name))
        }
        
        cat("  ✓ Fixed db_name for", nrow(wrong_db_records), "samp_names in existing results\n")
        
        # Save corrected results immediately
        cat("  Saving corrected results...\n")
        write_csv(existing_results, output_file)
        cat("  ✓ Corrected results saved to:", output_file, "\n")
    } else {
        cat("  ✓ No db_name corrections needed\n")
    }
}

# ============================================================================
# STEP 4: Process new files
# ============================================================================
if(length(new_files) == 0) {
    cat("\n✓ All files have already been processed. Exiting.\n")
} else {
    # Process new files incrementally
    cat("Processing", length(new_files), "files (saving incrementally)...\n")
    flush.console()  # Ensure message is visible immediately
    
    # Function to validate file before reading
    validate_file = function(file_path) {
        if(!file.exists(file_path)) {
            return(list(valid = FALSE, reason = "File does not exist"))
        }
        
        file_info = file.info(file_path)
        
        # Check if file is empty
        if(file_info$size == 0) {
            return(list(valid = FALSE, reason = "File is empty"))
        }
        
        # Skip very large files (> 2 GB) that exceed R's character string limits
        # These files cause "R character strings are limited to 2^31-1 bytes" errors
        # even with chunked reading
        file_size_gb = file_info$size / 1024^3
        if(file_size_gb > 2.0) {
            return(list(valid = FALSE, reason = sprintf("File too large (%.2f GB) - exceeds R processing limits", file_size_gb)))
        }
        
        # Check if file size is suspiciously round (multiple of 4096 bytes suggests truncation)
        if(file_info$size > 0 && file_info$size %% 4096 == 0 && file_info$size >= 4096) {
            # This could indicate truncation, but not always - check if file ends with newline
            # Read last few bytes to check
            con = file(file_path, "rb")
            seek(con, -min(100, file_info$size), origin = "end")
            last_bytes = readBin(con, "raw", min(100, file_info$size))
            close(con)
            
            # Check if file ends with newline
            if(length(last_bytes) > 0 && last_bytes[length(last_bytes)] != as.raw(10)) {
                return(list(valid = FALSE, reason = "File may be truncated (size is multiple of 4096 and doesn't end with newline)"))
            }
        }
        
        # Check if file is readable
        if(!file.access(file_path, mode = 4) == 0) {
            return(list(valid = FALSE, reason = "File is not readable"))
        }
        
        return(list(valid = TRUE, reason = ""))
    }
    
    start_time = Sys.time()
    
    # Load existing results once at the start (keep in memory for fast lookups)
    current_results = NULL
    existing_samp_names_set = character(0)
    if(file.exists(output_file)) {
        cat("Loading existing results into memory...\n")
        current_results = data.table::fread(output_file, showProgress = FALSE)
        if("samp_name" %in% names(current_results)) {
            existing_samp_names_set = unique(current_results$samp_name)
        }
        cat("  Loaded", nrow(current_results), "existing records\n")
    }
    
    # Process files incrementally, saving periodically instead of after each file
    save_interval = 10  # Save every 10 files (reduced I/O)
    print_interval = 5  # Print progress every 5 files
    processed_count = 0
    processed_files_list = character(0)
    new_summaries_list = list()  # Accumulate new summaries in memory
    
    cat("Processing", length(new_files), "files (saving every", save_interval, "files, progress every", print_interval, "files)...\n")
    
    # Process files one at a time and append incrementally
    for(i in 1:length(new_files)) {
        file_path = new_files[i]
        
        # Process file
        tryCatch({
            # Validate file before reading
            validation = validate_file(file_path)
            if(!validation$valid) {
                cat("WARNING: Skipping", basename(file_path), "-", validation$reason, "\n")
                next
            }
            
            # Read in the file with chunked reading for large files
            samp_name_from_file = normalize_samp_name(sub("_scores.output", "", basename(file_path)))
            file_size_gb = file.info(file_path)$size / 1024^3
            file_summary = NULL
            
            # Skip files > 2 GB that exceed R's processing limits
            if(file_size_gb > 2.0) {
                cat("WARNING: Skipping", basename(file_path), sprintf("- file too large (%.2f GB) for R to process\n", file_size_gb))
                next
            }
            
            # Use chunked reading for files > 500MB or if normal read fails
            if(file_size_gb > 0.5) {
                # Large file - use chunked reading
                cat("  Reading large file", basename(file_path), sprintf("(%.2f GB) in chunks...\n", file_size_gb))
                file_summary = read_and_summarize_chunked(file_path, samp_name_from_file, seq_depth_df)
                if(is.null(file_summary)) {
                    cat("WARNING: Failed to process", basename(file_path), "with chunked reading\n")
                    next
                }
            } else {
                # Small file - try normal reading first
                tryCatch({
                    df_in = suppressWarnings(data.table::fread(file_path, showProgress = FALSE))
                    if(is.null(df_in) || nrow(df_in) == 0) {
                        next
                    }
                    df_in$samp_name = samp_name_from_file
                    file_summary <- summarize_filter_scores(df_in, seq_depth_df)
                }, error = function(e) {
                    if(grepl("character strings are limited to 2\\^31-1 bytes", e$message)) {
                        # Fall back to chunked reading
                        cat("  File hit character limit, switching to chunked reading...\n")
                        file_summary <<- read_and_summarize_chunked(file_path, samp_name_from_file, seq_depth_df)
                        if(is.null(file_summary)) {
                            cat("WARNING: Failed to process", basename(file_path), "even with chunked reading\n")
                        }
                    } else {
                        cat("WARNING: Error reading", basename(file_path), ":", e$message, "\n")
                    }
                })
            }
            
            if(is.null(file_summary) || nrow(file_summary) == 0) {
                next
            }
            
            # Fix sampleID and db_name
            file_summary = file_summary %>%
                mutate(
                    sampleID = sub("_(gtdb_207_unfiltered|gtdb_207|soil_microbe_db|pluspf)$", "", samp_name),
                    db_name = case_when(
                        grepl("gtdb_207_unfiltered", samp_name) ~ "gtdb_207_unfiltered",
                        grepl("gtdb_207", samp_name) & !grepl("gtdb_207_unfiltered", samp_name) ~ "gtdb_207",
                        grepl("soil_microbe_db", samp_name) ~ "soil_microbe_db",
                        grepl("pluspf", samp_name) ~ "pluspf",
                        TRUE ~ db_name
                    )
                )
            
            # Check if this samp_name is already in results (fast in-memory lookup)
            samp_name_from_file = unique(file_summary$samp_name)
            already_in_results = FALSE
            if(length(samp_name_from_file) > 0 && samp_name_from_file[1] %in% existing_samp_names_set) {
                already_in_results = TRUE
            }
            
            if(already_in_results) {
                # Skip files that are already in results
                if(processed_count %% print_interval == 0 || i == length(new_files)) {
                    cat(sprintf("  [%d/%d] Skipped (already in results), %d new files processed\n", 
                               i, length(new_files), processed_count))
                    flush.console()
                }
                next
            }
            
            # Add to in-memory list (will write in batches)
            new_summaries_list[[length(new_summaries_list) + 1]] = file_summary
            existing_samp_names_set = c(existing_samp_names_set, samp_name_from_file)
            
            # Track processed file
            processed_files_list = c(processed_files_list, normalizePath(file_path))
            processed_count = processed_count + 1
            
            # Save periodically (every save_interval files) instead of after each file
            if(processed_count %% save_interval == 0 || i == length(new_files)) {
                # Combine all new summaries
                if(length(new_summaries_list) > 0) {
                    new_data = bind_rows(new_summaries_list)
                    
                    # Remove any existing records for these samp_names (in case of reprocessing)
                    if(!is.null(current_results)) {
                        current_results = current_results %>%
                            filter(!samp_name %in% new_data$samp_name)
                        combined_results = bind_rows(current_results, new_data)
                    } else {
                        combined_results = new_data
                    }
                    
                    # Write to file using faster fwrite
                    data.table::fwrite(combined_results, output_file)
                    
                    # Update in-memory current_results
                    current_results = combined_results
                    
                    # Clear the list
                    new_summaries_list = list()
                    
                    cat(sprintf("  [%d/%d] Saved %d NEW files (%d records total)\n", 
                               i, length(new_files), processed_count, nrow(combined_results)))
                }
            }
            
            # Update processed log periodically (every save_interval) instead of after each file
            if(processed_count %% save_interval == 0 || i == length(new_files)) {
                normalized_processed = if(length(processed_files) > 0) {
                    sapply(processed_files, function(x) {
                        if(file.exists(x)) normalizePath(x) else x
                    })
                } else {
                    character(0)
                }
                all_processed = unique(c(normalized_processed, processed_files_list))
                writeLines(all_processed, processed_log_file)
            }
            
            # Print progress every print_interval files
            if(processed_count %% print_interval == 0 || i == length(new_files)) {
                cat(sprintf("  [%d/%d] Processed %d files\n", 
                           i, length(new_files), processed_count))
                flush.console()
            }
            
        }, error = function(e) {
            cat("ERROR processing", basename(file_path), ":", e$message, "\n")
        })
    }
    
    end_time = Sys.time()
    cat("\n=== Processing completed ===\n")
    cat("Total processing time:", round(as.numeric(difftime(end_time, start_time, units = "mins")), 2), "minutes\n")
    cat("Successfully processed", processed_count, "files\n")
    
    # Load final results for compatibility with rest of script
    if(file.exists(output_file)) {
        new_summaries = data.table::fread(output_file, showProgress = FALSE)
        # Filter to only newly processed samples
        if(length(processed_files_list) > 0) {
            processed_samp_names = normalize_samp_name(sub("_scores.output$", "", basename(processed_files_list)))
            new_summaries = new_summaries %>%
                filter(samp_name %in% processed_samp_names)
        }
    } else {
        new_summaries = tibble()
    }
    
}  # End of if(length(new_files) > 0) block

# Results have already been merged incrementally into output_file
# Just load the final file for final corrections
cat("Loading final results for corrections...\n")
if(file.exists(output_file)) {
    all_summaries = data.table::fread(output_file, showProgress = FALSE)
    cat("  Loaded", nrow(all_summaries), "total records\n")
} else {
    # Fallback: use new_summaries if file doesn't exist (shouldn't happen)
    all_summaries = new_summaries
    if(!inherits(all_summaries, "tbl")) {
        all_summaries = as_tibble(all_summaries)
    }
    cat("  Using new_summaries (", nrow(all_summaries), " records)\n")
}

# Remove duplicates
if(nrow(all_summaries) > 0) {
    key_cols = intersect(c("samp_name", "metric", "db_name", "sampleID"), names(all_summaries))
    if(length(key_cols) > 0) {
        cat("  Removing duplicates using key columns:", paste(key_cols, collapse=", "), "\n")
        all_summaries = all_summaries %>%
            distinct(across(all_of(key_cols)), .keep_all = TRUE)
        cat("  Records after deduplication:", nrow(all_summaries), "\n")
    }
}

# Skip to final corrections - all the old diagnostic code below is not needed
# (Old code that was checking for missing files, etc. - removed for clarity)

# Final correction of database names in combined results
# Ensure gtdb_207_unfiltered and gtdb_207 are correctly assigned
if("samp_name" %in% names(all_summaries) && "db_name" %in% names(all_summaries)) {
    all_summaries = all_summaries %>%
        mutate(db_name = ifelse(grepl("gtdb_207_unfiltered", samp_name), "gtdb_207_unfiltered",
                               ifelse(db_name == "gtdb_207_filtered" & grepl("gtdb", samp_name), "gtdb_207", db_name)))
    cat("  Applied final database name corrections\n")
}

# Additional diagnostic: Check newest files and compare with results
cat("\n=== Diagnostic: Checking newest scores files ===\n")
if(length(filter_scores_list) > 0) {
    # Get file modification times
    file_info = file.info(filter_scores_list)
    file_info$file_path = rownames(file_info)
    file_info = file_info[order(file_info$mtime, decreasing = TRUE), ]
    
    # Show newest files
    newest_files = head(file_info, 20)
    cat("  Newest 20 _scores.output files (by modification time):\n")
    for(i in 1:min(20, nrow(newest_files))) {
        samp_name = normalize_samp_name(sub("_scores.output$", "", basename(newest_files$file_path[i])))
        in_results = if(!is.null(existing_results)) samp_name %in% existing_samp_names else FALSE
        status = if(in_results) "✓" else "✗ MISSING"
        cat(sprintf("    %d. %s %s (modified: %s)\n", 
                   i, status, basename(newest_files$file_path[i]), 
                   format(newest_files$mtime[i], "%Y-%m-%d %H:%M:%S")))
    }
    
    # Count how many newest files are missing
    newest_samp_names = normalize_samp_name(sub("_scores.output$", "", basename(newest_files$file_path)))
    if(!is.null(existing_results)) {
        newest_missing = sum(!newest_samp_names %in% existing_samp_names)
        cat("  Of the newest 20 files,", newest_missing, "are missing from results\n")
        
        if(newest_missing > 0) {
            missing_newest = newest_files[!newest_samp_names %in% existing_samp_names, ]
            cat("  ➕ Adding", nrow(missing_newest), "newest missing files to processing queue\n")
            new_files = c(new_files, missing_newest$file_path)
            new_files = unique(new_files)
        }
    }
}

# Additional diagnostic: Check if there are files that should be in results but aren't
cat("\n=== Diagnostic: Checking for missing samples ===\n")
cat("  Total files found:", length(filter_scores_list), "\n")
cat("  Files after processed log check:", length(new_files), "\n")

if(!is.null(existing_results) && length(existing_samp_names) > 0) {
    cat("  Existing results contain", length(existing_samp_names), "unique samp_names\n")
    
    # Get all samp_names from existing files
    all_file_samp_names = normalize_samp_name(sub("_scores.output$", "", filter_scores_basenames))
    cat("  Total unique samp_names in files:", length(unique(all_file_samp_names)), "\n")
    
    # Find files that exist but aren't in results
    missing_samp_names = setdiff(unique(all_file_samp_names), existing_samp_names)
    
    if(length(missing_samp_names) > 0) {
        cat("  ⚠️  Found", length(missing_samp_names), "samp_names in files but NOT in results\n")
        cat("  Example missing samp_names:", paste(head(missing_samp_names, 5), collapse=", "), "\n")
        if(length(missing_samp_names) > 5) {
            cat("  ... and", length(missing_samp_names) - 5, "more\n")
        }
        
        # Check if these are in the processed log
        missing_basenames = paste0(missing_samp_names, "_scores.output")
        in_processed_log = sum(missing_basenames %in% processed_basenames)
        cat("  Of these,", in_processed_log, "are incorrectly marked as processed in the log\n")
        
        # Find the actual files that need to be processed
        missing_files = filter_scores_list[filter_scores_basenames %in% missing_basenames]
        if(length(missing_files) > 0) {
            cat("  ➕ Adding", length(missing_files), "missing files to processing queue\n")
            new_files = c(new_files, missing_files)
            new_files = unique(new_files)  # Remove duplicates
            cat("  Total files to process now:", length(new_files), "\n")
        } else {
            cat("  ⚠️  WARNING: Found", length(missing_samp_names), "missing samp_names but couldn't find their files\n")
            cat("  This suggests files may be in a location not being searched\n")
        }
    } else {
        cat("  ✓ All file samp_names are in results\n")
        
        # But check if they're in results with correct db_name
        if("db_name" %in% names(existing_results)) {
            # Check if any files have samp_names that exist but with wrong db_name
            file_db_check = tibble(
                samp_name = all_file_samp_names,
                file_path = filter_scores_list
            ) %>%
                mutate(
                    file_db_name = case_when(
                        grepl("soil_microbe_db", samp_name) ~ "soil_microbe_db",
                        grepl("pluspf", samp_name) ~ "pluspf",
                        grepl("gtdb_207_unfiltered", samp_name) ~ "gtdb_207_unfiltered",
                        grepl("gtdb", samp_name) ~ "gtdb_207",
                        TRUE ~ NA_character_
                    )
                ) %>%
                filter(!is.na(file_db_name))
            
            # Check if samp_name exists but with wrong db_name
            existing_by_db = existing_results %>%
                filter(db_name %in% c("soil_microbe_db", "pluspf", "gtdb_207", "gtdb_207_unfiltered")) %>%
                distinct(samp_name, db_name) %>%
                rename(db_name_existing = db_name)
            
            wrong_db_files = file_db_check %>%
                left_join(existing_by_db, by = "samp_name") %>%
                filter(
                    !is.na(db_name_existing) &  # samp_name exists in results
                    file_db_name != db_name_existing  # but with wrong db_name
                )
            
            if(nrow(wrong_db_files) > 0) {
                cat("  ⚠️  Found", nrow(wrong_db_files), "files whose samp_names exist but with wrong db_name\n")
                cat("  These need to be reprocessed with correct db_name\n")
                new_files = c(new_files, wrong_db_files$file_path)
                new_files = unique(new_files)
                cat("  ➕ Added", nrow(wrong_db_files), "files with wrong db_name to processing queue\n")
                cat("  Total files to process now:", length(new_files), "\n")
            }
        }
        
        # But check if they have complete data - count records per samp_name
        if("samp_name" %in% names(existing_results) && "db_name" %in% names(existing_results)) {
            records_per_sample = existing_results %>%
                group_by(samp_name, db_name) %>%
                summarize(n_records = n(), .groups = "drop")
            
            # Check if any samp_names have suspiciously few records (less than expected metrics)
            # Expected: should have multiple metrics per sample
            suspicious_samples = records_per_sample %>%
                filter(n_records < 3)  # Less than 3 records suggests incomplete data
            
            if(nrow(suspicious_samples) > 0) {
                cat("  ⚠️  Found", nrow(suspicious_samples), "samp_names with suspiciously few records (< 3)\n")
                cat("  These may have incomplete data and should be reprocessed\n")
                
                # Get the files for these suspicious samples
                suspicious_samp_names = suspicious_samples$samp_name
                suspicious_basenames = paste0(suspicious_samp_names, "_scores.output")
                suspicious_files = filter_scores_list[filter_scores_basenames %in% suspicious_basenames]
                
                if(length(suspicious_files) > 0) {
                    cat("  ➕ Adding", length(suspicious_files), "suspicious files to processing queue for reprocessing\n")
                    new_files = c(new_files, suspicious_files)
                    new_files = unique(new_files)
                }
            }
        }
    }
} else {
    cat("  No existing results to compare against\n")
}

# Final check: If we still have no files but there are files that exist, something is wrong
cat("\n=== Final Summary ===\n")
cat("  Total _scores.output files found:", length(filter_scores_list), "\n")
cat("  Files marked as processed in log:", length(processed_basenames), "\n")
if(!is.null(existing_results)) {
    cat("  Unique samp_names in existing results:", length(existing_samp_names), "\n")
    cat("  Total records in existing results:", nrow(existing_results), "\n")
}
cat("  Files in processing queue:", length(new_files), "\n")

if(length(new_files) == 0) {
    cat("\n⚠️  WARNING: No files to process!\n")
    
    # Double-check: Are there files that should be processed?
    if(length(filter_scores_list) > 0) {
        cat("  But", length(filter_scores_list), "files exist - checking for discrepancies...\n")
        
        # AGGRESSIVE: If we have files but none in queue, process ALL files that aren't in processed log
        # This catches cases where files exist but diagnostic logic missed them
        if(length(processed_basenames) < length(filter_scores_list)) {
            cat("  ⚠️  Found", length(filter_scores_list) - length(processed_basenames), 
                "files not in processed log - adding ALL to queue\n")
            new_files = filter_scores_list[!filter_scores_basenames %in% processed_basenames]
            cat("  ➕ Added", length(new_files), "files to processing queue\n")
        }
        
        if(!is.null(existing_results) && length(existing_samp_names) > 0) {
            all_file_samp_names = unique(normalize_samp_name(sub("_scores.output$", "", filter_scores_basenames)))
            missing_count = length(setdiff(all_file_samp_names, existing_samp_names))
            
            cat("  Unique samp_names in files:", length(all_file_samp_names), "\n")
            cat("  Unique samp_names in results:", length(existing_samp_names), "\n")
            cat("  Missing samp_names:", missing_count, "\n")
            
            if(missing_count > 0) {
                cat("  ⚠️  DISCREPANCY DETECTED: Found", missing_count, "samp_names in files but not in results!\n")
                cat("  This suggests files need to be processed but were incorrectly filtered out.\n")
                cat("  FORCING processing of missing files...\n")
                
                # Force add all missing files
                missing_samp_names = setdiff(all_file_samp_names, existing_samp_names)
                missing_basenames = paste0(missing_samp_names, "_scores.output")
                missing_files = filter_scores_list[filter_scores_basenames %in% missing_basenames]
                new_files = missing_files
                cat("  ➕ Added", length(new_files), "missing files to processing queue\n")
            } else {
                # Even if all samp_names are present, check if counts match and data is complete
                cat("  All samp_names present, but checking record counts and data completeness...\n")
                if(length(all_file_samp_names) != length(existing_samp_names)) {
                    cat("  ⚠️  Count mismatch: files have", length(all_file_samp_names), "unique samp_names,")
                    cat("  but results have", length(existing_samp_names), "\n")
                }
                
                # Check if results have complete data for each samp_name
                if("samp_name" %in% names(existing_results) && "db_name" %in% names(existing_results)) {
                    # Count records per samp_name per database
                    records_per_sample = existing_results %>%
                        group_by(samp_name, db_name) %>%
                        summarize(n_records = n(), .groups = "drop")
                    
                    # Expected: each samp_name should have multiple metrics (at least 3-4)
                    incomplete_samples = records_per_sample %>%
                        filter(n_records < 3) %>%
                        group_by(samp_name) %>%
                        summarize(n_databases = n(), .groups = "drop")
                    
                    if(nrow(incomplete_samples) > 0) {
                        cat("  ⚠️  Found", nrow(incomplete_samples), "samp_names with incomplete data (< 3 records)\n")
                        cat("  These may need to be reprocessed\n")
                        
                        # Get files for incomplete samples
                        incomplete_samp_names = incomplete_samples$samp_name
                        incomplete_basenames = paste0(incomplete_samp_names, "_scores.output")
                        incomplete_files = filter_scores_list[filter_scores_basenames %in% incomplete_basenames]
                        
                        if(length(incomplete_files) > 0) {
                            cat("  ➕ Adding", length(incomplete_files), "incomplete files to processing queue\n")
                            new_files = c(new_files, incomplete_files)
                            new_files = unique(new_files)
                        }
                    }
                    
                    # Check for incorrectly assigned database names (sample IDs as db_name)
                    valid_databases = c("soil_microbe_db", "pluspf", "gtdb_207", "gtdb_207_unfiltered")
                    invalid_db_names = existing_results %>%
                        filter(!db_name %in% valid_databases) %>%
                        distinct(samp_name, db_name)
                    
                    if(nrow(invalid_db_names) > 0) {
                        cat("  ⚠️  Found", nrow(invalid_db_names), "samp_names with invalid db_name assignments:\n")
                        print(invalid_db_names)
                        cat("  These need to be reprocessed with correct database assignment\n")
                        
                        # Get files for these samples
                        invalid_samp_names = invalid_db_names$samp_name
                        invalid_basenames = paste0(invalid_samp_names, "_scores.output")
                        invalid_files = filter_scores_list[filter_scores_basenames %in% invalid_basenames]
                        
                        if(length(invalid_files) > 0) {
                            cat("  ➕ Adding", length(invalid_files), "files with invalid db_name to processing queue\n")
                            new_files = c(new_files, invalid_files)
                            new_files = unique(new_files)
                        }
                    }
                    
                    # Check database distribution
                    cat("  Records per database in results:\n")
                    db_counts = existing_results %>%
                        filter(db_name %in% valid_databases) %>%
                        group_by(db_name) %>%
                        summarize(
                            n_samples = length(unique(samp_name)),
                            n_records = n(),
                            .groups = "drop"
                        )
                    print(db_counts)
                    
                    # Check file distribution by database
                    cat("  Files per database:\n")
                    # Use full file list, not unique samp_names (since there are multiple files per samp_name)
                    file_db_map = tibble(
                        file_path = filter_scores_list,
                        basename = filter_scores_basenames
                    ) %>%
                        mutate(
                            samp_name = normalize_samp_name(sub("_scores.output$", "", basename)),
                            db_name = case_when(
                                grepl("soil_microbe_db", samp_name) ~ "soil_microbe_db",
                                grepl("pluspf", samp_name) ~ "pluspf",
                                grepl("gtdb_207_unfiltered", samp_name) ~ "gtdb_207_unfiltered",
                                grepl("gtdb", samp_name) ~ "gtdb_207",
                                TRUE ~ NA_character_
                            )
                        ) %>%
                        filter(!is.na(db_name))
                    
                    file_db_dist = file_db_map %>%
                        group_by(db_name) %>%
                        summarize(n_files = n(), .groups = "drop")
                    print(file_db_dist)
                    
                    # Compare file counts vs result counts per database
                    comparison = left_join(file_db_dist, db_counts, by = "db_name") %>%
                        mutate(
                            n_samples = ifelse(is.na(n_samples), 0, n_samples),
                            missing = n_files - n_samples
                        )
                    
                    missing_by_db = comparison %>% filter(missing > 0)
                    if(nrow(missing_by_db) > 0) {
                        cat("  ⚠️  Database count mismatches (files vs results):\n")
                        print(missing_by_db %>% select(db_name, n_files, n_samples, missing))
                        
                        # Find missing samp_names per database
                        for(i in 1:nrow(missing_by_db)) {
                            db = missing_by_db$db_name[i]
                            n_missing = missing_by_db$missing[i]
                            
                            # Get samp_names in files for this database
                            db_file_samp_names = file_db_map %>%
                                filter(db_name == db) %>%
                                pull(samp_name)
                            
                            # Get samp_names in results for this database
                            db_result_samp_names = existing_results %>%
                                filter(db_name == db) %>%
                                distinct(samp_name) %>%
                                pull(samp_name)
                            
                            # Find missing
                            db_missing_samp_names = setdiff(db_file_samp_names, db_result_samp_names)
                            
                            if(length(db_missing_samp_names) > 0) {
                                cat("  ➕ Found", length(db_missing_samp_names), "missing samp_names for", db, "\n")
                                if(length(db_missing_samp_names) <= 5) {
                                    cat("    Missing:", paste(db_missing_samp_names, collapse=", "), "\n")
                                } else {
                                    cat("    Examples:", paste(head(db_missing_samp_names, 3), collapse=", "), "...\n")
                                }
                                
                                db_missing_basenames = paste0(db_missing_samp_names, "_scores.output")
                                db_missing_files = filter_scores_list[filter_scores_basenames %in% db_missing_basenames]
                                
                                if(length(db_missing_files) > 0) {
                                    cat("  ➕ Adding", length(db_missing_files), "missing files for", db, "to processing queue\n")
                                    new_files = c(new_files, db_missing_files)
                                    new_files = unique(new_files)
                                } else {
                                    cat("  ⚠️  Could not find files for", length(db_missing_samp_names), "missing samp_names\n")
                                    cat("  Files may be in a location not being searched (e.g., on cluster)\n")
                                }
                            } else if(n_missing > 0) {
                                # Count mismatch but no missing samp_names found - this is suspicious
                                cat("  ⚠️  Count mismatch (", n_missing, "missing) but no missing samp_names found\n")
                                cat("  This suggests files may not be accessible or matching logic is incorrect\n")
                            }
                        }
                    }
                }
            }
        } else {
            cat("  ⚠️  No existing results to compare against\n")
        }
    }
}

# Final aggressive check: If we still have very few files but know there should be more
if(length(new_files) > 0 && length(new_files) < 50) {
    cat("\n⚠️  WARNING: Only", length(new_files), "files in processing queue, but diagnostic shows many missing samples\n")
    cat("  This suggests most files are not accessible locally.\n")
    cat("  Script should be run on the cluster where all files are located.\n")
    
    # AGGRESSIVE: Process ALL files that aren't in processed log
    # This ensures we catch any files that should be processed
    all_unprocessed = filter_scores_list[!filter_scores_basenames %in% processed_basenames]
    if(length(all_unprocessed) > length(new_files)) {
        cat("  ➕ Found", length(all_unprocessed), "total files not in processed log\n")
        cat("  Adding all", length(all_unprocessed), "unprocessed files to queue\n")
        new_files = all_unprocessed
        cat("  Total files to process now:", length(new_files), "\n")
    }
}

# Note: All file processing is handled in the first processing loop above
# This section is intentionally left empty to avoid duplicate processing

# Results have already been merged incrementally into output_file
# Just load the final file for final corrections
cat("Loading final results for corrections...\n")
if(file.exists(output_file)) {
    all_summaries = data.table::fread(output_file, showProgress = FALSE)
    cat("  Loaded", nrow(all_summaries), "total records\n")
} else {
    # Fallback: use new_summaries if file doesn't exist (shouldn't happen)
    all_summaries = new_summaries
    if(!inherits(all_summaries, "tbl")) {
        all_summaries = as_tibble(all_summaries)
    }
    cat("  Using new_summaries (", nrow(all_summaries), " records)\n")
}

# Remove duplicates
if(nrow(all_summaries) > 0) {
    key_cols = intersect(c("samp_name", "metric", "db_name", "sampleID"), names(all_summaries))
    if(length(key_cols) > 0) {
        cat("  Removing duplicates using key columns:", paste(key_cols, collapse=", "), "\n")
        all_summaries = all_summaries %>%
            distinct(across(all_of(key_cols)), .keep_all = TRUE)
        cat("  Records after deduplication:", nrow(all_summaries), "\n")
    }
}

# Final correction of database names in combined results
# Ensure gtdb_207_unfiltered and gtdb_207 are correctly assigned
if("samp_name" %in% names(all_summaries) && "db_name" %in% names(all_summaries)) {
    all_summaries = all_summaries %>%
        mutate(db_name = ifelse(grepl("gtdb_207_unfiltered", samp_name), "gtdb_207_unfiltered",
                               ifelse(db_name == "gtdb_207_filtered" & grepl("gtdb", samp_name), "gtdb_207", db_name)))
    cat("  Applied final database name corrections\n")
}

# Final correction of sampleID values - extract from samp_name to ensure correctness
# This fixes any cases where sampleID was incorrectly set to database names or other values
if("samp_name" %in% names(all_summaries) && "sampleID" %in% names(all_summaries)) {
    # Count how many have incorrect sampleID values before fixing
    invalid_sampleIDs_before = all_summaries %>%
        filter(sampleID %in% c("gtdb_207", "gtdb_207_unfiltered", "soil_microbe_db", "pluspf") |
               !grepl("^[A-Z]{4}_", sampleID))  # Sample IDs should start with 4 uppercase letters and underscore
    
    if(nrow(invalid_sampleIDs_before) > 0) {
        cat("  Found", nrow(invalid_sampleIDs_before), "rows with incorrect sampleID values - fixing...\n")
    } else {
        cat("  Checking sampleID values...\n")
    }
    
    # Always fix sampleID from samp_name to ensure correctness
    all_summaries = all_summaries %>%
        mutate(
            # Extract sampleID from samp_name (remove database suffix)
            # Database names are: gtdb_207, gtdb_207_unfiltered, soil_microbe_db, pluspf
            sampleID = sub("_(gtdb_207_unfiltered|gtdb_207|soil_microbe_db|pluspf)$", "", samp_name)
        )
    
    if(nrow(invalid_sampleIDs_before) > 0) {
        cat("  ✓ Fixed", nrow(invalid_sampleIDs_before), "rows with incorrect sampleID values\n")
    } else {
        cat("  ✓ All sampleID values are correct\n")
    }
}

# Save final corrected results (overwrites file with deduplicated and corrected version)
cat("Saving final corrected results to:", output_file, "\n")
cat("  Writing", nrow(all_summaries), "rows...\n")
data.table::fwrite(all_summaries, output_file)
cat("  ✓ File written successfully\n")

# Processed log is already updated incrementally during processing
# Just verify it's up to date
cat("Processed files log updated incrementally during processing\n")
cat("✓ Successfully processed", processed_count, "files\n")
cat("✓ Total records in output:", nrow(all_summaries), "\n")
cat("✓ Results saved to:", output_file, "\n")

# Diagnostic: Check for samples in lineage files that are missing from filter_results_summary
cat("\n=== Diagnostic: Checking for missing samples ===\n")
lineage_dir_new = "data/NEON_metagenome_classification/summary_files"
lineage_dir_old = "data/classification/taxonomic_rank_summaries/species"
databases_to_check = c("soil_microbe_db", "pluspf", "gtdb_207", "gtdb_207_unfiltered")

for(db in databases_to_check) {
    filename = paste0(db, "_species_merged_lineage.csv")
    lineage_file = NULL
    
    # Check new location
    if(file.exists(file.path(lineage_dir_new, filename))) {
        lineage_file = file.path(lineage_dir_new, filename)
    } else if(file.exists(file.path(lineage_dir_old, filename))) {
        lineage_file = file.path(lineage_dir_old, filename)
    }
    
    if(!is.null(lineage_file) && file.exists(lineage_file)) {
        tryCatch({
            lineage_df = read_csv(lineage_file, col_select = "sample_id", show_col_types = FALSE)
            if("sample_id" %in% names(lineage_df) && nrow(lineage_df) > 0) {
                lineage_samples = unique(lineage_df$sample_id)
                
                lineage_sampleIDs = sapply(lineage_samples, function(samp_name) {
                    # Remove database suffix to get sampleID
                    # Database names are: gtdb_207, gtdb_207_unfiltered, soil_microbe_db, pluspf
                    # Also remove _filtered if present (it's part of filename pattern, not database name)
                    samp_name = sub("_filtered$", "", samp_name)  # Remove _filtered first
                    sub("_(soil_microbe_db|pluspf|gtdb_207_unfiltered|gtdb_207)$", "", samp_name)
                })
                
                # Get samples in filter_results_summary for this database
                if(!is.null(all_summaries) && "sampleID" %in% names(all_summaries) && "db_name" %in% names(all_summaries)) {
                    summary_samples = all_summaries %>%
                        filter(db_name == db, metric == "percent_classified") %>%
                        distinct(sampleID) %>%
                        pull(sampleID)
                    
                    missing_samples = setdiff(unique(lineage_sampleIDs), summary_samples)
                    
                    if(length(missing_samples) > 0) {
                        cat(sprintf("\n  %s: %d samples in lineage but missing from filter_results_summary\n", db, length(missing_samples)))
                        
                        # Check if their _scores.output files exist
                        missing_with_scores = character(0)
                        missing_without_scores = character(0)
                        
                        for(samp_id in head(missing_samples, 20)) {  # Check first 20
                            # Try to find matching scores file
                            # samp_name format: {sampleID}-COMP_{db}
                            if(db == "gtdb_207_unfiltered") {
                                expected_samp_name = paste0(samp_id, "-COMP_gtdb_207_unfiltered")
                            } else if(db == "gtdb_207") {
                                expected_samp_name = paste0(samp_id, "-COMP_gtdb_207")
                            } else {
                                expected_samp_name = paste0(samp_id, "-COMP_", db)
                            }
                            
                            expected_scores_file = paste0(expected_samp_name, "_scores.output")
                            
                            # Check all directories
                            scores_found = FALSE
                            for(dir_path in dirs_to_search) {
                                scores_path = file.path(dir_path, expected_scores_file)
                                if(file.exists(scores_path)) {
                                    missing_with_scores = c(missing_with_scores, samp_id)
                                    scores_found = TRUE
                                    break
                                }
                            }
                            
                            if(!scores_found) {
                                missing_without_scores = c(missing_without_scores, samp_id)
                            }
                        }
                        
                        if(length(missing_with_scores) > 0) {
                            cat(sprintf("    ⚠️  %d missing samples have _scores.output files that weren't processed:\n", length(missing_with_scores)))
                            cat(sprintf("      Examples: %s\n", paste(head(missing_with_scores, 5), collapse = ", ")))
                            if(length(missing_with_scores) > 5) cat(sprintf("      ... and %d more\n", length(missing_with_scores) - 5))
                            
                            # Find and add all missing files with scores (not just first 20)
                            cat("    Searching for all missing files with scores...\n")
                            all_missing_with_scores = character(0)
                            for(samp_id in missing_samples) {
                                if(db == "gtdb_207_unfiltered") {
                                    expected_samp_name = paste0(samp_id, "-COMP_gtdb_207_unfiltered")
                                } else if(db == "gtdb_207") {
                                    expected_samp_name = paste0(samp_id, "-COMP_gtdb_207")
                                } else {
                                    expected_samp_name = paste0(samp_id, "-COMP_", db)
                                }
                                expected_scores_file = paste0(expected_samp_name, "_scores.output")
                                
                                for(dir_path in dirs_to_search) {
                                    scores_path = file.path(dir_path, expected_scores_file)
                                    if(file.exists(scores_path)) {
                                        all_missing_with_scores = c(all_missing_with_scores, scores_path)
                                        break
                                    }
                                }
                            }
                            
                            if(length(all_missing_with_scores) > 0) {
                                cat(sprintf("    ➕ Found %d missing files with scores - these should be processed in next run\n", length(all_missing_with_scores)))
                                # Note: We can't add them to new_files here because processing has already completed
                                # But we can write them to a file for the next run, or the user can re-run the script
                            }
                        }
                        
                        if(length(missing_without_scores) > 0) {
                            cat(sprintf("    ⚠️  %d missing samples have NO _scores.output files (may have been deleted):\n", length(missing_without_scores)))
                            cat(sprintf("      Examples: %s\n", paste(head(missing_without_scores, 5), collapse = ", ")))
                            if(length(missing_without_scores) > 5) cat(sprintf("      ... and %d more\n", length(missing_without_scores) - 5))
                        }
                        
                        if(length(missing_samples) > 20) {
                            cat(sprintf("    (Only checked first 20 of %d missing samples for file existence)\n", length(missing_samples)))
                        }
                    } else {
                        cat(sprintf("  %s: ✓ All lineage samples are in filter_results_summary\n", db))
                    }
                }
            }
        }, error = function(e) {
            cat(sprintf("  %s: Error reading lineage file: %s\n", db, e$message))
        })
    }
}
cat("\n")
}  # End of if(length(new_files) > 0) block
