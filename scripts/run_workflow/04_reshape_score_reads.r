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

# Helper function to extract database name from samp_name
extract_db_name = function(samp_name) {
    case_when(
        grepl("gtdb_207_unfiltered", samp_name) ~ "gtdb_207_unfiltered",
        grepl("gtdb_207", samp_name) ~ "gtdb_207",
        grepl("soil_microbe_db", samp_name) ~ "soil_microbe_db",
        grepl("pluspf", samp_name) ~ "pluspf",
        TRUE ~ NA_character_
    )
}

# Function to aggregate very large files using command-line tools (for files > 2GB)
# Uses awk to calculate statistics without loading entire file into R
# Returns summary dataframe or NULL if error
aggregate_with_cli = function(file_path, samp_name, seq_depth_df,
                              max_entropy = 0.1, max_multiplicity = 2, min_consistency = 0.9) {
    samp_name = normalize_samp_name(samp_name)
    
    tryCatch({
        # First, check if awk is available
        awk_check = system2("which", "awk", stdout = TRUE, stderr = TRUE)
        if(length(awk_check) == 0 || awk_check == "") {
            cat("    ERROR: awk not found - cannot use CLI aggregation\n")
            return(NULL)
        }
        
        # Read header to detect delimiter and get column positions
        header_line = readLines(file_path, n = 1)
        
        # Detect delimiter (comma or tab)
        if(grepl(",", header_line) && !grepl("\t", header_line)) {
            delimiter = ","
            FS_char = ","
        } else if(grepl("\t", header_line)) {
            delimiter = "\t"
            FS_char = "\t"
        } else {
            # Try comma first, then tab
            delimiter = ","
            FS_char = ","
        }
        
        header_cols = strsplit(header_line, delimiter)[[1]]
        
        # Find column indices (0-indexed for awk)
        n_kmers_idx = which(header_cols == "n_kmers") - 1
        consistency_idx = which(header_cols == "consistency") - 1
        multiplicity_idx = which(header_cols == "multiplicity") - 1
        entropy_idx = which(header_cols == "entropy") - 1
        confidence_idx = which(header_cols == "confidence") - 1
        
        # Check all required columns are found
        missing_cols = character(0)
        if(length(n_kmers_idx) == 0 || n_kmers_idx < 0) missing_cols = c(missing_cols, "n_kmers")
        if(length(consistency_idx) == 0 || consistency_idx < 0) missing_cols = c(missing_cols, "consistency")
        if(length(multiplicity_idx) == 0 || multiplicity_idx < 0) missing_cols = c(missing_cols, "multiplicity")
        if(length(entropy_idx) == 0 || entropy_idx < 0) missing_cols = c(missing_cols, "entropy")
        if(length(confidence_idx) == 0 || confidence_idx < 0) missing_cols = c(missing_cols, "confidence")
        
        if(length(missing_cols) > 0) {
            cat("    ERROR: Required columns not found in header:", paste(missing_cols, collapse=", "), "\n")
            cat("    Header columns found:", paste(header_cols, collapse=", "), "\n")
            return(NULL)
        }
        
        # Use awk to calculate aggregates in one pass
        # Strategy: Calculate sums and counts without loading data into R
        # Note: awk uses 1-based indexing for fields
        n_kmers_col = n_kmers_idx + 1
        consistency_col = consistency_idx + 1
        multiplicity_col = multiplicity_idx + 1
        entropy_col = entropy_idx + 1
        confidence_col = confidence_idx + 1
        
        # Escape delimiter for awk if it's a comma (need to use -F flag instead)
        if(FS_char == ",") {
            FS_flag = "-F,"
        } else {
            FS_flag = ""
        }
        
        # Build awk script - if using -F flag, don't set FS in BEGIN (let -F handle it)
        # Build awk script using paste to avoid sprintf argument counting issues
        if(FS_flag != "") {
            # Using -F flag, so don't set FS in BEGIN block
            awk_script = paste0('
BEGIN {
    sum_consistency = 0
    sum_multiplicity = 0
    sum_entropy = 0
    sum_confidence = 0
    n_valid = 0
    n_passing = 0
}
')
        } else {
            # Not using -F flag, so set FS in BEGIN block
            awk_script = paste0('
BEGIN {
    FS = "', FS_char, '"
    sum_consistency = 0
    sum_multiplicity = 0
    sum_entropy = 0
    sum_confidence = 0
    n_valid = 0
    n_passing = 0
}
')
        }
        
        # Append the rest of the script
        awk_script = paste0(awk_script, '
NR > 1 && $', n_kmers_col, ' > 0 {
    sum_consistency += $', consistency_col, '
    sum_multiplicity += $', multiplicity_col, '
    sum_entropy += $', entropy_col, '
    sum_confidence += $', confidence_col, '
    n_valid++
    if ($', multiplicity_col, ' <= ', max_multiplicity, ' && $', consistency_col, ' >= ', min_consistency, ' && $', entropy_col, ' <= ', max_entropy, ') {
        n_passing++
    }
}
END {
    if (n_valid > 0) {
        printf "%.10f %.10f %.10f %.10f %d %d\\n", 
            sum_consistency/n_valid, 
            sum_multiplicity/n_valid, 
            sum_entropy/n_valid, 
            sum_confidence/n_valid,
            n_valid,
            n_passing
    } else {
        printf "0 0 0 0 0 0\\n"
    }
}')
        
        # Write awk script to temp file
        awk_file = tempfile(fileext = ".awk")
        writeLines(awk_script, awk_file)
        
        # Run awk with appropriate field separator flag
        # Note: system2 with stderr=TRUE combines stdout and stderr
        # We'll capture stderr separately to see actual errors
        stderr_file = tempfile(fileext = ".stderr")
        
        if(FS_flag != "") {
            awk_args = c(FS_flag, "-f", awk_file, file_path)
        } else {
            awk_args = c("-f", awk_file, file_path)
        }
        
        # Run awk and capture stderr separately
        awk_result = tryCatch({
            system2("awk", awk_args, stdout = TRUE, stderr = stderr_file, wait = TRUE)
        }, error = function(e) {
            cat("    ERROR running awk:", e$message, "\n")
            return(NULL)
        })
        
        # Check for errors in stderr
        if(file.exists(stderr_file) && file.info(stderr_file)$size > 0) {
            stderr_content = readLines(stderr_file, warn = FALSE)
            if(length(stderr_content) > 0 && any(nchar(stderr_content) > 0)) {
                cat("    awk stderr:", paste(stderr_content, collapse="; "), "\n")
            }
        }
        unlink(stderr_file)
        
        # Clean up
        unlink(awk_file)
        
        if(is.null(awk_result) || length(awk_result) == 0) {
            cat("    ERROR: awk aggregation failed (no output)\n")
            return(NULL)
        }
        
        if(length(awk_result) > 1) {
            # If stderr was captured, it might be in the result
            # Check if any line looks like an error
            error_lines = grep("error|Error|ERROR|warning|Warning", awk_result, ignore.case = TRUE, value = TRUE)
            if(length(error_lines) > 0) {
                cat("    ERROR: awk reported:", paste(error_lines, collapse="; "), "\n")
                return(NULL)
            }
            # If multiple lines but no errors, take the last one (should be the printf output)
            awk_result = awk_result[length(awk_result)]
        }
        
        # Parse results
        stats = as.numeric(strsplit(awk_result, " ")[[1]])
        if(length(stats) != 6) {
            cat("    ERROR: Unexpected awk output format\n")
            return(NULL)
        }
        
        mean_consistency = stats[1]
        mean_multiplicity = stats[2]
        mean_entropy = stats[3]
        mean_confidence = stats[4]
        n_valid_rows = as.integer(stats[5])
        n_passing = as.integer(stats[6])
        
        if(n_valid_rows == 0) {
            return(NULL)
        }
        
        # Extract sampleID and db_name
        db_name_val = extract_db_name(samp_name)
        sampleID = sub("_(gtdb_207_unfiltered|gtdb_207|soil_microbe_db|pluspf)$", "", samp_name)
        
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
        
        # Create summary dataframe
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
        cat("    ERROR in CLI aggregation:", e$message, "\n")
        return(NULL)
    })
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
        
        # Use file connection for efficient chunked reading (avoids re-scanning from top)
        # Only load required numeric columns to avoid 2^31-1 byte limit on string columns
        required_cols = c("n_kmers", "consistency", "multiplicity", "entropy", "confidence")
        
        # Read header to get column names and detect delimiter
        header_line = readLines(file_path, n = 1)
        
        # Detect delimiter
        if(grepl(",", header_line) && !grepl("\t", header_line)) {
            delimiter = ","
        } else if(grepl("\t", header_line)) {
            delimiter = "\t"
        } else {
            delimiter = ","
        }
        
        header_cols = strsplit(header_line, delimiter)[[1]]
        
        # Check required columns exist
        missing_cols = setdiff(required_cols, header_cols)
        if(length(missing_cols) > 0) {
            cat("    ERROR: Missing required columns:", paste(missing_cols, collapse=", "), "\n")
            return(NULL)
        }
        
        # Open file connection for efficient sequential reading
        con = file(file_path, "r")
        on.exit(close(con), add = TRUE)
        
        # Read and skip header
        readLines(con, n = 1)
        
        # Process chunks using connection (much faster than skip-based reading)
        chunk_num = 0
        lines_processed = 0
        
        while(TRUE) {
            chunk_num = chunk_num + 1
            
            # Read chunk of lines
            chunk_lines = readLines(con, n = chunk_size, warn = FALSE)
            
            if(length(chunk_lines) == 0) {
                break  # End of file
            }
            
            # Process chunk using fread on text (only load required numeric columns)
            chunk_df = NULL
            tryCatch({
                # Create temporary text with header + chunk lines
                chunk_text = c(header_line, chunk_lines)
                
                # Use fread with select to only load numeric columns (avoids 2^31-1 byte limit)
                chunk_df = suppressWarnings(data.table::fread(
                    text = chunk_text,
                    select = required_cols,  # Only load needed columns
                    colClasses = list(
                        numeric = c("n_kmers", "consistency", "multiplicity", "entropy", "confidence")
                    ),
                    showProgress = FALSE
                ))
            }, error = function(e) {
                if(grepl("character strings are limited to 2\\^31-1 bytes", e$message)) {
                    cat("    ERROR: File contains lines too long for R to process at chunk", chunk_num, "\n")
                    cat("    This usually means a single column value exceeds 2GB - file may be corrupted\n")
                } else {
                    cat("    ERROR reading chunk", chunk_num, ":", e$message, "\n")
                }
                chunk_df <<- NULL
            })
            
            if(is.null(chunk_df) || nrow(chunk_df) == 0) {
                break  # End of file or error
            }
            
            # Aggregate statistics from this chunk
            valid_chunk = chunk_df[n_kmers > 0]
            if(nrow(valid_chunk) > 0) {
                sum_consistency = sum_consistency + sum(valid_chunk$consistency, na.rm = TRUE)
                sum_multiplicity = sum_multiplicity + sum(valid_chunk$multiplicity, na.rm = TRUE)
                sum_entropy = sum_entropy + sum(valid_chunk$entropy, na.rm = TRUE)
                sum_confidence = sum_confidence + sum(valid_chunk$confidence, na.rm = TRUE)
                n_valid_rows = n_valid_rows + nrow(valid_chunk)
                
                # Count passing reads (using data.table syntax for speed)
                n_passing = n_passing + nrow(valid_chunk[
                    multiplicity <= max_multiplicity & 
                    consistency >= min_consistency & 
                    entropy <= max_entropy
                ])
            }
            
            total_rows = total_rows + nrow(chunk_df)
            chunk_lines_length = length(chunk_lines)
            lines_processed = lines_processed + chunk_lines_length
            
            # Clean up chunk from memory
            rm(chunk_df, valid_chunk, chunk_lines, chunk_text)
            gc(verbose = FALSE)
            
            if(chunk_lines_length < chunk_size) {
                break  # Last chunk
            }
            
            if(chunk_num %% 20 == 0) {
                cat("    Processed", chunk_num, "chunks (", format(total_rows, big.mark=","), "rows)...\n")
            }
        }  # end while
        
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
    # Use data.table::melt instead of pivot_longer for better memory efficiency
    joined_df = inner_join(score_summary_df_unique, filter_summary_df_unique, by=join_by(samp_name))
    
    # Convert to data.table for efficient melting
    joined_dt = as.data.table(joined_df)
    
    # Identify metric columns (all columns except db_name, sampleID, samp_name)
    id_cols = c("db_name", "sampleID", "samp_name")
    measure_cols = setdiff(names(joined_dt), id_cols)
    
    # Use data.table::melt (more memory-efficient than pivot_longer)
    out_df = melt(
        joined_dt,
        id.vars = id_cols,
        measure.vars = measure_cols,
        variable.name = "metric",
        value.name = "value",
        na.rm = TRUE
    ) %>%
        distinct() %>%
        as_tibble()
    
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
# STEP 1: Build inventory of expected vs actual state
# ============================================================================
cat("\n=== Building file inventory ===\n")

# Load lineage files to determine what samples SHOULD exist
# NOTE: This is for DIAGNOSTICS ONLY - we process ALL _scores.output files that exist,
# regardless of whether they're in lineage files. Some samples may have _scores.output
# files but not be in lineage files yet (e.g., if script 3 skipped them due to missing .b2 files)
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
                
                # After script 3 normalization, sample_id format is: {sampleID}-COMP_{db}_filtered
                # Actual _scores.output files have format: {sampleID}-COMP_{db}_scores.output
                # So we just need to remove _filtered to get the expected samp_name
                expected_samp_names = sub("_filtered$", "", lineage_samples)
                # Normalize to remove any trailing underscores
                expected_samp_names = normalize_samp_name(expected_samp_names)
                
                # Extract sampleID for reference (remove database suffix and -COMP)
                lineage_sampleIDs = sapply(expected_samp_names, function(samp_name) {
                    # Remove database suffix (handle both -COMP_ and _ formats)
                    base_id = sub("[-_]COMP?_(soil_microbe_db|pluspf|gtdb_207_unfiltered|gtdb_207)$", "", samp_name)
                    # Also try removing just the database suffix if -COMP_ removal didn't work
                    if(base_id == samp_name) {
                        base_id = sub("_(soil_microbe_db|pluspf|gtdb_207_unfiltered|gtdb_207)$", "", samp_name)
                    }
                    return(base_id)
                })
                
                expected_samples = bind_rows(
                    expected_samples,
                    tibble(
                        sampleID = lineage_sampleIDs,
                        samp_name = expected_samp_names,
                        db_name = db,
                        expected_scores_file = paste0(expected_samp_names, "_scores.output"),
                        expected_scores_file_alt = paste0(expected_samp_names, "__scores.output")
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
    samp_name = normalize_samp_name(sub("_+scores.output$", "", filter_scores_basenames))
) %>%
    mutate(db_name = extract_db_name(samp_name))

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
    samp_name = normalize_samp_name(sub("_+scores.output$", "", processed_basenames))
) %>%
    mutate(db_name = extract_db_name(samp_name))

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
# STEP 2: Make processing decisions
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
new_files_count = sum(!files_to_process$in_results, na.rm = TRUE)
wrong_db_count = sum(files_to_process$in_results & !is.na(files_to_process$db_name) & 
                     !is.na(files_to_process$results_db_name) & 
                     files_to_process$db_name != files_to_process$results_db_name, na.rm = TRUE)

cat("  Files to process (not in results):", new_files_count, "\n")
if(wrong_db_count > 0) {
    cat("  Files with wrong db_name in results:", wrong_db_count, " (will fix in results)\n")
}

# Extract files to process - only files that are NOT in results
# Files with wrong db_name will be fixed in existing results, not reprocessed
# IMPORTANT: Process ALL files that exist and aren't in results, regardless of lineage files
# Lineage files are only for diagnostics - some samples may have _scores.output files
# but not be in lineage files yet (e.g., if script 3 skipped them due to missing .b2 files)
new_files = files_to_process %>%
    filter(needs_processing & !in_results) %>%
    pull(file_path)

cat("\n=== Processing queue ===\n")
cat("  Files to process:", length(new_files), "\n")

# Diagnostic: Compare with expected samples from lineage (for info only)
if(nrow(expected_samples) > 0) {
    cat("\n  Diagnostic: Expected samples from lineage:", nrow(expected_samples), "\n")
    expected_samp_names = expected_samples$samp_name
    found_samp_names = unique(files_to_process$samp_name)
    
    missing_expected = setdiff(expected_samp_names, found_samp_names)
    found_but_not_expected = setdiff(found_samp_names, expected_samp_names)
    
    if(length(missing_expected) > 0) {
        cat("    ⚠️  ", length(missing_expected), "expected samples from lineage don't have _scores.output files\n")
    }
    if(length(found_but_not_expected) > 0) {
        cat("    ℹ️  ", length(found_but_not_expected), "files exist but aren't in lineage (will still be processed)\n")
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
        
        # Check if file is readable
        if(!file.access(file_path, mode = 4) == 0) {
            return(list(valid = FALSE, reason = "File is not readable"))
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
        
        # Determine processing strategy based on file size
        # Note: Very large files (> 2 GB) will use command-line aggregation strategy
        # instead of R's chunked reading to avoid character string limits
        file_size_gb = file_info$size / 1024^3
        use_cli = file_size_gb > 2.0
        
        return(list(valid = TRUE, reason = "", use_cli_aggregation = use_cli))
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
    skipped_count = 0
    processed_files_list = character(0)
    skipped_files_list = character(0)  # Track skipped files to log them
    skip_reasons = character(0)  # Track why files were skipped
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
                skipped_files_list = c(skipped_files_list, normalizePath(file_path))
                skip_reasons = c(skip_reasons, validation$reason)
                skipped_count = skipped_count + 1
                next
            }
            
            # Read in the file with appropriate strategy based on size
            samp_name_from_file = normalize_samp_name(sub("_+scores.output", "", basename(file_path)))
            file_size_gb = file.info(file_path)$size / 1024^3
            file_summary = NULL
            
            # Strategy selection based on file size:
            # - > 2GB: Use command-line aggregation (awk) to avoid R memory limits
            # - > 500MB: Use R chunked reading
            # - < 500MB: Try normal reading first, fall back to chunked if needed
            if(validation$use_cli_aggregation) {
                # Very large file - use command-line aggregation
                cat("  Processing very large file", basename(file_path), sprintf("(%.2f GB) with CLI aggregation...\n", file_size_gb))
                file_summary = aggregate_with_cli(file_path, samp_name_from_file, seq_depth_df)
                if(is.null(file_summary)) {
                    reason = "CLI aggregation failed"
                    cat("WARNING: Skipping", basename(file_path), "-", reason, "\n")
                    skipped_files_list = c(skipped_files_list, normalizePath(file_path))
                    skip_reasons = c(skip_reasons, reason)
                    skipped_count = skipped_count + 1
                    next
                }
            } else if(file_size_gb > 0.5) {
                # Large file - use chunked reading
                cat("  Reading large file", basename(file_path), sprintf("(%.2f GB) in chunks...\n", file_size_gb))
                file_summary = read_and_summarize_chunked(file_path, samp_name_from_file, seq_depth_df)
                if(is.null(file_summary)) {
                    reason = "failed to process with chunked reading"
                    cat("WARNING: Skipping", basename(file_path), "-", reason, "\n")
                    skipped_files_list = c(skipped_files_list, normalizePath(file_path))
                    skip_reasons = c(skip_reasons, reason)
                    skipped_count = skipped_count + 1
                    next
                }
            } else {
                # Small file - try normal reading first (only load required columns to avoid 2^31-1 limit)
                tryCatch({
                    # Only load numeric columns needed for aggregation (avoids loading large string columns)
                    required_cols = c("n_kmers", "consistency", "multiplicity", "entropy", "confidence")
                    df_in = suppressWarnings(data.table::fread(
                        file_path, 
                        select = required_cols,  # Only load needed columns
                        colClasses = list(
                            numeric = required_cols
                        ),
                        showProgress = FALSE
                    ))
                    if(is.null(df_in) || nrow(df_in) == 0) {
                        reason = "file is empty or has no data rows"
                        cat("WARNING: Skipping", basename(file_path), "-", reason, "\n")
                        skipped_files_list = c(skipped_files_list, normalizePath(file_path))
                        skip_reasons = c(skip_reasons, reason)
                        skipped_count = skipped_count + 1
                        next
                    }
                    df_in$samp_name = samp_name_from_file
                    file_summary <- summarize_filter_scores(df_in, seq_depth_df)
                    
                    # Clean up immediately
                    rm(df_in)
                    gc(verbose = FALSE)
                }, error = function(e) {
                    if(grepl("character strings are limited to 2\\^31-1 bytes", e$message)) {
                        # Fall back to chunked reading
                        cat("  File hit character limit, switching to chunked reading...\n")
                        file_summary <<- read_and_summarize_chunked(file_path, samp_name_from_file, seq_depth_df)
                        if(is.null(file_summary)) {
                            reason = "failed to process even with chunked reading"
                            cat("WARNING: Skipping", basename(file_path), "-", reason, "\n")
                            skipped_files_list = c(skipped_files_list, normalizePath(file_path))
                            skip_reasons = c(skip_reasons, reason)
                            skipped_count = skipped_count + 1
                        }
                    } else {
                        cat("WARNING: Error reading", basename(file_path), ":", e$message, "\n")
                    }
                })
            }
            
            if(is.null(file_summary) || nrow(file_summary) == 0) {
                reason = "no valid data extracted from file"
                cat("WARNING: Skipping", basename(file_path), "-", reason, "\n")
                skipped_files_list = c(skipped_files_list, normalizePath(file_path))
                skip_reasons = c(skip_reasons, reason)
                skipped_count = skipped_count + 1
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
            
            # Clean up file_summary from memory
            rm(file_summary)
            
            # Force garbage collection periodically (every 10 files or after large files)
            if(processed_count %% 10 == 0 || file_size_gb > 0.5) {
                gc(verbose = FALSE)
            }
            
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
                # Load existing processed files from log
                existing_processed = if(file.exists(processed_log_file)) {
                    readLines(processed_log_file)
                } else {
                    character(0)
                }
                all_processed = unique(c(existing_processed, processed_files_list))
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
    if(skipped_count > 0) {
        cat("Skipped", skipped_count, "files (empty, too large, or failed to process)\n")
        # Add skipped files to processed log so they don't get re-queued
        if(length(skipped_files_list) > 0) {
            # Load existing processed files from log
            existing_processed = if(file.exists(processed_log_file)) {
                readLines(processed_log_file)
            } else {
                character(0)
            }
            all_processed = unique(c(existing_processed, processed_files_list, skipped_files_list))
            writeLines(all_processed, processed_log_file)
            cat("  Added skipped files to processed log to prevent re-queuing\n")
            
            # Show summary of skip reasons
            skip_summary = table(skip_reasons)
            cat("  Skip reasons:\n")
            for(i in 1:length(skip_summary)) {
                cat(sprintf("    %s: %d files\n", names(skip_summary)[i], skip_summary[i]))
            }
        }
    }
    
    # Load final results for compatibility with rest of script
    if(file.exists(output_file)) {
        new_summaries = data.table::fread(output_file, showProgress = FALSE)
        # Filter to only newly processed samples
        if(length(processed_files_list) > 0) {
            processed_samp_names = normalize_samp_name(sub("_+scores.output$", "", basename(processed_files_list)))
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


# Simple diagnostic summary
cat("\n=== Summary ===\n")
cat("  Total _scores.output files found:", length(filter_scores_list), "\n")
if(!is.null(existing_results) && length(existing_samp_names) > 0) {
    all_file_samp_names = unique(normalize_samp_name(sub("_+scores.output$", "", filter_scores_basenames)))
    missing_count = length(setdiff(all_file_samp_names, existing_samp_names))
    cat("  Files already in results:", sum(files_to_process$in_results, na.rm = TRUE), "\n")
    cat("  Files to process:", length(new_files), "\n")
    if(missing_count > 0) {
        cat("  Missing samp_names:", missing_count, "\n")
    }
} else {
    cat("  Files to process:", length(new_files), "\n")
}

# Final safety check: If no files to process but files exist, process all unprocessed
if(length(new_files) == 0 && length(filter_scores_list) > 0) {
    cat("\n⚠️  No files in queue, but", length(filter_scores_list), "files exist\n")
    cat("  Processing all files not in processed log...\n")
    new_files = filter_scores_list[!filter_scores_basenames %in% processed_basenames]
    cat("  ➕ Added", length(new_files), "files to processing queue\n")
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

# Helper function to find scores file for a sampleID
find_scores_file = function(samp_id, db_name, search_dirs) {
    # Normalize samp_id: remove trailing underscores and any database suffix
    samp_id_clean = sub("_+$", "", samp_id)
    samp_id_clean = sub("_(gtdb_207_unfiltered|gtdb_207|soil_microbe_db|pluspf)_*$", "", samp_id_clean)
    
    # Construct base pattern (sampleID should be like "BARR_002-O-20170808")
    # Try multiple patterns to account for variations
    patterns = c(
        paste0("^", samp_id_clean, "-COMP_", db_name, "_scores.output$"),
        paste0("^", samp_id_clean, "-COMP_", db_name, "__scores.output$"),
        paste0("^", samp_id_clean, "-COMP_", db_name, "_+scores.output$")
    )
    
    # Also try if samp_id already includes -COMP
    if(grepl("-COMP", samp_id_clean)) {
        base_name = sub("-COMP.*$", "", samp_id_clean)
        patterns = c(patterns,
            paste0("^", samp_id_clean, "_scores.output$"),
            paste0("^", samp_id_clean, "__scores.output$"),
            paste0("^", samp_id_clean, "_+scores.output$")
        )
    }
    
    # Search all directories with pattern matching
    for(dir_path in search_dirs) {
        if(!dir.exists(dir_path)) next
        
        # List all scores files in this directory
        all_scores = list.files(dir_path, pattern = "_scores.output$", full.names = TRUE)
        
        # Try each pattern
        for(pattern in patterns) {
            matching = grep(pattern, basename(all_scores), value = TRUE)
            if(length(matching) > 0) {
                return(file.path(dir_path, matching[1]))
            }
        }
        
        # Fallback: try fuzzy matching on sampleID (first part before -COMP or _)
        samp_base = sub("[-_].*$", "", samp_id_clean)
        fuzzy_pattern = paste0("^", samp_base, ".*", db_name, ".*_scores.output$")
        matching = grep(fuzzy_pattern, basename(all_scores), value = TRUE, ignore.case = TRUE)
        if(length(matching) > 0) {
            # Check if it's a reasonable match (contains the sampleID pattern)
            for(match_file in matching) {
                if(grepl(samp_id_clean, match_file) || grepl(samp_base, match_file)) {
                    return(file.path(dir_path, match_file))
                }
            }
        }
    }
    
    return(NULL)
}

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
                    # After script 3 normalization, all sample_ids have _filtered (rank suffixes removed)
                    # Remove _filtered (always present after normalization) then database suffix
                    samp_name = sub("_filtered$", "", samp_name)  # Remove _filtered first
                    # Remove trailing underscores
                    samp_name = sub("_+$", "", samp_name)
                    # Remove database suffix (handle both with and without trailing underscores)
                    samp_name = sub("_(soil_microbe_db|pluspf|gtdb_207_unfiltered|gtdb_207)_*$", "", samp_name)
                    # Remove any remaining trailing underscores
                    sub("_+$", "", samp_name)
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
                            scores_file = find_scores_file(samp_id, db, dirs_to_search)
                            if(!is.null(scores_file) && file.exists(scores_file)) {
                                missing_with_scores = c(missing_with_scores, samp_id)
                            } else {
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
                                scores_file = find_scores_file(samp_id, db, dirs_to_search)
                                if(!is.null(scores_file) && file.exists(scores_file)) {
                                    all_missing_with_scores = c(all_missing_with_scores, scores_file)
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