#!/usr/bin/env Rscript
# Extract and reshape scoring information from Architeuthis _scores.output files
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

# Function to summarize the "scores" output from Architeuthis
# Format is one row per classified read
summarize_filter_scores = function(scores_df, seq_depth_df, max_entropy = .1,
                                   max_multiplicity = 2,
                                   min_consistency = .9) {
    
    scores_df = scores_df %>% mutate(db_name = ifelse(grepl("soil_microbe_db", samp_name), "soil_microbe_db",
                                                      ifelse(grepl("pluspf", samp_name), "pluspf",
                                                             ifelse(grepl("unfiltered", samp_name), "gtdb_207_unfiltered",
                                                                    ifelse(grepl("gtdb", samp_name), "gtdb_207_filtered", NA
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
        add_tally(name = "n_total_classified_reads") %>%
        group_by(db_name, samp_name, pass_filter) %>%
        add_tally(name = "n_reads") %>%
        group_by(db_name, samp_name, pass_filter,n_total_classified_reads,n_reads) %>%
        reframe(pct_of_classified_passing = n_reads/n_total_classified_reads) %>% 
        filter(pass_filter == 1) %>% ungroup %>%
        select(-pass_filter) %>% distinct() %>% 
        
        # Parse sample IDs into grouping information
        separate(samp_name,
                 into = c("sampleID","db_name"),
                 sep = "COMP_", remove = F, extra = "merge") %>%
        mutate(sampleID = str_remove(samp_name, paste0("_",db_name))) %>% distinct
    
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
seq_depth_df <- readRDS(seq_depth_file) %>% 
    select(-c(db_name, identified_reads))

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

# Load list of already processed files
processed_files = character(0)
if(file.exists(processed_log_file)) {
    processed_files = readLines(processed_log_file)
    cat("Found", length(processed_files), "already processed files\n")
}

# Identify new files to process
new_files = setdiff(filter_scores_list, processed_files)
cat("Found", length(new_files), "new files to process\n")

if(length(new_files) == 0) {
    cat("✓ All files have already been processed. Exiting.\n")
    quit(status = 0)
}

# Load existing results if output file exists
existing_results = NULL
if(file.exists(output_file)) {
    cat("Loading existing results from:", output_file, "\n")
    
    # Check file size first
    file_size_mb = file.info(output_file)$size / 1024 / 1024
    cat("  File size:", round(file_size_mb, 2), "MB\n")
    
    if(file_size_mb > 500) {
        cat("  WARNING: Large file detected. Using data.table for faster reading...\n")
        existing_results = data.table::fread(output_file, showProgress = TRUE)
        existing_results = as_tibble(existing_results)
    } else {
        existing_results = read_csv(output_file, show_col_types = FALSE)
    }
    
    cat("Found", nrow(existing_results), "existing records\n")
    
    # Get list of samp_names already in results
    existing_samp_names = unique(existing_results$samp_name)
    
    # Filter out files that are already in results
    new_files = new_files[!basename(new_files) %in% paste0(existing_samp_names, "_scores.output")]
    cat("After checking existing results:", length(new_files), "new files to process\n")
}

if(length(new_files) == 0) {
    cat("✓ All files have already been processed. Exiting.\n")
    quit(status = 0)
}

# Process new files
cat("Processing", length(new_files), "new scoring files...\n")
cat("This may take several minutes...\n")

# Set up parallel processing
plan(multisession, workers = min(18, length(new_files)))

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
summaries_list = future_lapply(new_files, function(x) {
    tryCatch({
        # Validate file before reading
        validation = validate_file(x)
        if(!validation$valid) {
            cat("WARNING: Skipping", basename(x), "-", validation$reason, "\n")
            return(NULL)
        }
        
        # Read in the (large) files with warning suppression
        # Suppress warnings about missing newlines - we handle validation separately
        df_in = suppressWarnings(data.table::fread(x, showProgress = FALSE))
        
        # Check if file was read successfully
        if(is.null(df_in) || nrow(df_in) == 0) {
            cat("WARNING: File", basename(x), "is empty or could not be read\n")
            return(NULL)
        }
        
        df_in$samp_name = sub("_scores.output", "", basename(x))
        out <- summarize_filter_scores(df_in, seq_depth_df)
        return(out)
    }, error = function(e) {
        cat("ERROR processing", basename(x), ":", e$message, "\n")
        return(NULL)
    })
})

end_time = Sys.time()
cat("Processing completed in", round(as.numeric(difftime(end_time, start_time, units = "secs")), 2), "seconds\n")

# Remove NULL results (from errors)
summaries_list = summaries_list[!sapply(summaries_list, is.null)]

if(length(summaries_list) == 0) {
    stop("❌ ERROR: No files were successfully processed!")
}

# Combine new results
cat("Combining new results...\n")
new_summaries = data.table::rbindlist(summaries_list)
cat("Processed", nrow(new_summaries), "new records\n")

# Combine with existing results
if(!is.null(existing_results)) {
    cat("Merging with existing results (", nrow(existing_results), " existing records)...\n")
    
    # Convert to data.table for faster operations
    if(!inherits(existing_results, "data.table")) {
        existing_dt = as.data.table(existing_results)
    } else {
        existing_dt = existing_results
    }
    
    if(!inherits(new_summaries, "data.table")) {
        new_dt = as.data.table(new_summaries)
    } else {
        new_dt = new_summaries
    }
    
    cat("  Combining data.tables...\n")
    flush.console()  # Ensure output is visible
    
    # Use data.table rbindlist for faster combination
    # Check if columns match, if not use fill=TRUE
    existing_cols = names(existing_dt)
    new_cols = names(new_dt)
    
    if(!setequal(existing_cols, new_cols)) {
        cat("  WARNING: Column mismatch detected. Using fill=TRUE\n")
        cat("  Existing cols:", length(existing_cols), "New cols:", length(new_cols), "\n")
        combined_dt = rbindlist(list(existing_dt, new_dt), fill = TRUE)
    } else {
        combined_dt = rbindlist(list(existing_dt, new_dt))
    }
    
    cat("  Combined records:", nrow(combined_dt), "\n")
    flush.console()
    
    cat("  Removing duplicates...\n")
    flush.console()
    
    # Use data.table unique for faster deduplication
    # Check which columns exist for unique check
    key_cols = intersect(c("samp_name", "metric", "db_name", "sampleID"), names(combined_dt))
    if(length(key_cols) > 0) {
        cat("  Using key columns for deduplication:", paste(key_cols, collapse=", "), "\n")
        flush.console()
        all_summaries = unique(combined_dt, by = key_cols)
    } else {
        cat("  WARNING: No key columns found, using all columns for deduplication\n")
        flush.console()
        all_summaries = unique(combined_dt)
    }
    
    cat("Total records after combining and deduplication:", nrow(all_summaries), "\n")
    flush.console()
    
    # Convert back to tibble for write_csv compatibility
    all_summaries = as_tibble(all_summaries)
} else {
    all_summaries = new_summaries
    if(!inherits(all_summaries, "tbl")) {
        all_summaries = as_tibble(all_summaries)
    }
}

# Save combined results
cat("Saving results to:", output_file, "\n")
cat("  Writing", nrow(all_summaries), "rows...\n")
write_csv(all_summaries, output_file)
cat("  ✓ File written successfully\n")

# Update processed files log
cat("Updating processed files log:", processed_log_file, "\n")
all_processed_files = c(processed_files, new_files)
writeLines(all_processed_files, processed_log_file)

cat("✓ Successfully processed", length(new_files), "files\n")
cat("✓ Total records in output:", nrow(all_summaries), "\n")
cat("✓ Results saved to:", output_file, "\n")
