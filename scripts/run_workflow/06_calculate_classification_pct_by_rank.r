# Calculate classification percentage at each taxonomic rank
# Shows what % of total sequencing depth was classified to each specific rank
# Calculates for BOTH before and after Architeuthis filtering
# Output: Table for supplement showing classification % at phylum → class → order → family → genus → species → strain
# Note: Calculates CUMULATIVE proportion (reads classified at this rank OR HIGHER/more specific)
# Denominator: total sequencing depth (seq_depth)
# Uses kreport files from hard drive

library(tidyverse)
library(data.table)
library(future.apply)
library(pavian)
source("scripts/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")

# Set up parallel processing
n_workers <- min(18, availableCores() - 1)
plan(multisession, workers = n_workers)
cat("Using", n_workers, "parallel workers\n")

# Read sequencing depth for normalization
seq_depth_file <- "data/NEON_metagenome_classification/seq_depth_df.rds"
if(!file.exists(seq_depth_file)) {
    seq_depth_file <- "data/classification/analysis_files/seq_depth_df.rds"
}
if(!file.exists(seq_depth_file)) {
    stop("❌ MISSING FILE: seq_depth_df.rds not found!\n",
         "   Checked: data/NEON_metagenome_classification/seq_depth_df.rds\n",
         "   Checked: data/classification/analysis_files/seq_depth_df.rds\n",
         "   Please run calculate_sequencing_depth.r first.")
}
seq_depth_df <- readRDS(seq_depth_file)

# Define fungal phyla for filtering
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")

# Set up paths - check both local and HARDDRIVE locations
# Before filtering: original Kraken2 kreport files
kreport_before_local = "data/NEON_metagenome_classification/01_kraken_output"
kreport_before_harddrive = "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/01_kraken_output"

# After filtering: filtered kreport files
kreport_after_local = "data/NEON_metagenome_classification/02_bracken_output"
kreport_after_harddrive = "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/02_bracken_output"

# Find kreport files BEFORE filtering (original Kraken2 output)
kreport_before_files <- list.files(kreport_before_local, pattern = "_kraken.kreport", full.names = TRUE, recursive = TRUE)
if(dir.exists(kreport_before_harddrive)) {
    kreport_before_harddrive_files <- list.files(kreport_before_harddrive, pattern = "_kraken.kreport", full.names = TRUE, recursive = TRUE)
    kreport_before_files <- c(kreport_before_files, kreport_before_harddrive_files)
    cat("✓ Found", length(kreport_before_harddrive_files), "kreport files (before filtering) on HARDDRIVE\n")
}

# Find kreport files AFTER filtering (Architeuthis filtered)
kreport_after_files <- list.files(kreport_after_local, pattern = "_filtered_kraken.kreport", full.names = TRUE, recursive = TRUE)
if(dir.exists(kreport_after_harddrive)) {
    kreport_after_harddrive_files <- list.files(kreport_after_harddrive, pattern = "_filtered_kraken.kreport", full.names = TRUE, recursive = TRUE)
    kreport_after_files <- c(kreport_after_files, kreport_after_harddrive_files)
    cat("✓ Found", length(kreport_after_harddrive_files), "kreport files (after filtering) on HARDDRIVE\n")
}

cat("Found", length(kreport_before_files), "kreport files BEFORE filtering\n")
cat("Found", length(kreport_after_files), "kreport files AFTER filtering\n")

if(length(kreport_before_files) == 0 && length(kreport_after_files) == 0) {
    stop("No kreport files found! Checked:\n",
         "  Before (local): ", kreport_before_local, "\n",
         "  Before (HARDDRIVE): ", kreport_before_harddrive, "\n",
         "  After (local): ", kreport_after_local, "\n",
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
read_kreport_all_ranks <- function(kreport_file, filter_status, fungal_phyla) {
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
        report$name <- gsub("^ *", "", report$name)
        report$name <- paste(tolower(report$taxRank), report$name, sep = "_")
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
    setDT(report)
    
    # Identify fungi once
    report[, is_fungi := grepl(paste(fungal_phyla, collapse = "|"), name, ignore.case = TRUE)]
    
    # Process all ranks
    results_list <- list()
    
    for(rank_name in names(rank_map)) {
        rank_letter <- rank_map[rank_name]
        
        # All-domain: sum cladeReads for this rank
        all_domain_reads <- report[taxRank == rank_letter, sum(cladeReads, na.rm = TRUE)]
        
        # Fungi-specific: sum cladeReads for this rank where is_fungi == TRUE
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

# Set up output files
output_file <- "data/NEON_metagenome_classification/summary_files/classification_pct_by_rank_per_sample.csv"
processed_log_file <- "data/NEON_metagenome_classification/summary_files/classification_pct_by_rank_processed_files.txt"

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
    
    # Get list of files already in results (by sampleID and filter_status)
    # Use sampleID + filter_status as unique identifier since samp_name format may vary
    existing_combos <- existing_results %>%
        distinct(sampleID, filter_status) %>%
        mutate(file_key = paste0(sampleID, "_", filter_status))
    
    # Create file keys for new files by parsing sampleID from filename
    file_keys_new <- sapply(seq_along(new_files), function(i) {
        kreport_file <- new_files[i]
        filter_status <- names(new_files)[i]
        sample_info <- parse_sample_id(kreport_file, is_filtered = (filter_status == "after"))
        paste0(sample_info$sampleID, "_", filter_status)
    })
    
    # Filter out files that are already in results
    new_files <- new_files[!file_keys_new %in% existing_combos$file_key]
    cat("After checking existing results:", length(new_files), "new files to process\n")
}

if(length(new_files) == 0) {
    cat("✓ All files have already been processed. Exiting.\n")
    quit(status = 0)
}

# Process new files in parallel
cat("Processing", length(new_files), "new kreport files in parallel...\n")

# Adjust number of workers based on files to process
n_workers <- min(n_workers, length(new_files))
plan(multisession, workers = n_workers)

# Process files in parallel
results <- future_lapply(seq_along(new_files), function(i) {
    kreport_file <- new_files[i]
    filter_status <- names(new_files)[i]
    
    tryCatch({
        read_kreport_all_ranks(kreport_file, filter_status, fungal_phyla)
    }, error = function(e) {
        cat("    Error processing", basename(kreport_file), ":", e$message, "\n")
        return(NULL)
    })
}, future.seed = TRUE)

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
if(nrow(all_results) > 0) {
    all_results <- all_results %>%
        left_join(seq_depth_df %>% select(sampleID, db_name, seq_depth), 
                 by = c("sampleID", "db_name")) %>%
        mutate(
            reads_classified = replace_na(reads_classified, 0),
            pct_classified = (reads_classified / seq_depth) * 100
        )
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

write_csv(classification_by_rank, "data/classification/analysis_files/classification_pct_by_rank.csv")
cat("✅ Saved classification percentages by rank to: data/classification/analysis_files/classification_pct_by_rank.csv\n")
cat("   Note: Percentages are % of total sequencing depth (seq_depth)\n")
cat("   Values are CUMULATIVE (reads classified at this rank OR HIGHER/more specific)\n")
cat("   filter_status indicates 'before' or 'after' Architeuthis filtering\n")

# Also save per-sample results
output_dir <- "data/NEON_metagenome_classification/summary_files"
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
