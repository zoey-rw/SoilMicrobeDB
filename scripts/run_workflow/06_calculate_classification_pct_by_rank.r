# Calculate classification percentage at each taxonomic rank
# Shows what % of total sequencing depth was classified to each specific rank
# Calculates for BOTH before and after Architeuthis filtering
# Output: Table for supplement showing classification % at phylum → class → order → family → genus → species → strain
# Note: Calculates proportion classified EXACTLY to each rank (not cumulative)
# Denominator: total sequencing depth (seq_depth)
# Uses kreport files from hard drive

library(tidyverse)
library(data.table)
library(pavian)
source("scripts/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")

# Read sequencing depth for normalization
seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds")

# Define fungal phyla for filtering
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")

# Set up paths - check both local and HARDDRIVE locations
# Before filtering: original Kraken2 kreport files
kreport_before_local = "data/classification/01_kraken_output"
kreport_before_harddrive = "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/01_kraken_output"

# After filtering: filtered kreport files
kreport_after_local = "data/classification/02_bracken_output"
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

# Filter to soil_microbe_db only
kreport_before_files <- kreport_before_files[grepl("soil_microbe_db", kreport_before_files)]
kreport_after_files <- kreport_after_files[grepl("soil_microbe_db", kreport_after_files)]

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

# Function to calculate classification % at a specific rank from kreport
# kreport files contain all taxonomic ranks, so we can extract reads at each rank
calculate_pct_from_kreport <- function(kreport_file, rank_code, rank_name, filter_fungi = FALSE, filter_status = "unknown") {
    # Read kreport file - handle both 6-column and 8-column formats
    # Try reading with explicit column detection
    first_line <- readLines(kreport_file, n = 1)
    has_header <- grepl("^[a-zA-Z]", first_line)
    
    if(has_header) {
        report <- read_report3(kreport_file, keep_taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S", "S1"), has_header = TRUE)
    } else {
        # Check number of columns
        test_line <- read.table(kreport_file, sep = "\t", nrows = 1, header = FALSE)
        n_cols <- ncol(test_line)
        
        if(n_cols == 8) {
            # 8-column format: percentage, cladeReads, taxonReads, nKmers, n_unique_kmers, taxRank, taxID, name
            report <- read.table(kreport_file, sep = "\t", header = FALSE,
                               col.names = c("percentage", "cladeReads", "taxonReads", "nKmers", "n_unique_kmers", 
                                            "taxRank", "taxID", "name"),
                               quote = "", stringsAsFactors = FALSE, comment.char = "#")
        } else if(n_cols == 6) {
            # 6-column format: percentage, cladeReads, taxonReads, taxRank, taxID, name
            report <- read.table(kreport_file, sep = "\t", header = FALSE,
                               col.names = c("percentage", "cladeReads", "taxonReads", "taxRank", "taxID", "name"),
                               quote = "", stringsAsFactors = FALSE, comment.char = "#")
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
    
    # Filter to specific rank
    rank_map <- c("phylum" = "P", "class" = "C", "order" = "O", 
                  "family" = "F", "genus" = "G", "species" = "S", "strain" = "S1")
    
    if(!rank_name %in% names(rank_map)) {
        stop("Unknown rank: ", rank_name)
    }
    
    rank_letter <- rank_map[rank_name]
    
    # Get reads classified exactly at this rank (not including children)
    # In kreport format, taxonReads is reads at this specific taxon
    rank_data <- report %>%
        filter(taxRank == rank_letter)
    
    # Filter fungi if requested
    if(filter_fungi) {
        # Check if lineage contains fungal phyla
        rank_data <- rank_data %>%
            filter(grepl(paste(fungal_phyla, collapse = "|"), name, ignore.case = TRUE))
    }
    
    # Sum reads at this rank
    reads_classified <- sum(rank_data$taxonReads, na.rm = TRUE)
    
    # Parse sample ID
    is_filtered <- grepl("_filtered_kraken.kreport", kreport_file)
    sample_info <- parse_sample_id(kreport_file, is_filtered = is_filtered)
    
    return(data.frame(
        sampleID = sample_info$sampleID,
        db_name = sample_info$db_name,
        samp_name = sample_info$samp_name,
        taxonomic_rank = rank_name,
        filter_status = filter_status,
        reads_classified = reads_classified,
        stringsAsFactors = FALSE
    ))
}

# Process all kreport files for each rank
cat("Processing kreport files to extract classification by rank...\n")

all_domain_results_list <- list()
fungi_results_list <- list()

ranks <- c("phylum" = "P", "class" = "C", "order" = "O", 
           "family" = "F", "genus" = "G", "species" = "S", "strain" = "S1")

# Process BEFORE filtering files
if(length(kreport_before_files) > 0) {
    cat("Processing BEFORE filtering files...\n")
    for(rank_name in names(ranks)) {
        cat("  Processing", rank_name, "(before filtering)...\n")
        rank_code <- ranks[rank_name]
        
        # Process all-domain
        all_domain_rank_results <- lapply(kreport_before_files, function(f) {
            tryCatch({
                calculate_pct_from_kreport(f, rank_name, rank_name, filter_fungi = FALSE, filter_status = "before")
            }, error = function(e) {
                cat("    Error processing", basename(f), ":", e$message, "\n")
                return(NULL)
            })
        })
        all_domain_rank_results <- bind_rows(Filter(Negate(is.null), all_domain_rank_results))
        if(nrow(all_domain_rank_results) > 0) {
            all_domain_results_list[[paste0(rank_name, "_before")]] <- all_domain_rank_results
        }
        
        # Process fungi-specific
        fungi_rank_results <- lapply(kreport_before_files, function(f) {
            tryCatch({
                calculate_pct_from_kreport(f, rank_name, rank_name, filter_fungi = TRUE, filter_status = "before")
            }, error = function(e) {
                cat("    Error processing", basename(f), ":", e$message, "\n")
                return(NULL)
            })
        })
        fungi_rank_results <- bind_rows(Filter(Negate(is.null), fungi_rank_results))
        if(nrow(fungi_rank_results) > 0) {
            fungi_results_list[[paste0(rank_name, "_before")]] <- fungi_rank_results
        }
    }
}

# Process AFTER filtering files
if(length(kreport_after_files) > 0) {
    cat("Processing AFTER filtering files...\n")
    for(rank_name in names(ranks)) {
        cat("  Processing", rank_name, "(after filtering)...\n")
        rank_code <- ranks[rank_name]
        
        # Process all-domain
        all_domain_rank_results <- lapply(kreport_after_files, function(f) {
            tryCatch({
                calculate_pct_from_kreport(f, rank_name, rank_name, filter_fungi = FALSE, filter_status = "after")
            }, error = function(e) {
                cat("    Error processing", basename(f), ":", e$message, "\n")
                return(NULL)
            })
        })
        all_domain_rank_results <- bind_rows(Filter(Negate(is.null), all_domain_rank_results))
        if(nrow(all_domain_rank_results) > 0) {
            all_domain_results_list[[paste0(rank_name, "_after")]] <- all_domain_rank_results
        }
        
        # Process fungi-specific
        fungi_rank_results <- lapply(kreport_after_files, function(f) {
            tryCatch({
                calculate_pct_from_kreport(f, rank_name, rank_name, filter_fungi = TRUE, filter_status = "after")
            }, error = function(e) {
                cat("    Error processing", basename(f), ":", e$message, "\n")
                return(NULL)
            })
        })
        fungi_rank_results <- bind_rows(Filter(Negate(is.null), fungi_rank_results))
        if(nrow(fungi_rank_results) > 0) {
            fungi_results_list[[paste0(rank_name, "_after")]] <- fungi_rank_results
        }
    }
}

# Combine all results
all_domain_results <- bind_rows(all_domain_results_list)
fungi_results <- bind_rows(fungi_results_list)

# Merge with sequencing depth to calculate percentage
if(nrow(all_domain_results) > 0) {
    all_domain_results <- all_domain_results %>%
        left_join(seq_depth_df %>% select(sampleID, db_name, seq_depth), 
                 by = c("sampleID", "db_name")) %>%
        mutate(
            reads_classified = replace_na(reads_classified, 0),
            pct_classified = (reads_classified / seq_depth) * 100
        )
} else {
    all_domain_results <- seq_depth_df %>%
        select(sampleID, db_name, seq_depth) %>%
        mutate(
            taxonomic_rank = NA_character_,
            filter_status = NA_character_,
            reads_classified = 0,
            pct_classified = 0
        )
}

if(nrow(fungi_results) > 0) {
    fungi_results <- fungi_results %>%
        left_join(seq_depth_df %>% select(sampleID, db_name, seq_depth), 
                 by = c("sampleID", "db_name")) %>%
        mutate(
            reads_classified = replace_na(reads_classified, 0),
            pct_classified = (reads_classified / seq_depth) * 100
        )
} else {
    fungi_results <- seq_depth_df %>%
        select(sampleID, db_name, seq_depth) %>%
        mutate(
            taxonomic_rank = NA_character_,
            filter_status = NA_character_,
            reads_classified = 0,
            pct_classified = 0
        )
}

# ============================================================================
# Aggregate results
# ============================================================================

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
cat("   Values represent reads classified EXACTLY to each rank\n")
cat("   filter_status indicates 'before' or 'after' Architeuthis filtering\n")

# Also save per-sample results
all_results <- bind_rows(
    all_domain_results %>% mutate(taxon_group = "all_domain"),
    fungi_results %>% mutate(taxon_group = "fungi")
)

write_csv(all_results, "data/classification/analysis_files/classification_pct_by_rank_per_sample.csv")
cat("✅ Saved per-sample classification percentages by rank to: data/classification/analysis_files/classification_pct_by_rank_per_sample.csv\n")
