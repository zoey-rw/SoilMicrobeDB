# Calculate % of fungal reads classified to each taxonomic rank
# Input: 
#   - Filtered kreport files after Bracken reassignment (02_bracken_output/*_filtered_kraken.kreport)
#     These files contain full KPCOFGS taxonomic information
# Output: Table showing % of fungal reads classified at each taxonomic rank
#         after Bracken reassignment (final processed data)
# Note: Uses filtered kreport files which have complete taxonomic hierarchy (KPCOFGS)
#       and represent the final processed classification results after Bracken reassignment
#       
#       Cumulative percentages: For each rank, calculates % of fungal reads classified
#       at that rank OR MORE SPECIFIC. Percentages should decrease with more specific ranks.
#       If domain and phylum are identical, it means there are no reads classified
#       specifically at domain level - all fungal reads are at phylum or below.

library(tidyverse)
library(data.table)
library(pavian)
source("scripts/helper_functions.r")

fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota",
                 "Chytridiomycota", "Cryptomycota", "Mucoromycota",
                 "Microsporidia", "Olpidiomycota", "Zoopagomycota")
fungal_pattern <- paste(fungal_phyla, collapse = "|")

# Rank mapping
rank_map <- c("domain" = "D", "kingdom" = "K", "phylum" = "P", "class" = "C", 
              "order" = "O", "family" = "F", "genus" = "G", "species" = "S", "strain" = "S1")

# Databases to process (exclude pluspf)
databases_to_check <- c("soil_microbe_db", "gtdb_207", "gtdb_207_unfiltered")

# Helper function to normalize sample_id and extract db_name
normalize_sample_ids <- function(sample_ids, db_default) {
    # Parse sample_ids
    parts_list <- strsplit(sample_ids, "COMP_", fixed = TRUE)
    
    # Extract components
    sampleID <- sapply(parts_list, function(x) str_remove(x[1], "-$"))
    db_name_from_sample <- sapply(parts_list, function(x) if(length(x) >= 2) x[2] else "")
    
    # Normalize sample_id by removing rank-specific suffixes
    sample_id_normalized <- str_remove(sample_ids, "_[a-z]+_filtered$|_filtered$")
    
    # Extract db_name by removing rank-specific suffixes
    db_name <- ifelse(db_name_from_sample == "", 
                     db_default,
                     str_remove(db_name_from_sample, "_[a-z]+_filtered$|_filtered$"))
    db_name <- ifelse(is.na(db_name) | db_name == "", db_default, db_name)
    
    return(data.frame(
        sampleID = sampleID,
        db_name = db_name,
        sample_id_normalized = sample_id_normalized,
        stringsAsFactors = FALSE
    ))
}

# Build taxonomy_id -> is_fungi mapping from merged lineage files (with caching)
cache_file <- "data/classification/analysis_files/taxonomy_is_fungi_cache.rds"
lineage_dir <- "data/classification/taxonomic_rank_summaries"

if(file.exists(cache_file)) {
    taxonomy_is_fungi <- readRDS(cache_file)
} else {
    taxonomy_is_fungi <- data.frame(taxonomy_id = integer(), is_fungi = logical(), stringsAsFactors = FALSE)
    for(db in databases_to_check) {
        # Check all ranks to get comprehensive mapping (species, genus, family, order, class, phylum)
        for(rank in c("species", "genus", "family", "order", "class", "phylum")) {
            filename <- paste0(db, "_", rank, "_merged_lineage.csv")
            lineage_file <- file.path(lineage_dir, filename)
            
            if(file.exists(lineage_file)) {
                tryCatch({
                    lineage_df <- read_csv(lineage_file, col_select = c("taxonomy_id", "lineage"), show_col_types = FALSE)
                    if(nrow(lineage_df) > 0) {
                        lineage_df <- lineage_df %>%
                            distinct(taxonomy_id, .keep_all = TRUE) %>%
                            mutate(is_fungi = grepl(fungal_pattern, lineage, ignore.case = TRUE)) %>%
                            select(taxonomy_id, is_fungi)
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
    
    if(nrow(taxonomy_is_fungi) > 0) {
        taxonomy_is_fungi <- taxonomy_is_fungi %>%
            group_by(taxonomy_id) %>%
            summarize(is_fungi = any(is_fungi), .groups = "drop") %>%
            arrange(taxonomy_id)
        setDT(taxonomy_is_fungi)
        
        dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
        saveRDS(taxonomy_is_fungi, cache_file)
    } else {
        taxonomy_is_fungi <- NULL
    }
}

# Cache for read_report3 results
kreport_cache <- new.env()
kreport_cache_dir <- "data/classification/analysis_files/kreport_cache"
dir.create(kreport_cache_dir, recursive = TRUE, showWarnings = FALSE)

# Function to read and cache kreport files
read_kreport_cached <- function(kreport_file) {
    # Create cache key from file path and modification time
    file_mtime <- file.mtime(kreport_file)
    cache_key <- paste0(basename(kreport_file), "_", as.numeric(file_mtime))
    cache_file <- file.path(kreport_cache_dir, paste0(cache_key, ".rds"))
    
    # Check if cached version exists
    if(file.exists(cache_file)) {
        return(readRDS(cache_file))
    }
    
    # Read the file
    first_line <- readLines(kreport_file, n = 1)
    has_header <- grepl("^[a-zA-Z]", first_line)
    
    # Read with all ranks needed (R, U for overall; D-K-P-C-O-F-G-S-S1 for both)
    if(has_header) {
        report <- read_report3(kreport_file, keep_taxRanks = c("R", "U", "D", "K", "P", "C", "O", "F", "G", "S", "S1"), has_header = TRUE)
    } else {
        report <- read_report3(kreport_file, keep_taxRanks = c("R", "U", "D", "K", "P", "C", "O", "F", "G", "S", "S1"), has_header = FALSE)
    }
    
    if(is.null(report) || nrow(report) == 0) {
        return(NULL)
    }
    
    if(!inherits(report, "data.table")) {
        setDT(report)
    }
    
    # Build taxLineage for fungi identification (needed for fungi calculations)
    if("depth" %in% names(report)) {
        original_name <- gsub("^[a-z]+_", "", report$name)
        report$taxLineage <- original_name
        if(nrow(report) > 1) {
            parent_stack <- list()
            for(i in seq_len(nrow(report))) {
                current_depth <- report$depth[i]
                while(length(parent_stack) > 0 && parent_stack[[length(parent_stack)]]$depth >= current_depth) {
                    parent_stack <- parent_stack[-length(parent_stack)]
                }
                if(length(parent_stack) > 0) {
                    parent_lineage <- parent_stack[[length(parent_stack)]]$lineage
                    report$taxLineage[i] <- paste(parent_lineage, original_name[i], sep = "|")
                }
                parent_stack[[length(parent_stack) + 1]] <- list(depth = current_depth, lineage = report$taxLineage[i])
            }
        }
    }
    
    # Cache the result
    saveRDS(report, cache_file)
    return(report)
}

# Function to parse sample ID from kreport filename
# Handles both patterns: _filtered_kraken.kreport and __filtered_kraken.kreport
parse_sample_id_from_kreport <- function(filename) {
    # Remove both single and double underscore patterns
    samp_name <- gsub("_filtered_kraken.kreport|__filtered_kraken.kreport", "", basename(filename))
    
    parts <- strsplit(samp_name, "COMP_", fixed = TRUE)[[1]]
    if(length(parts) >= 2) {
        sampleID_temp <- parts[1]
        db_name <- parts[2]
        # Remove trailing underscores and suffixes
        sampleID <- sub(paste0("_", db_name), "", samp_name, fixed = TRUE)
        sampleID <- gsub("_+$", "", sampleID)  # Remove trailing underscores
        db_name <- gsub("_filtered|_kraken|^_+|_+$", "", db_name)  # Clean up db_name
    } else {
        if(grepl("soil_microbe_db", samp_name)) {
            sampleID <- gsub("_soil_microbe_db.*", "", samp_name)
            db_name <- "soil_microbe_db"
        } else if(grepl("gtdb_207", samp_name)) {
            sampleID <- gsub("_gtdb_207.*", "", samp_name)
            db_name <- "gtdb_207"
        } else if(grepl("pluspf", samp_name)) {
            sampleID <- gsub("_pluspf.*", "", samp_name)
            db_name <- "pluspf"
        } else {
            sampleID <- samp_name
            db_name <- "unknown"
        }
    }
    # Clean up any remaining double underscores
    sampleID <- gsub("_+", "_", sampleID)
    sampleID <- gsub("^_|_$", "", sampleID)
    
    return(list(sampleID = sampleID, db_name = db_name, samp_name = samp_name))
}

# Function to read filtered kreport and extract fungal reads by rank
read_filtered_kreport_fungi_by_rank <- function(report, sample_info, taxonomy_is_fungi_map) {
    if(is.null(report) || nrow(report) == 0) {
        return(NULL)
    }
    
    # Check if taxLineage exists (should be built by read_kreport_cached)
    if(!"taxLineage" %in% names(report)) {
        warning("taxLineage not found in report. Skipping.")
        return(NULL)
    }
    
    # Identify fungi: use ONLY taxLineage (consistent method for all ranks)
    # taxLineage is built from depth information and includes full hierarchy
    # This ensures consistent identification across all ranks (phylum, class, order, family, genus, species)
    
    # Check if entry is in Fungi kingdom lineage using taxLineage
    # This works for all ranks since taxLineage includes full hierarchy
    report[, is_fungi := grepl("(^|\\|)Fungi(\\||$)", taxLineage, ignore.case = TRUE)]
    
    # Calculate total fungal reads: sum all taxonReads (reads specifically assigned at any rank)
    # This gives us the total number of fungal reads across all ranks
    # Using taxonReads avoids double-counting that would occur with cladeReads
    rank_hierarchy_letters <- c("D", "K", "P", "C", "O", "F", "G", "S", "S1")
    total_fungal_reads <- 0
    
    for(rank_letter in rank_hierarchy_letters) {
        fungi_at_rank <- report[taxRank == rank_letter & is_fungi == TRUE]
        if(nrow(fungi_at_rank) > 0) {
            total_fungal_reads <- total_fungal_reads + sum(fungi_at_rank$taxonReads, na.rm = TRUE)
        }
    }
    
    # Require taxonReads - no fallback
    if(total_fungal_reads == 0) {
        return(NULL)
    }
    
    # Calculate cumulative reads at each rank (reads at this rank OR MORE SPECIFIC)
    # For cumulative percentages, we need to sum taxonReads from this rank and all more specific ranks
    # taxonReads = reads assigned specifically at this rank (not including children)
    # cladeReads = reads at this rank AND all descendants (cumulative)
    # 
    # For cumulative: sum taxonReads from target rank down to most specific
    # This avoids double-counting that would occur from summing cladeReads
    rank_hierarchy <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
    rank_hierarchy_letters <- c("D", "K", "P", "C", "O", "F", "G", "S", "S1")
    
    results_list <- list()
    for(rank_name in names(rank_map)) {
        rank_letter <- rank_map[rank_name]
        
        # Get all ranks at this level or more specific
        rank_idx <- which(rank_hierarchy == rank_name)
        if(length(rank_idx) == 0) next
        
        ranks_to_sum <- rank_hierarchy[rank_idx:length(rank_hierarchy)]
        ranks_to_sum_letters <- rank_hierarchy_letters[rank_idx:length(rank_hierarchy_letters)]
        
        # Sum taxonReads (reads specifically at each rank) across all ranks from target down
        # This gives cumulative reads at target rank or more specific
        fungi_reads_at_rank <- 0
        for(r_letter in ranks_to_sum_letters) {
            fungi_at_r <- report[taxRank == r_letter & is_fungi == TRUE]
            if(nrow(fungi_at_r) > 0) {
                fungi_reads_at_rank <- fungi_reads_at_rank + sum(fungi_at_r$taxonReads, na.rm = TRUE)
            }
        }
        
        if(fungi_reads_at_rank > 0) {
            pct_at_rank <- (fungi_reads_at_rank / total_fungal_reads) * 100
            results_list[[rank_name]] <- data.frame(
                sampleID = sample_info$sampleID,
                db_name = sample_info$db_name,
                samp_name = sample_info$samp_name,
                taxonomic_rank = rank_name,
                filter_status = "after_bracken_reassignment",
                fungal_reads = fungi_reads_at_rank,
                total_fungal_reads = total_fungal_reads,
                pct_fungal_reads = pct_at_rank,
                stringsAsFactors = FALSE
            )
        }
    }
    
    if(length(results_list) == 0) {
        return(NULL)
    }
    
    return(bind_rows(results_list))
}

# Load sequencing depth data (original total read counts, same across databases)
seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds")

# Function to get original total reads from seq_depth_df.rds
get_original_total_reads <- function(sampleID, db_name) {
    if(is.null(sampleID) || is.null(db_name)) {
        return(NULL)
    }
    
    match_row <- seq_depth_df[seq_depth_df$sampleID == sampleID & seq_depth_df$db_name == db_name, ]
    if(nrow(match_row) > 0 && !is.na(match_row$seq_depth[1]) && match_row$seq_depth[1] > 0) {
        return(match_row$seq_depth[1])
    }
    
    return(NULL)
}

# Function to read filtered kreport and extract overall classification percentages by rank
read_filtered_kreport_overall_by_rank <- function(report, sample_info, kreport_file) {
    if(is.null(report) || nrow(report) == 0) {
        return(NULL)
    }
    
    # Get original total reads from seq_depth_df.rds (original sequencing depth, same across databases)
    original_total <- get_original_total_reads(sample_info$sampleID, sample_info$db_name)
    
    if(!is.null(original_total) && original_total > 0) {
        total_reads <- original_total
    } else {
        # Fallback: use root's cladeReads from Bracken file
        # This should rarely happen if seq_depth_df.rds is complete
        root_entry <- report[taxRank == "R" | name == "root" | name == "r_root"]
        if(nrow(root_entry) > 0) {
            total_reads <- root_entry$cladeReads[1]
        } else {
            # Last resort: sum classified reads
            unclassified_entry <- report[taxRank == "U" | grepl("unclassified", name, ignore.case = TRUE)]
            rank_hierarchy_letters <- c("D", "K", "P", "C", "O", "F", "G", "S", "S1")
            total_reads <- 0
            for(rank_letter in rank_hierarchy_letters) {
                reads_at_rank <- report[taxRank == rank_letter]
                if(nrow(reads_at_rank) > 0) {
                    total_reads <- total_reads + sum(reads_at_rank$taxonReads, na.rm = TRUE)
                }
            }
            if(nrow(unclassified_entry) > 0) {
                total_reads <- total_reads + sum(unclassified_entry$cladeReads, na.rm = TRUE)
            }
        }
    }
    
    if(total_reads == 0) {
        return(NULL)
    }
    
    # Calculate cumulative reads at each rank (reads at this rank OR MORE SPECIFIC)
    # Use cladeReads which already includes all descendants (cumulative)
    # cladeReads = reads at this rank AND all more specific ranks
    rank_hierarchy <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
    rank_hierarchy_letters <- c("D", "K", "P", "C", "O", "F", "G", "S", "S1")
    
    results_list <- list()
    for(rank_name in names(rank_map)) {
        rank_letter <- rank_map[rank_name]
        
        # Get all entries at this rank and sum their cladeReads
        # cladeReads already includes all descendants, so this gives cumulative reads
        reads_at_rank <- report[taxRank == rank_letter, sum(cladeReads, na.rm = TRUE)]
        
        if(reads_at_rank > 0) {
            pct_at_rank <- (reads_at_rank / total_reads) * 100
            results_list[[rank_name]] <- data.frame(
                sampleID = sample_info$sampleID,
                db_name = sample_info$db_name,
                samp_name = sample_info$samp_name,
                taxonomic_rank = rank_name,
                filter_status = "after_bracken_reassignment",
                classified_reads = reads_at_rank,
                total_reads = total_reads,
                pct_classified = pct_at_rank,
                stringsAsFactors = FALSE
            )
        }
    }
    
    if(length(results_list) == 0) {
        return(NULL)
    }
    
    return(bind_rows(results_list))
}

# Process Bracken reassignment from filtered kreport files
bracken_results <- list()
overall_results <- list()

# Find filtered kreport files (check multiple locations)
kreport_dirs <- c(
    "data/classification/02_bracken_output",
    "data/NEON_metagenome_classification/02_bracken_output",
    "data/classification/analysis_files"
)

filtered_kreport_files <- character()
for(kreport_dir in kreport_dirs) {
    if(dir.exists(kreport_dir)) {
        # Look for main filtered kreport files (not rank-specific ones)
        # Pattern 1: _filtered_kraken.kreport (original pattern, single underscore)
        # Pattern 2: __filtered_kraken.kreport (new pattern with double underscore)
        # Exclude rank-specific files like *_bracken_phylums.kreport, *_bracken_species.kreport, etc.
        # IMPORTANT: Check double underscore pattern FIRST to avoid matching it with single underscore pattern
        files2 <- list.files(kreport_dir, pattern = "__filtered_kraken\\.kreport$", full.names = TRUE, recursive = TRUE)
        files1 <- list.files(kreport_dir, pattern = "_filtered_kraken\\.kreport$", full.names = TRUE, recursive = TRUE)
        
        # Remove double underscore files from single underscore results to avoid duplicates
        files1 <- files1[!files1 %in% files2]
        
        # Exclude rank-specific bracken files
        all_files <- c(files1, files2)
        # Remove files that have rank-specific suffixes (bracken_phylums, bracken_species, etc.)
        all_files <- all_files[!grepl("_bracken_(phylums|classes|orders|families|genuses|species|domains)\\.kreport$", all_files)]
        
        filtered_kreport_files <- c(filtered_kreport_files, all_files)
    }
}

if(length(filtered_kreport_files) == 0) {
    stop("No filtered kreport files found in ", kreport_dir)
}

# Process each filtered kreport file
# Read each file once and use for both fungi and overall calculations
for(i in seq_along(filtered_kreport_files)) {
    kreport_file <- filtered_kreport_files[i]
    
    # Read and cache the report (includes taxLineage building)
    report <- read_kreport_cached(kreport_file)
    if(is.null(report)) {
        next
    }
    
    # Get sample info once
    sample_info <- parse_sample_id_from_kreport(kreport_file)
    
    # Process for fungi
    result_fungi <- read_filtered_kreport_fungi_by_rank(report, sample_info, taxonomy_is_fungi)
    if(!is.null(result_fungi)) {
        bracken_results[[length(bracken_results) + 1]] <- result_fungi
    }
    
    # Process for overall (using same report)
    result_overall <- read_filtered_kreport_overall_by_rank(report, sample_info, kreport_file)
    if(!is.null(result_overall)) {
        overall_results[[length(overall_results) + 1]] <- result_overall
    }
}

bracken_df <- bind_rows(bracken_results)

# Aggregate Bracken results by sampleID + db_name + filter_status + taxonomic_rank
# (in case there are duplicates - shouldn't happen with kreport files, but just in case)
if(nrow(bracken_df) > 0) {
    bracken_df <- bracken_df %>%
        group_by(sampleID, db_name, filter_status, taxonomic_rank) %>%
        summarize(
            samp_name = first(samp_name),
            fungal_reads = sum(fungal_reads, na.rm = TRUE),
            total_fungal_reads = first(total_fungal_reads),
            pct_fungal_reads = ifelse(first(total_fungal_reads) > 0,
                                    (sum(fungal_reads, na.rm = TRUE) / first(total_fungal_reads)) * 100, 0),
            .groups = "drop"
        )
}

# Use Bracken results (final processed data)
all_results <- bracken_df

if(nrow(all_results) == 0) {
    stop("No results generated! Check that input files exist and contain fungal reads.")
}

# Create output directory
output_dir <- "data/classification/analysis_files"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Write output table
output_file <- file.path(output_dir, "fungal_reads_pct_by_rank.csv")
write_csv(all_results, output_file)

# Process overall classification results
overall_df <- bind_rows(overall_results)

# Aggregate overall results by sampleID + db_name + filter_status + taxonomic_rank
if(nrow(overall_df) > 0) {
    overall_df <- overall_df %>%
        group_by(sampleID, db_name, filter_status, taxonomic_rank) %>%
        summarize(
            samp_name = first(samp_name),
            classified_reads = sum(classified_reads, na.rm = TRUE),
            total_reads = first(total_reads),
            pct_classified = ifelse(first(total_reads) > 0,
                                  (sum(classified_reads, na.rm = TRUE) / first(total_reads)) * 100, 0),
            .groups = "drop"
        )
}

# Write overall classification table
if(nrow(overall_df) > 0) {
    overall_output_file <- file.path(output_dir, "overall_classification_pct_by_rank.csv")
    write_csv(overall_df, overall_output_file)
}

# Generate summary table with all rank comparisons
all_results_filtered <- all_results %>% filter(db_name != "pluspf")

# Generate all rank comparisons: each rank compared to all lower ranks
rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
rank_comparisons <- list()

for(i in seq_along(rank_hierarchy)) {
    higher_rank <- rank_hierarchy[i]
    if(i < length(rank_hierarchy)) {
        for(j in (i+1):length(rank_hierarchy)) {
            lower_rank <- rank_hierarchy[j]
            comp_name <- paste0(stringr::str_to_title(higher_rank), " not ", lower_rank)
            rank_comparisons[[comp_name]] <- c(higher_rank, lower_rank)
        }
    }
}

summary_table <- data.frame()

for(comp_name in names(rank_comparisons)) {
    ranks <- rank_comparisons[[comp_name]]
    comparison <- all_results_filtered %>%
        filter(taxonomic_rank %in% ranks) %>%
        select(sampleID, db_name, taxonomic_rank, pct_fungal_reads) %>%
        pivot_wider(names_from = taxonomic_rank, values_from = pct_fungal_reads)
    
    # Check if both rank columns exist
    if(ranks[1] %in% names(comparison) && ranks[2] %in% names(comparison)) {
        comparison <- comparison %>%
            mutate(diff = .data[[ranks[1]]] - .data[[ranks[2]]])
        
        summary_table <- rbind(summary_table, data.frame(
            Comparison = comp_name,
            Mean_pct = round(mean(comparison$diff, na.rm = TRUE), 2),
            Median_pct = round(median(comparison$diff, na.rm = TRUE), 2),
            SD_pct = round(sd(comparison$diff, na.rm = TRUE), 2)
        ))
    }
}

# Order by rank hierarchy (kingdom first, then phylum, etc.)
rank_order <- c()
for(i in seq_along(rank_hierarchy)) {
    higher_rank <- rank_hierarchy[i]
    if(i < length(rank_hierarchy)) {
        for(j in (i+1):length(rank_hierarchy)) {
            lower_rank <- rank_hierarchy[j]
            comp_name <- paste0(stringr::str_to_title(higher_rank), " not ", lower_rank)
            rank_order <- c(rank_order, comp_name)
        }
    }
}

summary_table <- summary_table %>%
    mutate(Order = match(Comparison, rank_order)) %>%
    arrange(Order) %>%
    select(-Order)

# Write summary table
summary_file <- file.path(output_dir, "fungal_reads_pct_by_rank_summary.csv")
write_csv(summary_table, summary_file)

# Generate summary table for overall classification percentages
if(nrow(overall_df) > 0) {
    overall_df_filtered <- overall_df %>% filter(db_name != "pluspf")
    
    # Calculate summary statistics for each rank
    rank_hierarchy_full <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
    
    overall_summary_table <- overall_df_filtered %>%
        filter(taxonomic_rank %in% rank_hierarchy_full) %>%
        group_by(taxonomic_rank) %>%
        summarize(
            Mean_pct = round(mean(pct_classified, na.rm = TRUE), 2),
            Median_pct = round(median(pct_classified, na.rm = TRUE), 2),
            SD_pct = round(sd(pct_classified, na.rm = TRUE), 2),
            .groups = "drop"
        ) %>%
        mutate(Order = match(taxonomic_rank, rank_hierarchy_full)) %>%
        arrange(Order) %>%
        select(Taxonomic_rank = taxonomic_rank, Mean_pct, Median_pct, SD_pct)
    
    # Write overall summary table
    overall_summary_file <- file.path(output_dir, "overall_classification_pct_by_rank_summary.csv")
    write_csv(overall_summary_table, overall_summary_file)
}
