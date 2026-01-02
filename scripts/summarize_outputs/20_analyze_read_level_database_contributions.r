#!/usr/bin/env Rscript
# 20: Analyze read-level database contributions
# Classify reads by source using read-level scoring files
# Goal: Match each read's taxid with genome table to get source database
#       For taxids unique to one source: direct attribution
#       For taxids in multiple sources: proportional attribution based on genome counts
#
# Usage: Rscript scripts/summarize_outputs/20_analyze_read_level_database_contributions.r
#
# Input:  *_scores.output files from 02_bracken_output
#         Genome table with source information
# Output: Summary statistics showing read counts by source for each kingdom
#         Efficiency metrics, rank distribution, unique contributions,
#         overlap analysis, and cross-kingdom contributions
#
# Note: Scoring files only contain taxid-level information, not genome-level.
#       For taxids present in multiple sources, reads are attributed proportionally
#       based on the number of genomes from each source in the database.

library(data.table)

# File paths
# Check multiple possible locations (local and remote mount)
possible_scores_dirs <- c(
    "data/classification/02_bracken_output",
    "/Users/zoeywerbin/remote_soil_microbe_db/data/NEON_metagenome_classification/02_bracken_output",
    "/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/02_bracken_output",
    "data/NEON_metagenome_classification/02_bracken_output"
)

scores_dir <- NULL
for(dir in possible_scores_dirs) {
    if(dir.exists(dir)) {
        score_files_check <- list.files(dir, pattern = "_+scores.output$", full.names = TRUE, recursive = TRUE)
        # Filter to only soil_microbe_db outputs (exclude gtdb_207, pluspf, etc.)
        score_files_check <- score_files_check[grepl("_soil_microbe_db", basename(score_files_check))]
        if(length(score_files_check) > 0) {
            scores_dir <- dir
            cat("Using scores directory:", scores_dir, "\n")
            cat("Found", length(score_files_check), "soil_microbe_db score files\n")
            break
        }
    }
}

if(is.null(scores_dir)) {
    stop("No soil_microbe_db score files found in any of the checked directories:\n", 
         paste(possible_scores_dirs, collapse = "\n"))
}

genome_table_file <- "data/genome_database/soil_microbe_db_genome_table.csv"
output_dir <- "data/classification/analysis_files"

# Fungal phyla for identification
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")
fungal_pattern <- paste0("(", paste(fungal_phyla, collapse = "|"), ")")

# Read genome table
genome_table <- fread(genome_table_file, nThread = 8)

# Pre-calculate proportional weights for every taxid
# This creates a static lookup table: "For taxid X, Y% goes to Source A, Z% goes to Source B"
weights_dt <- genome_table[
    !is.na(ncbi_species_taxid),
    .(n_genomes = .N, kingdom = kingdom[1]),
    by = .(ncbi_species_taxid, source)
]

# Calculate total genomes per taxid and proportion
weights_dt[, total_genomes := sum(n_genomes), by = ncbi_species_taxid]
weights_dt[, proportion := n_genomes / total_genomes]

# Create kingdom lookup (one row per taxid)
taxid_to_kingdom <- weights_dt[, .(kingdom = kingdom[1]), by = ncbi_species_taxid]
setkey(taxid_to_kingdom, ncbi_species_taxid)

# Set key for efficient joins
setkey(weights_dt, ncbi_species_taxid)

# Find and process scoring files (recursively search subdirectories)
# Note: Values in output tables are AVERAGED across samples (not summed)
# Filter to only soil_microbe_db outputs (exclude other databases)
score_files <- list.files(scores_dir, pattern = "_+scores.output$", full.names = TRUE, recursive = TRUE)
score_files <- score_files[grepl("_soil_microbe_db", basename(score_files))]
cat("Processing", length(score_files), "soil_microbe_db score files from", scores_dir, "\n")

if(length(score_files) == 0) {
    stop("No scoring files found in ", scores_dir)
}

# Check if outputs are up to date
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_files <- c(
    "read_level_database_efficiency_metrics.csv",
    "read_level_database_unique_contributions.csv",
    "read_level_database_overlap_matrix.csv",
    "read_level_database_rank_distribution.csv",
    "read_level_database_cross_kingdom.csv",
    "read_level_all_kingdoms_read_counts_by_source.csv"
)

all_outputs_exist <- all(file.exists(file.path(output_dir, output_files)))

if(all_outputs_exist) {
    output_mtime <- max(file.mtime(file.path(output_dir, output_files)), na.rm = TRUE)
    input_mtime <- max(file.mtime(score_files), na.rm = TRUE)
    
    if(input_mtime <= output_mtime) {
        quit(status = 0)
    }
}

# ============================================================================
# MAP-REDUCE APPROACH: Process files individually and aggregate summaries
# ============================================================================

# Initialize lists to store summary results (not raw data)
summary_list <- list()
rank_summary_list <- list()
source_taxid_list <- list()
taxid_counts_list <- list()

# Process each file individually
for(file_idx in seq_along(score_files)) {
    file <- score_files[file_idx]
    samp_name <- sub("_+scores.output$", "", basename(file))
    # Extract sample ID (remove database suffix, should only be soil_microbe_db at this point)
    sampleID <- sub("_+soil_microbe_db.*$", "", samp_name)
    
    # OPTIMIZED: Aggregate first, join once
    # 1. Load only needed columns (saves massive RAM)
    tryCatch({
        scores_raw <- fread(file, 
                            select = c("taxid", "rank", "name"),
                            nThread = 4,
                            fill = TRUE,
                            showProgress = FALSE)
    }, error = function(e) {
        cat("ERROR reading file:", file, "\n")
        cat("Error message:", e$message, "\n")
        stop("Failed to read file: ", file)
    })
    
    # 2. AGGREGATE IMMEDIATELY: Count how many times each taxid appears
    # This turns millions of reads into thousands of unique TaxIDs
    sample_counts <- scores_raw[
        !is.na(taxid) & taxid > 0,
        .(read_count = .N),
        by = .(taxid, rank, name)
    ]
    
    # Clean up raw data immediately
    rm(scores_raw)
    gc()
    
    # 3. JOIN with weights table (handles proportional attribution automatically)
    # weights_dt has multiple rows per taxid (one per source), so this creates
    # the proportional expansion we need
    attributed <- merge(sample_counts, 
                        weights_dt, 
                        by.x = "taxid", 
                        by.y = "ncbi_species_taxid", 
                        allow.cartesian = TRUE)
    
    # 4. Calculate weighted reads (proportional attribution)
    attributed[, weighted_n := read_count * proportion]
    
    # 5. Identify kingdoms
    attributed[, kingdom_classified := fcase(
        grepl(fungal_pattern, name, perl = TRUE) | (!is.na(kingdom) & kingdom == "Fungi"), "Fungi",
        (rank == "k" & grepl("Bacteria", name, fixed = TRUE)) |
        grepl("^k__Bacteria|;k__Bacteria", name, fixed = FALSE) |
        (!is.na(kingdom) & kingdom == "Bacteria"), "Bacteria",
        (rank == "k" & grepl("Archaea", name, fixed = TRUE)) |
        grepl("^k__Archaea|;k__Archaea", name, fixed = FALSE) |
        (!is.na(kingdom) & kingdom == "Archaea"), "Archaea",
        default = "Other"
    )]
    
    # 6. REDUCE: Collapse to summaries immediately
    # Summary by source and kingdom
    samp_summary <- attributed[
        !is.na(source),
        .(n_reads = sum(weighted_n, na.rm = TRUE)),
        by = .(source, kingdom_classified)
    ]
    samp_summary[, sampleID := sampleID]
    summary_list[[file_idx]] <- samp_summary
    
    # Summary by source and rank
    samp_rank <- attributed[
        !is.na(source) & !is.na(rank),
        .(n_reads = sum(weighted_n, na.rm = TRUE)),
        by = .(source, rank)
    ]
    samp_rank[, sampleID := sampleID]
    rank_summary_list[[file_idx]] <- samp_rank
    
    # Track unique taxids per source per sample (for efficiency metrics)
    samp_taxids <- attributed[
        !is.na(source) & !is.na(taxid) & weighted_n > 0,
        .(unique_taxids = list(unique(taxid))),
        by = source
    ]
    samp_taxids[, sampleID := sampleID]
    source_taxid_list[[file_idx]] <- samp_taxids
    
    # Track taxid counts per source for unique contributions
    samp_taxid_counts <- attributed[
        !is.na(source) & !is.na(taxid),
        .(n_reads = sum(weighted_n, na.rm = TRUE)),
        by = .(source, taxid)
    ]
    samp_taxid_counts[, sampleID := sampleID]
    taxid_counts_list[[file_idx]] <- samp_taxid_counts
    
    # Clean up
    rm(sample_counts, attributed, samp_summary, samp_rank, samp_taxids, samp_taxid_counts)
    gc()
}

# Combine all summary tables (these are small, so this is fast)
final_summary_dt <- rbindlist(summary_list)
final_rank_dt <- rbindlist(rank_summary_list)
final_taxid_dt <- rbindlist(source_taxid_list)
final_taxid_counts_dt <- rbindlist(taxid_counts_list)

rm(summary_list, rank_summary_list, source_taxid_list, taxid_counts_list)
gc()

# ============================================================================
# AGGREGATE SUMMARIES
# Note: Values are AVERAGED across samples (not summed)
# ============================================================================

# Overall source distribution
# Calculate average reads per sample (not total)
all_sources <- final_summary_dt[
    ,
    .(
        n_reads_avg = mean(n_reads),
        n_reads_total = sum(n_reads),  # Keep total for percentage calculations
        n_unique_samples = uniqueN(sampleID)
    ),
    by = source
]
# Use average for main reporting
all_sources[, n_reads := n_reads_avg]

# Get unique taxids per source across all samples
source_taxids_all <- final_taxid_dt[
    ,
    .(all_taxids = list(unique(unlist(unique_taxids)))),
    by = source
]
source_taxids_all[, n_unique_taxids := lengths(all_taxids)]

all_sources <- all_sources[source_taxids_all[, .(source, n_unique_taxids)], on = "source"]
setorder(all_sources, -n_reads)

all_sources_file <- file.path(output_dir, "read_level_all_kingdoms_read_counts_by_source.csv")
fwrite(all_sources, all_sources_file)

# Efficiency metrics
# Note: n_reads is averaged across samples (may be fractional for overlapping taxids due to proportional attribution)
efficiency_metrics <- all_sources[, `:=`(
    reads_per_taxid = n_reads / n_unique_taxids,
    reads_per_sample = n_reads,  # Already averaged, so this is the same
    taxids_per_sample = n_unique_taxids / n_unique_samples
)]
setorder(efficiency_metrics, -n_reads)

efficiency_file <- file.path(output_dir, "read_level_database_efficiency_metrics.csv")
fwrite(efficiency_metrics[, .(source, n_reads, n_unique_taxids, n_unique_samples, 
                               reads_per_taxid, reads_per_sample, taxids_per_sample)], 
       efficiency_file)

# Taxonomic rank distribution
# Average reads per sample by rank
rank_totals <- final_rank_dt[
    ,
    .(n_reads = mean(n_reads)),
    by = .(source, rank)
]

source_totals_rank <- final_rank_dt[
    ,
    .(total_reads = mean(n_reads)),
    by = source
]

rank_totals <- rank_totals[source_totals_rank, on = "source"]
rank_totals[, pct_reads := (n_reads / total_reads) * 100]

rank_distribution_wide <- dcast(
    rank_totals,
    source ~ rank,
    value.var = "pct_reads",
    fill = 0
)

rank_cols <- setdiff(names(rank_distribution_wide), "source")
for(col in rank_cols) {
    rank_distribution_wide[[col]] <- round(rank_distribution_wide[[col]], 2)
}

setkey(rank_distribution_wide, source)
setkey(source_totals_rank, source)
rank_distribution_wide <- rank_distribution_wide[source_totals_rank, nomatch = 0]
setorder(rank_distribution_wide, -total_reads)
rank_distribution_wide[, total_reads := NULL]

rank_file <- file.path(output_dir, "read_level_database_rank_distribution.csv")
fwrite(rank_distribution_wide, rank_file)

# Cross-kingdom analysis
# Average reads per sample
cross_kingdom <- final_summary_dt[
    ,
    .(n_reads = mean(n_reads)),
    by = .(source, kingdom_classified)
]

source_totals <- final_summary_dt[
    ,
    .(total_reads = mean(n_reads)),
    by = source
]

cross_kingdom <- cross_kingdom[source_totals, on = "source"]
cross_kingdom[, pct_of_source_reads := (n_reads / total_reads) * 100]

cross_kingdom_wide <- dcast(
    cross_kingdom,
    source ~ kingdom_classified,
    value.var = "pct_of_source_reads",
    fill = 0
)

setorder(source_totals, -total_reads)
setkey(cross_kingdom_wide, source)
setkey(source_totals, source)
cross_kingdom_wide <- cross_kingdom_wide[source_totals, nomatch = 0]
setorder(cross_kingdom_wide, -total_reads)
cross_kingdom_wide[, total_reads := NULL]

kingdom_cols <- setdiff(names(cross_kingdom_wide), "source")
for(col in kingdom_cols) {
    cross_kingdom_wide[[col]] <- round(cross_kingdom_wide[[col]], 2)
}

cross_kingdom_file <- file.path(output_dir, "read_level_database_cross_kingdom.csv")
fwrite(cross_kingdom_wide, cross_kingdom_file)

# Kingdom-specific summaries (consolidated into single file)
# Average reads per sample by kingdom
# First average by source and kingdom, then sum across sources to get kingdom totals
kingdom_source_summary <- final_summary_dt[
    ,
    .(
        n_reads = mean(n_reads),
        n_unique_samples = uniqueN(sampleID)
    ),
    by = .(source, kingdom_classified)
]

# Calculate kingdom totals by summing the averaged source values
kingdom_totals <- kingdom_source_summary[
    ,
    .(total_reads_kingdom = sum(n_reads)),
    by = kingdom_classified
]

# Add percentage of reads within each kingdom
kingdom_source_summary <- kingdom_source_summary[
    kingdom_totals,
    on = "kingdom_classified"
]
kingdom_source_summary[, pct_of_kingdom_reads := (n_reads / total_reads_kingdom) * 100]
kingdom_source_summary[, total_reads_kingdom := NULL]

setorder(kingdom_source_summary, kingdom_classified, -n_reads)

kingdom_source_file <- file.path(output_dir, "read_level_read_counts_by_source_and_kingdom.csv")
fwrite(kingdom_source_summary, kingdom_source_file)

# Grouped comparison by kingdom
# Group sources into meaningful categories
final_summary_dt[, source_group := fcase(
    source == "GTDB", "GTDB",
    source %in% c("SMAG", "SPIRE_MAGs", "GEM catalog"), "MAG sources",
    source == "JGI GOLD", "JGI GOLD",
    source == "Mycocosm", "Mycocosm",
    source == "RefSoil", "RefSoil",
    default = "Other"
)]

kingdom_grouped_summary <- final_summary_dt[
    ,
    .(
        n_reads = mean(n_reads),
        n_unique_samples = uniqueN(sampleID)
    ),
    by = .(source_group, kingdom_classified)
]

# Calculate kingdom totals for grouped summary (sum of averaged group values)
kingdom_totals_grouped <- kingdom_grouped_summary[
    ,
    .(total_reads_kingdom = sum(n_reads)),
    by = kingdom_classified
]

kingdom_grouped_summary <- kingdom_grouped_summary[
    kingdom_totals_grouped,
    on = "kingdom_classified"
]
kingdom_grouped_summary[, pct_of_kingdom_reads := (n_reads / total_reads_kingdom) * 100]
kingdom_grouped_summary[, total_reads_kingdom := NULL]

setorder(kingdom_grouped_summary, kingdom_classified, -n_reads)

kingdom_grouped_file <- file.path(output_dir, "read_level_read_counts_by_source_group_and_kingdom.csv")
fwrite(kingdom_grouped_summary, kingdom_grouped_file)

# Unique contributions analysis (using genome_table, not raw reads)
taxid_source_count <- genome_table[
    !is.na(ncbi_species_taxid),
    .(n_sources = uniqueN(source)),
    by = ncbi_species_taxid
]

unique_taxids_only <- taxid_source_count[n_sources == 1, ncbi_species_taxid]

unique_taxid_source_map <- genome_table[
    ncbi_species_taxid %in% unique_taxids_only,
    .(unique_source = source[1]),
    by = ncbi_species_taxid
]
setkey(unique_taxid_source_map, ncbi_species_taxid)

# Get taxids actually used in classifications from our summaries
taxids_in_reads <- unique(unlist(source_taxids_all$all_taxids))

# Count reads from unique taxids per source (using aggregated taxid counts)
unique_contributions <- rbindlist(lapply(unique(efficiency_metrics$source), function(src) {
    unique_to_source <- unique_taxid_source_map[
        unique_source == src & ncbi_species_taxid %in% taxids_in_reads,
        ncbi_species_taxid
    ]
    
    if(length(unique_to_source) == 0) {
        return(data.table(
            source = src,
            n_unique_taxids = 0L,
            n_reads_classified_by_unique_taxids = 0L,
            pct_of_source_reads = 0
        ))
    }
    
    # Count reads for unique taxids from aggregated counts (average across samples)
    unique_reads <- final_taxid_counts_dt[
        source == src & taxid %in% unique_to_source,
        mean(n_reads)
    ]
    
    total_reads_source <- efficiency_metrics[source == src, n_reads]
    
    data.table(
        source = src,
        n_unique_taxids = length(unique_to_source),
        n_reads_classified_by_unique_taxids = unique_reads,
        pct_of_source_reads_from_unique_taxids = ifelse(total_reads_source > 0, 
                                                         (unique_reads / total_reads_source) * 100, 
                                                         0)
    )
}))

setorder(unique_contributions, -n_unique_taxids)

unique_file <- file.path(output_dir, "read_level_database_unique_contributions.csv")
fwrite(unique_contributions, unique_file)

# Overlap analysis (using genome_table, not raw reads)
major_sources <- efficiency_metrics[n_reads > 1000, source]
if(length(major_sources) > 1) {
    taxids_per_source <- genome_table[
        !is.na(ncbi_species_taxid) & source %in% major_sources,
        .(taxids = list(unique(ncbi_species_taxid))),
        by = source
    ]
    
    overlap_matrix <- matrix(0, nrow = length(major_sources), ncol = length(major_sources))
    rownames(overlap_matrix) <- major_sources
    colnames(overlap_matrix) <- major_sources
    
    for(i in seq_along(major_sources)) {
        src_i <- major_sources[i]
        taxids_i <- taxids_per_source[source == src_i, unlist(taxids)]
        
        for(j in seq_along(major_sources)) {
            if(i == j) {
                overlap_matrix[i, j] <- 100
            } else {
                src_j <- major_sources[j]
                taxids_j <- taxids_per_source[source == src_j, unlist(taxids)]
                
                if(length(taxids_i) > 0) {
                    overlap_pct <- (length(intersect(taxids_i, taxids_j)) / length(taxids_i)) * 100
                    overlap_matrix[i, j] <- round(overlap_pct, 2)
                }
            }
        }
    }
    
    overlap_dt <- as.data.table(overlap_matrix, keep.rownames = TRUE)
    setnames(overlap_dt, "rn", "source")
    
    overlap_file <- file.path(output_dir, "read_level_database_overlap_matrix.csv")
    fwrite(overlap_dt, overlap_file)
}
