#!/usr/bin/env Rscript
# 20: Analyze read-level database contributions
# Classify reads by source using read-level scoring files
# Goal: Match each read's taxid with genome table to get source database
#       This gives us the actual source used for each read classification
#
# Usage: Rscript scripts/summarize_outputs/20_analyze_read_level_database_contributions.r
#
# Input:  *_scores.output files from 02_bracken_output
#         Genome table with source information
# Output: Summary statistics showing read counts by source for each kingdom
#         Efficiency metrics, rank distribution, unique contributions,
#         overlap analysis, and cross-kingdom contributions

library(data.table)

# File paths
scores_dir <- "data/classification/02_bracken_output"
genome_table_file <- "data/genome_database/soil_microbe_db_genome_table.csv"
output_dir <- "data/classification/analysis_files"

# Fungal phyla for identification
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")
fungal_pattern <- paste0("(", paste(fungal_phyla, collapse = "|"), ")")

# Read genome table
genome_table <- fread(genome_table_file, nThread = 8)

# Handle many-to-many: if a taxid has multiple sources, take the most common one
source_mapping <- genome_table[
    !is.na(ncbi_species_taxid),
    .(source = names(sort(table(source), decreasing = TRUE))[1],
      kingdom = names(sort(table(kingdom), decreasing = TRUE))[1]),
    by = ncbi_species_taxid
]
setkey(source_mapping, ncbi_species_taxid)

# Find and process scoring files (recursively search subdirectories)
score_files <- list.files(scores_dir, pattern = "_+scores.output$", full.names = TRUE, recursive = TRUE)

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
    sampleID <- sub("_+(soil_microbe_db|pluspf|gtdb_207_unfiltered|gtdb_207)$", "", samp_name)
    
    # Read file - use fill=Inf to handle lines with extra fields
    # Expected format: sample_id,read_id,taxid,name,rank,n_kmers,consistency,confidence,multiplicity,entropy (10 fields)
    # 
    # NOTE: Rare data corruption: Some files may contain isolated corrupted lines where
    # multiple taxonomic assignments are concatenated (15+ fields instead of 10).
    # Investigation shows:
    # - All files have identical 10-column headers
    # - Only 1 line in 1 file (KONZ) has extra fields (data corruption, not format difference)
    # - The corrupted line contains two complete assignments concatenated
    # 
    # We handle this robustly by:
    # 1. Using fill=Inf to read all fields without error
    # 2. Extracting only the first assignment (columns: taxid, name, rank)
    # 3. Ignoring any duplicate assignments in extra columns
    scores_raw <- fread(file, 
                        nThread = 4, 
                        fill = Inf,
                        showProgress = FALSE)
    
    # Extract only the first taxonomic assignment (use named columns)
    # For corrupted lines with multiple assignments, we only use the first one
    if(all(c("taxid", "name", "rank") %in% names(scores_raw))) {
        scores <- scores_raw[, .(taxid = as.numeric(taxid), name = name, rank = rank)]
    } else {
        # Fallback: use positional columns if names don't match
        if(ncol(scores_raw) >= 5) {
            scores <- scores_raw[, .(taxid = as.numeric(get(names(scores_raw)[3])), 
                                name = get(names(scores_raw)[4]), 
                                rank = get(names(scores_raw)[5]))]
        } else {
            stop("Cannot parse file: ", file)
        }
    }
    
    # Remove rows with invalid taxids (NA or 0)
    scores <- scores[!is.na(taxid) & taxid > 0]
    
    # Clean up
    rm(scores_raw)
    
    # Join with source_mapping immediately (keeps memory footprint small)
    setkey(scores, taxid)
    scores <- source_mapping[scores, on = .(ncbi_species_taxid = taxid)]
    setnames(scores, "ncbi_species_taxid", "taxid")
    
    # Identify kingdoms (inside loop to process smaller chunks)
    scores[, kingdom_classified := fcase(
        grepl(fungal_pattern, name, perl = TRUE) | (!is.na(kingdom) & kingdom == "Fungi"), "Fungi",
        (rank == "k" & grepl("Bacteria", name, fixed = TRUE)) |
        grepl("^k__Bacteria|;k__Bacteria", name, fixed = FALSE) |
        (!is.na(kingdom) & kingdom == "Bacteria"), "Bacteria",
        (rank == "k" & grepl("Archaea", name, fixed = TRUE)) |
        grepl("^k__Archaea|;k__Archaea", name, fixed = FALSE) |
        (!is.na(kingdom) & kingdom == "Archaea"), "Archaea",
        default = "Other"
    )]
    
    # REDUCE: Collapse to summaries immediately (turns millions of rows into ~10-100 rows)
    # Summary by source and kingdom
    samp_summary <- scores[
        !is.na(source),
        .(n_reads = .N),
        by = .(source, kingdom_classified)
    ]
    samp_summary[, sampleID := sampleID]
    summary_list[[file_idx]] <- samp_summary
    
    # Summary by source and rank
    samp_rank <- scores[
        !is.na(source) & !is.na(rank),
        .(n_reads = .N),
        by = .(source, rank)
    ]
    samp_rank[, sampleID := sampleID]
    rank_summary_list[[file_idx]] <- samp_rank
    
    # Track unique taxids per source per sample (for efficiency metrics)
    samp_taxids <- scores[
        !is.na(source) & !is.na(taxid),
        .(unique_taxids = list(unique(taxid))),
        by = source
    ]
    samp_taxids[, sampleID := sampleID]
    source_taxid_list[[file_idx]] <- samp_taxids
    
    # Track taxid counts per source for unique contributions (aggregate immediately)
    samp_taxid_counts <- scores[
        !is.na(source) & !is.na(taxid),
        .(n_reads = .N),
        by = .(source, taxid)
    ]
    samp_taxid_counts[, sampleID := sampleID]
    taxid_counts_list[[file_idx]] <- samp_taxid_counts
    
    # Clean up
    rm(scores, samp_summary, samp_rank, samp_taxids, samp_taxid_counts)
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
# ============================================================================

# Overall source distribution
all_sources <- final_summary_dt[
    ,
    .(
        n_reads = sum(n_reads),
        n_unique_samples = uniqueN(sampleID)
    ),
    by = source
]

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
efficiency_metrics <- all_sources[, `:=`(
    reads_per_taxid = n_reads / n_unique_taxids,
    reads_per_sample = n_reads / n_unique_samples,
    taxids_per_sample = n_unique_taxids / n_unique_samples
)]
setorder(efficiency_metrics, -n_reads)

efficiency_file <- file.path(output_dir, "read_level_database_efficiency_metrics.csv")
fwrite(efficiency_metrics[, .(source, n_reads, n_unique_taxids, n_unique_samples, 
                               reads_per_taxid, reads_per_sample, taxids_per_sample)], 
       efficiency_file)

# Taxonomic rank distribution
rank_totals <- final_rank_dt[
    ,
    .(n_reads = sum(n_reads)),
    by = .(source, rank)
]

source_totals_rank <- final_rank_dt[
    ,
    .(total_reads = sum(n_reads)),
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
cross_kingdom <- final_summary_dt[
    ,
    .(n_reads = sum(n_reads)),
    by = .(source, kingdom_classified)
]

source_totals <- final_summary_dt[
    ,
    .(total_reads = sum(n_reads)),
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

# Kingdom-specific summaries
for(kingdom_name in c("Fungi", "Bacteria", "Archaea")) {
    kingdom_data <- final_summary_dt[kingdom_classified == kingdom_name]
    
    if(nrow(kingdom_data) == 0) next
    
    total_reads <- sum(kingdom_data$n_reads)
    
    source_summary <- kingdom_data[
        ,
        .(
            n_reads = sum(n_reads),
            n_unique_samples = uniqueN(sampleID),
            pct_reads = (sum(n_reads) / total_reads) * 100
        ),
        by = source
    ]
    setorder(source_summary, -n_reads)
    
    kingdom_lower <- tolower(kingdom_name)
    source_file <- file.path(output_dir, paste0("read_level_", kingdom_lower, "_read_counts_by_source.csv"))
    fwrite(source_summary, source_file)
    
    # Grouped comparison
    kingdom_data[, source_group := fcase(
        source == "Mycocosm", "Mycocosm",
        source == "GTDB", "GTDB",
        source == "JGI GOLD", "JGI GOLD",
        source == "SPIRE_MAGs", "SPIRE MAGs",
        source == "SMAG", "SMAG",
        source == "GEM catalog", "GEM catalog",
        source == "RefSoil", "RefSoil",
        default = "Other"
    )]
    
    comparison <- kingdom_data[
        ,
        .(
            n_reads = sum(n_reads),
            n_unique_samples = uniqueN(sampleID),
            pct_reads = (sum(n_reads) / total_reads) * 100
        ),
        by = source_group
    ]
    setorder(comparison, -n_reads)
    
    comparison_file <- file.path(output_dir, paste0("read_level_", kingdom_lower, "_read_counts_by_source_group.csv"))
    fwrite(comparison, comparison_file)
}

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
    
    # Count reads for unique taxids from aggregated counts
    unique_reads <- final_taxid_counts_dt[
        source == src & taxid %in% unique_to_source,
        sum(n_reads)
    ]
    
    total_reads_source <- efficiency_metrics[source == src, n_reads]
    
    data.table(
        source = src,
        n_unique_taxids = length(unique_to_source),
        n_reads_classified_by_unique_taxids = unique_reads,
        pct_of_source_reads = ifelse(total_reads_source > 0, 
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
