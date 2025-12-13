# Calculate classification percentage at each taxonomic rank
# Shows how classification % changes as taxonomic resolution increases
# Output: Table for supplement showing classification % at phylum → class → order → family → genus → species → strain

library(tidyverse)
library(data.table)

source("scripts/helper_functions.r")

# Read sequencing depth for normalization
seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds")

# Define fungal phyla for filtering
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")

# Read rank-specific files (similar approach to calculate_pct_phylum_without_genus.r)
# Each file contains reads classified at that rank or higher
cat("Reading phylum data...\n")
phylum_data <- fread("data/classification/taxonomic_rank_summaries/phylum/soil_microbe_db_filtered_phylum_merged_lineage.csv", nThread = 8)

cat("Reading genus data...\n")
genus_data <- fread("data/classification/taxonomic_rank_summaries/genus/soil_microbe_db_filtered_genus_merged_lineage.csv", nThread = 8)

cat("Reading species data (this may take a moment - file is large)...\n")
species_data <- fread("data/classification/taxonomic_rank_summaries/species/soil_microbe_db_filtered_species_merged_lineage.csv", nThread = 8)
cat("Species data loaded. Parsing sample IDs and lineage...\n")

# Parse sample IDs for all datasets using data.table for speed
parse_sample_ids_dt <- function(df) {
    setDT(df)
    df[, c("sampleID", "db_name") := tstrsplit(sample_id, "COMP_", fixed = TRUE, keep = 1:2)]
    df[, sampleID := sub(paste0("_", db_name), "", sample_id, fixed = TRUE)]
    df[, db_name := gsub("_filtered|_genus_filtered|_phylum_filtered", "", db_name)]
    as_tibble(df)
}

phylum_data <- parse_sample_ids_dt(phylum_data)
genus_data <- parse_sample_ids_dt(genus_data)
cat("Parsing sample IDs for species data...\n")
species_data <- parse_sample_ids_dt(species_data)
cat("Sample ID parsing complete.\n")

# Identify fungi in each dataset
phylum_data <- phylum_data %>%
    mutate(is_fungi = name %in% fungal_phyla)

genus_data <- genus_data %>%
    mutate(is_fungi = grepl(paste(fungal_phyla, collapse = "|"), lineage))

# Parse all ranks from species lineage once using data.table for speed
cat("Parsing taxonomic ranks from species lineage (this may take a moment)...\n")
setDT(species_data)
species_data[, is_fungi := grepl(paste(fungal_phyla, collapse = "|"), lineage)]
species_data[, has_class := grepl("c__[^;]+", lineage) & !grepl("c__;", lineage)]
species_data[, has_order := grepl("o__[^;]+", lineage) & !grepl("o__;", lineage)]
species_data[, has_family := grepl("f__[^;]+", lineage) & !grepl("f__;", lineage)]
species_data[, has_strain := grepl("s1__[^;]+", lineage)]
species_data <- as_tibble(species_data)
cat("Lineage parsing complete.\n")

# Function to calculate classification % at a specific rank
# Uses the appropriate rank-level file (phylum, genus, or species)
# For intermediate ranks (class, order, family), parses lineage from species file
calculate_pct_at_rank <- function(data, rank_name, filter_fungi = FALSE, use_lineage_parsing = FALSE) {
    if (filter_fungi) {
        data <- data %>% filter(is_fungi == TRUE)
    }
    
    if (use_lineage_parsing) {
        # Use pre-parsed flags for efficiency (already done when loading species_data)
        if (rank_name == "class") {
            data <- data %>% filter(has_class == TRUE)
        } else if (rank_name == "order") {
            data <- data %>% filter(has_order == TRUE)
        } else if (rank_name == "family") {
            data <- data %>% filter(has_family == TRUE)
        } else if (rank_name == "strain") {
            data <- data %>% filter(has_strain == TRUE)
        }
    }
    
    # Sum reads classified at this rank (file contains reads at this rank or higher)
    classified <- data %>%
        group_by(sampleID, db_name) %>%
        summarize(
            reads_classified = sum(kraken_assigned_reads, na.rm = TRUE),
            .groups = "drop"
        )
    
    # Merge with sequencing depth to calculate percentage
    summary_df <- seq_depth_df %>%
        select(sampleID, db_name, seq_depth) %>%
        left_join(classified, by = c("sampleID", "db_name")) %>%
        mutate(
            reads_classified = replace_na(reads_classified, 0),
            pct_classified = (reads_classified / seq_depth) * 100,
            taxonomic_rank = rank_name
        )
    
    return(summary_df)
}

# ============================================================================
# All-domain classification by rank
# ============================================================================

# Calculate for each rank
# Use rank-specific files where available, parse lineage for intermediate ranks
all_domain_results <- bind_rows(
    calculate_pct_at_rank(phylum_data, "phylum", filter_fungi = FALSE),
    calculate_pct_at_rank(species_data, "class", filter_fungi = FALSE, use_lineage_parsing = TRUE),
    calculate_pct_at_rank(species_data, "order", filter_fungi = FALSE, use_lineage_parsing = TRUE),
    calculate_pct_at_rank(species_data, "family", filter_fungi = FALSE, use_lineage_parsing = TRUE),
    calculate_pct_at_rank(genus_data, "genus", filter_fungi = FALSE),
    calculate_pct_at_rank(species_data, "species", filter_fungi = FALSE),
    calculate_pct_at_rank(species_data, "strain", filter_fungi = FALSE, use_lineage_parsing = TRUE)
)

# Aggregate by database and rank
all_domain_summary <- all_domain_results %>%
    group_by(db_name, taxonomic_rank) %>%
    summarize(
        mean_pct_classified = mean(pct_classified, na.rm = TRUE),
        median_pct_classified = median(pct_classified, na.rm = TRUE),
        sd_pct_classified = sd(pct_classified, na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
    ) %>%
    mutate(taxon_group = "all_domain") %>%
    select(taxon_group, db_name, taxonomic_rank, mean_pct_classified, median_pct_classified, sd_pct_classified, n_samples)

# ============================================================================
# Fungi-specific classification by rank
# ============================================================================

# Calculate for fungi at each rank
fungi_results <- bind_rows(
    calculate_pct_at_rank(phylum_data, "phylum", filter_fungi = TRUE),
    calculate_pct_at_rank(species_data, "class", filter_fungi = TRUE, use_lineage_parsing = TRUE),
    calculate_pct_at_rank(species_data, "order", filter_fungi = TRUE, use_lineage_parsing = TRUE),
    calculate_pct_at_rank(species_data, "family", filter_fungi = TRUE, use_lineage_parsing = TRUE),
    calculate_pct_at_rank(genus_data, "genus", filter_fungi = TRUE),
    calculate_pct_at_rank(species_data, "species", filter_fungi = TRUE),
    calculate_pct_at_rank(species_data, "strain", filter_fungi = TRUE, use_lineage_parsing = TRUE)
)

# Aggregate by database and rank
fungi_summary <- fungi_results %>%
    group_by(db_name, taxonomic_rank) %>%
    summarize(
        mean_pct_classified = mean(pct_classified, na.rm = TRUE),
        median_pct_classified = median(pct_classified, na.rm = TRUE),
        sd_pct_classified = sd(pct_classified, na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
    ) %>%
    mutate(taxon_group = "fungi") %>%
    select(taxon_group, db_name, taxonomic_rank, mean_pct_classified, median_pct_classified, sd_pct_classified, n_samples)

# ============================================================================
# Combine and save results
# ============================================================================

classification_by_rank <- bind_rows(all_domain_summary, fungi_summary) %>%
    arrange(taxon_group, db_name, taxonomic_rank)

write_csv(classification_by_rank, "data/classification/analysis_files/classification_pct_by_rank.csv")
cat("✅ Saved classification percentages by rank to: data/classification/analysis_files/classification_pct_by_rank.csv\n")

# Also save per-sample results for detailed analysis
all_results <- bind_rows(
    all_domain_results %>% mutate(taxon_group = "all_domain"),
    fungi_results %>% mutate(taxon_group = "fungi")
)

write_csv(all_results, "data/classification/analysis_files/classification_pct_by_rank_per_sample.csv")
cat("✅ Saved per-sample classification percentages by rank to: data/classification/analysis_files/classification_pct_by_rank_per_sample.csv\n")
