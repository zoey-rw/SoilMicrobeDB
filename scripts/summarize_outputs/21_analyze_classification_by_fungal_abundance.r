#!/usr/bin/env Rscript
# 21: Analyze how database classification rates change with fungal abundances
# Compares % classified and % passing filters across all databases
# as a function of PLFA-estimated fungal abundance
#
# Usage: Rscript scripts/summarize_outputs/21_analyze_classification_by_fungal_abundance.r
#
# Input:  filter_results_summary.csv (from 04_reshape_score_reads.r)
#         plfa_comparison.csv (from 10_merge_neon_plfa.r)
# Output: classification_by_fungal_abundance.csv (summary)
#         classification_by_fungal_abundance_per_sample.csv (detailed)

library(tidyverse)
library(data.table)

# Set up paths - check both new and old locations
# Use species-level classification percentages from 06_calculate_classification_pct_by_rank.r
classification_file <- "data/NEON_metagenome_classification/analysis_files/classification_pct_by_rank_per_sample.csv"
if(!file.exists(classification_file)) {
    classification_file <- "data/classification/analysis_files/classification_pct_by_rank_per_sample.csv"
}

plfa_file <- "data/classification/analysis_files/plfa_comparison.csv"

# Check required files
required_files <- c(classification_file, plfa_file)
missing_files <- required_files[!file.exists(required_files)]
if(length(missing_files) > 0) {
    stop("Required files not found:\n  ", paste(missing_files, collapse = "\n  "))
}

cat("Loading classification results from:", classification_file, "\n")
classification_data <- read_csv(classification_file, show_col_types = FALSE)

cat("Loading PLFA data from:", plfa_file, "\n")
plfa_data <- read_csv(plfa_file, show_col_types = FALSE)

# Check if plfa_data has proportion_fungi column
if(!"proportion_fungi" %in% names(plfa_data)) {
    stop("PLFA data missing 'proportion_fungi' column. Available columns:\n  ", 
         paste(names(plfa_data), collapse = ", "))
}

# Filter to species-level, after filtering, all_domain results
# Convert pct_classified from percentage to proportion (divide by 100)
classification_wide <- classification_data %>%
    filter(taxonomic_rank == "species",
           filter_status == "after",
           taxon_group == "all_domain") %>%
    mutate(
        pct_classified = pct_classified / 100,  # Convert from percentage to proportion
        pct_passing = pct_classified  # For consistency with old format (after filtering = passing)
    ) %>%
    select(sampleID, db_name, pct_classified, pct_passing)

cat("Classification results:\n")
cat("  Samples:", length(unique(classification_wide$sampleID)), "\n")
cat("  Databases:", paste(unique(classification_wide$db_name), collapse = ", "), "\n")
cat("  Total records:", nrow(classification_wide), "\n")

# Prepare PLFA data - extract compositeSampleID and proportion_fungi
# Also keep other potentially useful columns
plfa_subset <- plfa_data %>%
    select(compositeSampleID, proportion_fungi, 
           siteID, biome, horizon) %>%
    filter(!is.na(proportion_fungi), 
           !is.na(compositeSampleID)) %>%
    distinct(compositeSampleID, .keep_all = TRUE)

cat("PLFA data:\n")
cat("  Samples with fungal abundance:", nrow(plfa_subset), "\n")
cat("  Fungal abundance range:", 
    round(min(plfa_subset$proportion_fungi, na.rm = TRUE), 4), "-",
    round(max(plfa_subset$proportion_fungi, na.rm = TRUE), 4), "\n")

# Join classification results with PLFA data
# sampleID in filter_results should match compositeSampleID in PLFA data
# (both end with "-COMP")
combined_data <- classification_wide %>%
    left_join(plfa_subset, by = c("sampleID" = "compositeSampleID"))

# Count samples with and without PLFA data
samples_with_plfa <- sum(!is.na(combined_data$proportion_fungi))
samples_without_plfa <- sum(is.na(combined_data$proportion_fungi))

cat("\nMerged data:\n")
cat("  Samples with both classification and PLFA:", samples_with_plfa, "\n")
cat("  Samples with classification but no PLFA:", samples_without_plfa, "\n")

# Filter to samples with PLFA data for analysis
analysis_data <- combined_data %>%
    filter(!is.na(proportion_fungi),
           !is.na(pct_classified),
           !is.na(pct_passing))

cat("  Samples in final analysis:", nrow(analysis_data), "\n")

if(nrow(analysis_data) == 0) {
    stop("No samples with both classification and PLFA data. Cannot proceed.")
}

# Create summary statistics by database and fungal abundance bins
# Bin fungal abundance for easier visualization
analysis_data <- analysis_data %>%
    mutate(
        fungal_abundance_bin = cut(proportion_fungi,
                                   breaks = quantile(proportion_fungi, 
                                                    probs = seq(0, 1, by = 0.25),
                                                    na.rm = TRUE),
                                   include.lowest = TRUE,
                                   labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"))
    )

# Summary statistics by database and fungal abundance bin
summary_by_bin <- analysis_data %>%
    group_by(db_name, fungal_abundance_bin) %>%
    summarize(
        n_samples = n(),
        mean_fungal_abundance = mean(proportion_fungi, na.rm = TRUE),
        median_fungal_abundance = median(proportion_fungi, na.rm = TRUE),
        mean_pct_classified = mean(pct_classified, na.rm = TRUE),
        median_pct_classified = median(pct_classified, na.rm = TRUE),
        sd_pct_classified = sd(pct_classified, na.rm = TRUE),
        mean_pct_passing = mean(pct_passing, na.rm = TRUE),
        median_pct_passing = median(pct_passing, na.rm = TRUE),
        sd_pct_passing = sd(pct_passing, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(db_name, fungal_abundance_bin)

# Summary statistics by database (overall)
summary_by_db <- analysis_data %>%
    group_by(db_name) %>%
    summarize(
        n_samples = n(),
        mean_fungal_abundance = mean(proportion_fungi, na.rm = TRUE),
        median_fungal_abundance = median(proportion_fungi, na.rm = TRUE),
        mean_pct_classified = mean(pct_classified, na.rm = TRUE),
        median_pct_classified = median(pct_classified, na.rm = TRUE),
        sd_pct_classified = sd(pct_classified, na.rm = TRUE),
        mean_pct_passing = mean(pct_passing, na.rm = TRUE),
        median_pct_passing = median(pct_passing, na.rm = TRUE),
        sd_pct_passing = sd(pct_passing, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(db_name)

# Calculate correlation coefficients by database (absolute values)
correlations_absolute <- analysis_data %>%
    group_by(db_name) %>%
    summarize(
        n_samples = n(),
        cor_classified_fungi = cor(pct_classified, proportion_fungi, 
                                   use = "complete.obs"),
        cor_passing_fungi = cor(pct_passing, proportion_fungi, 
                               use = "complete.obs"),
        .groups = "drop"
    ) %>%
    arrange(db_name)

cat("\nCorrelations between classification rates and fungal abundance (absolute):\n")
print(correlations_absolute)

# Calculate improvements: SoilMicrobeDB vs other databases
# Pivot to wide format to compare databases
comparison_data <- analysis_data %>%
    select(sampleID, db_name, pct_classified, pct_passing, proportion_fungi, 
           siteID, biome, horizon) %>%
    pivot_wider(names_from = db_name, 
                values_from = c(pct_classified, pct_passing),
                names_sep = "_") %>%
    filter(!is.na(proportion_fungi))

# Calculate improvements (SoilMicrobeDB - other databases)
if("pct_classified_soil_microbe_db" %in% names(comparison_data)) {
    if("pct_classified_gtdb_207" %in% names(comparison_data)) {
        comparison_data <- comparison_data %>%
            mutate(
                improvement_classified_vs_gtdb = pct_classified_soil_microbe_db - pct_classified_gtdb_207,
                improvement_passing_vs_gtdb = pct_passing_soil_microbe_db - pct_passing_gtdb_207
            )
    }
    if("pct_classified_gtdb_207_unfiltered" %in% names(comparison_data)) {
        comparison_data <- comparison_data %>%
            mutate(
                improvement_classified_vs_gtdb_unfiltered = pct_classified_soil_microbe_db - pct_classified_gtdb_207_unfiltered,
                improvement_passing_vs_gtdb_unfiltered = pct_passing_soil_microbe_db - pct_passing_gtdb_207_unfiltered
            )
    }
    if("pct_classified_pluspf" %in% names(comparison_data)) {
        comparison_data <- comparison_data %>%
            mutate(
                improvement_classified_vs_pluspf = pct_classified_soil_microbe_db - pct_classified_pluspf,
                improvement_passing_vs_pluspf = pct_passing_soil_microbe_db - pct_passing_pluspf
            )
    }
}

# Calculate correlations for improvements
improvement_correlations <- tibble(
    comparison = character(),
    metric = character(),
    n_samples = integer(),
    correlation = numeric()
)

if("improvement_classified_vs_gtdb" %in% names(comparison_data)) {
    n_samples <- sum(!is.na(comparison_data$improvement_classified_vs_gtdb))
    if(n_samples > 0) {
        cor_val <- cor(comparison_data$improvement_classified_vs_gtdb, 
                      comparison_data$proportion_fungi, use = "complete.obs")
        improvement_correlations <- bind_rows(
            improvement_correlations,
            tibble(comparison = "SoilMicrobeDB vs GTDB 207",
                   metric = "percent_classified",
                   n_samples = n_samples,
                   correlation = cor_val)
        )
    }
    if("improvement_passing_vs_gtdb" %in% names(comparison_data)) {
        n_samples <- sum(!is.na(comparison_data$improvement_passing_vs_gtdb))
        if(n_samples > 0) {
            cor_val <- cor(comparison_data$improvement_passing_vs_gtdb, 
                          comparison_data$proportion_fungi, use = "complete.obs")
            improvement_correlations <- bind_rows(
                improvement_correlations,
                tibble(comparison = "SoilMicrobeDB vs GTDB 207",
                       metric = "percent_passing",
                       n_samples = n_samples,
                       correlation = cor_val)
            )
        }
    }
}

if("improvement_classified_vs_gtdb_unfiltered" %in% names(comparison_data)) {
    n_samples <- sum(!is.na(comparison_data$improvement_classified_vs_gtdb_unfiltered))
    if(n_samples > 0) {
        cor_val <- cor(comparison_data$improvement_classified_vs_gtdb_unfiltered, 
                      comparison_data$proportion_fungi, use = "complete.obs")
        improvement_correlations <- bind_rows(
            improvement_correlations,
            tibble(comparison = "SoilMicrobeDB vs GTDB 207 unfiltered",
                   metric = "percent_classified",
                   n_samples = n_samples,
                   correlation = cor_val)
        )
    }
    if("improvement_passing_vs_gtdb_unfiltered" %in% names(comparison_data)) {
        n_samples <- sum(!is.na(comparison_data$improvement_passing_vs_gtdb_unfiltered))
        if(n_samples > 0) {
            cor_val <- cor(comparison_data$improvement_passing_vs_gtdb_unfiltered, 
                          comparison_data$proportion_fungi, use = "complete.obs")
            improvement_correlations <- bind_rows(
                improvement_correlations,
                tibble(comparison = "SoilMicrobeDB vs GTDB 207 unfiltered",
                       metric = "percent_passing",
                       n_samples = n_samples,
                       correlation = cor_val)
            )
        }
    }
}

if("improvement_classified_vs_pluspf" %in% names(comparison_data)) {
    n_samples <- sum(!is.na(comparison_data$improvement_classified_vs_pluspf))
    if(n_samples > 0) {
        cor_val <- cor(comparison_data$improvement_classified_vs_pluspf, 
                      comparison_data$proportion_fungi, use = "complete.obs")
        improvement_correlations <- bind_rows(
            improvement_correlations,
            tibble(comparison = "SoilMicrobeDB vs PlusPF",
                   metric = "percent_classified",
                   n_samples = n_samples,
                   correlation = cor_val)
        )
    }
    if("improvement_passing_vs_pluspf" %in% names(comparison_data)) {
        n_samples <- sum(!is.na(comparison_data$improvement_passing_vs_pluspf))
        if(n_samples > 0) {
            cor_val <- cor(comparison_data$improvement_passing_vs_pluspf, 
                          comparison_data$proportion_fungi, use = "complete.obs")
            improvement_correlations <- bind_rows(
                improvement_correlations,
                tibble(comparison = "SoilMicrobeDB vs PlusPF",
                       metric = "percent_passing",
                       n_samples = n_samples,
                       correlation = cor_val)
            )
        }
    }
}

cat("\nCorrelations: SoilMicrobeDB improvement vs fungal abundance:\n")
cat("(Positive = SoilMicrobeDB gets better relative to other DBs as fungal abundance increases)\n")
if(nrow(improvement_correlations) > 0) {
    print(improvement_correlations)
} else {
    cat("  No comparison data available (need samples with multiple databases)\n")
}

# Save output files
output_dir <- "data/classification/analysis_files"
if(!dir.exists(output_dir)) {
    output_dir <- "data/NEON_metagenome_classification/summary_files"
    if(!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
}

# Save per-sample data
per_sample_output <- file.path(output_dir, "classification_by_fungal_abundance_per_sample.csv")
write_csv(analysis_data, per_sample_output)
cat("\n✅ Saved per-sample data to:", per_sample_output, "\n")
cat("   Records:", nrow(analysis_data), "\n")

# Save summary by bin
summary_bin_output <- file.path(output_dir, "classification_by_fungal_abundance_by_bin.csv")
write_csv(summary_by_bin, summary_bin_output)
cat("✅ Saved summary by fungal abundance bin to:", summary_bin_output, "\n")

# Save summary by database
summary_db_output <- file.path(output_dir, "classification_by_fungal_abundance_by_db.csv")
write_csv(summary_by_db, summary_db_output)
cat("✅ Saved summary by database to:", summary_db_output, "\n")

# Save correlations (absolute)
correlations_output <- file.path(output_dir, "classification_by_fungal_abundance_correlations.csv")
write_csv(correlations_absolute, correlations_output)
cat("✅ Saved correlations (absolute) to:", correlations_output, "\n")

# Save improvement correlations
if(nrow(improvement_correlations) > 0) {
    improvement_correlations_output <- file.path(output_dir, "classification_by_fungal_abundance_improvements.csv")
    write_csv(improvement_correlations, improvement_correlations_output)
    cat("✅ Saved improvement correlations to:", improvement_correlations_output, "\n")
}

# Save comparison data (improvements)
if(nrow(comparison_data) > 0) {
    comparison_output <- file.path(output_dir, "classification_by_fungal_abundance_comparison.csv")
    write_csv(comparison_data, comparison_output)
    cat("✅ Saved database comparison data to:", comparison_output, "\n")
    cat("   Records:", nrow(comparison_data), "\n")
}

cat("\n✅ Script completed successfully\n")
cat("\nSummary:\n")
cat("  Databases analyzed:", paste(unique(analysis_data$db_name), collapse = ", "), "\n")
cat("  Total samples:", nrow(analysis_data), "\n")
cat("  Fungal abundance range:", 
    round(min(analysis_data$proportion_fungi, na.rm = TRUE), 4), "-",
    round(max(analysis_data$proportion_fungi, na.rm = TRUE), 4), "\n")

