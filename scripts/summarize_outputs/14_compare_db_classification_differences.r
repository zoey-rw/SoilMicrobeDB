#!/usr/bin/env Rscript
# 14: Compare classification differences between gtdb_207_unfiltered and soil_microbe_db
# Analyzes differences in classification percentages, filter results, and quality metrics
#
# Usage: Rscript scripts/summarize_outputs/14_compare_db_classification_differences.r
#
# Input:  filter_results_summary.csv (from 04_reshape_score_reads.r)
# Output: db_classification_comparison.csv, db_classification_comparison_summary.txt

library(tidyverse)

filter_results <- read_csv("data/classification/analysis_files/filter_results_summary.csv", 
                          show_col_types = FALSE)

target_dbs <- c("gtdb_207_unfiltered", "soil_microbe_db")
db_data <- filter_results %>%
    filter(db_name %in% target_dbs)

db_wide <- db_data %>%
    select(sampleID, db_name, metric, value) %>%
    pivot_wider(names_from = c(db_name, metric), 
                values_from = value,
                names_sep = "_")

metrics <- c("percent_classified", "percent_passing", "mean_consistency", 
             "mean_multiplicity", "mean_entropy", "mean_confidence", 
             "n_total_classified_reads", "n_reads", "pct_of_classified_passing")

comparison <- db_wide %>%
    select(sampleID, 
           contains("gtdb_207_unfiltered"), 
           contains("soil_microbe_db")) %>%
    filter(!is.na(gtdb_207_unfiltered_percent_classified) & 
           !is.na(soil_microbe_db_percent_classified))

excluded_samples <- comparison %>%
    filter(gtdb_207_unfiltered_percent_classified < 0.01) %>%
    pull(sampleID)

if(length(excluded_samples) > 0) {
    comparison <- comparison %>%
        filter(!sampleID %in% excluded_samples)
}

if("gtdb_207_unfiltered_percent_classified" %in% names(comparison) && 
   "soil_microbe_db_percent_classified" %in% names(comparison)) {
    
    comparison <- comparison %>%
        mutate(
            diff_percent_classified = gtdb_207_unfiltered_percent_classified - soil_microbe_db_percent_classified,
            ratio_percent_classified = gtdb_207_unfiltered_percent_classified / soil_microbe_db_percent_classified
        )
}

if("gtdb_207_unfiltered_n_total_classified_reads" %in% names(comparison) && 
   "soil_microbe_db_n_total_classified_reads" %in% names(comparison)) {
    
    comparison <- comparison %>%
        mutate(
            diff_n_classified = gtdb_207_unfiltered_n_total_classified_reads - soil_microbe_db_n_total_classified_reads,
            ratio_n_classified = gtdb_207_unfiltered_n_total_classified_reads / soil_microbe_db_n_total_classified_reads
        )
}

if("gtdb_207_unfiltered_percent_passing" %in% names(comparison) && 
   "soil_microbe_db_percent_passing" %in% names(comparison)) {
    
    comparison <- comparison %>%
        mutate(
            diff_percent_passing = gtdb_207_unfiltered_percent_passing - soil_microbe_db_percent_passing
        )
}

gtdb_higher <- comparison %>%
    filter(diff_percent_classified > 0) %>%
    arrange(desc(diff_percent_classified))

smdb_higher <- comparison %>%
    filter(diff_percent_classified < 0) %>%
    arrange(diff_percent_classified)

if("gtdb_207_unfiltered_seq_depth" %in% names(comparison)) {
    comparison <- comparison %>%
        mutate(seq_depth_category = case_when(
            gtdb_207_unfiltered_seq_depth < 1e6 ~ "< 1M",
            gtdb_207_unfiltered_seq_depth < 5e6 ~ "1-5M",
            gtdb_207_unfiltered_seq_depth < 10e6 ~ "5-10M",
            TRUE ~ "> 10M"
        ))
    
    depth_summary <- comparison %>%
        group_by(seq_depth_category) %>%
        summarize(
            n_samples = n(),
            mean_gtdb_pct = mean(gtdb_207_unfiltered_percent_classified, na.rm=TRUE) * 100,
            mean_smd_pct = mean(soil_microbe_db_percent_classified, na.rm=TRUE) * 100,
            mean_diff = mean(diff_percent_classified, na.rm=TRUE) * 100,
            n_gtdb_higher = sum(diff_percent_classified > 0, na.rm=TRUE),
            .groups = "drop"
        )
}

output_file <- "data/classification/analysis_files/db_classification_comparison.csv"
write_csv(comparison, output_file)

summary_file <- "data/classification/analysis_files/db_classification_comparison_summary.txt"
summary_text <- paste0(
    "Database Classification Comparison: GTDB 207 unfiltered vs SoilMicrobeDB\n",
    "======================================================================\n\n",
    "Total samples (after exclusions): ", nrow(comparison), "\n",
    "Excluded samples (GTDB < 1%): ", length(excluded_samples), "\n\n",
    "=== Overall Statistics ===\n",
    "GTDB 207 unfiltered - mean percent_classified: ",
    sprintf("%.4f%%", mean(comparison$gtdb_207_unfiltered_percent_classified, na.rm=TRUE) * 100), "\n",
    "SoilMicrobeDB - mean percent_classified: ",
    sprintf("%.4f%%", mean(comparison$soil_microbe_db_percent_classified, na.rm=TRUE) * 100), "\n",
    "Mean difference (GTDB - SMD): ",
    sprintf("%.4f%%", mean(comparison$diff_percent_classified, na.rm=TRUE) * 100), "\n",
    "Mean ratio (GTDB / SMD): ",
    sprintf("%.2fx", mean(comparison$ratio_percent_classified, na.rm=TRUE)), "\n\n",
    "Samples where GTDB > SMD: ", nrow(gtdb_higher), "\n",
    "Samples where SMD > GTDB: ", nrow(smdb_higher), "\n\n",
    "=== Pattern by Sequencing Depth ===\n",
    if(exists("depth_summary")) {
        paste(capture.output(print(depth_summary)), collapse = "\n")
    } else {
        "Sequencing depth data not available"
    },
    "\n\n",
    "Detailed comparison saved to: ", output_file, "\n"
)

writeLines(summary_text, summary_file)
