#!/usr/bin/env Rscript
# Visualize database correctness based on BLAST verification
# Shows which databases (if any) were correct for each BLASTed read
#
# Usage: Rscript scripts/visualize/visualize_database_correctness.r [sampleID]
#
# Input:  {sampleID}_blast_vs_databases_comparison.csv (from script 19)
# Output: {sampleID}_database_correctness.png

library(tidyverse)
library(ggplot2)
library(patchwork)

# Sample to analyze
sampleID <- if(length(commandArgs(trailingOnly = TRUE)) > 0) {
    commandArgs(trailingOnly = TRUE)[1]
} else {
    "ORNL_046-O-20170621-COMP"
}

cat("=== Visualizing Database Correctness ===\n")
cat("Sample:", sampleID, "\n\n")

# File paths
input_file <- file.path("data/classification/analysis_files", 
                        paste0(sampleID, "_blast_vs_databases_comparison.csv"))
output_dir <- "manuscript_figures"

if(!file.exists(input_file)) {
    stop("Comparison file not found: ", input_file, 
         "\nPlease run script 19 first: Rscript scripts/summarize_outputs/19_compare_blast_to_databases.r ", sampleID)
}

cat("Reading BLAST vs databases comparison...\n")
data <- read_csv(input_file, show_col_types = FALSE)

# Remove header row if present
data <- data %>% filter(read_id != "read_id")

# Categorize correctness patterns
data <- data %>%
    mutate(
        correctness_pattern = case_when(
            gtdb_correct & smdb_correct & pluspf_correct ~ "All 3 correct",
            gtdb_correct & smdb_correct ~ "GTDB + SMD correct",
            gtdb_correct & pluspf_correct ~ "GTDB + PlusPF correct",
            smdb_correct & pluspf_correct ~ "SMD + PlusPF correct",
            gtdb_correct ~ "Only GTDB correct",
            smdb_correct ~ "Only SMD correct",
            pluspf_correct ~ "Only PlusPF correct",
            TRUE ~ "None correct"
        ),
        conflict_type_label = case_when(
            conflict_type == "homo_sapiens" ~ "PlusPF: Homo sapiens\nvs\nSMD: Fungi/Eukaryota",
            conflict_type == "fungal_bacteria" ~ "SMD: Fungi\nvs\nGTDB: Bacteria",
            TRUE ~ conflict_type
        )
    )

# Create summary by conflict type
summary_data <- data %>%
    group_by(conflict_type_label, correctness_pattern) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(conflict_type_label) %>%
    mutate(
        pct = n / sum(n) * 100,
        total = sum(n)
    ) %>%
    ungroup()

# Order correctness patterns
pattern_order <- c("All 3 correct", "GTDB + SMD correct", "GTDB + PlusPF correct", 
                  "SMD + PlusPF correct", "Only GTDB correct", "Only SMD correct", 
                  "Only PlusPF correct", "None correct")
summary_data$correctness_pattern <- factor(summary_data$correctness_pattern, 
                                           levels = pattern_order)

# Create visualization
p1 <- ggplot(summary_data, aes(x = conflict_type_label, y = n, fill = correctness_pattern)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    geom_text(aes(label = ifelse(n > 0, paste0(n, "\n(", round(pct, 1), "%)"), "")),
             position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(
        name = "Database Correctness",
        values = c(
            "All 3 correct" = "#1B9E77",
            "GTDB + SMD correct" = "#66A61E",
            "GTDB + PlusPF correct" = "#E6AB02",
            "SMD + PlusPF correct" = "#A6761D",
            "Only GTDB correct" = "#7570B3",
            "Only SMD correct" = "#E7298A",
            "Only PlusPF correct" = "#D95F02",
            "None correct" = "#E31A1C"
        )
    ) +
    labs(
        x = "Conflict Type",
        y = "Number of Reads",
        title = "Database Correctness Based on BLAST Verification",
        subtitle = paste0("Sample: ", sampleID, " | Total reads with valid BLAST hits: ", 
                         sum(summary_data$n) / length(unique(summary_data$conflict_type_label)))
    ) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.text.x = element_text(size = 10)
    ) +
    guides(fill = guide_legend(nrow = 4, byrow = TRUE))

# Also create a summary table visualization
summary_table <- data %>%
    group_by(conflict_type_label) %>%
    summarize(
        total = n(),
        none_correct = sum(!gtdb_correct & !smdb_correct & !pluspf_correct, na.rm = TRUE),
        gtdb_only = sum(gtdb_correct & !smdb_correct & !pluspf_correct, na.rm = TRUE),
        smdb_only = sum(!gtdb_correct & smdb_correct & !pluspf_correct, na.rm = TRUE),
        pluspf_only = sum(!gtdb_correct & !smdb_correct & pluspf_correct, na.rm = TRUE),
        multiple = sum((gtdb_correct + smdb_correct + pluspf_correct) > 1, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(
        none_correct_pct = none_correct / total * 100,
        at_least_one_pct = (total - none_correct) / total * 100
    )

cat("\n=== Summary Statistics ===\n")
print(summary_table)

# Save visualization
output_file <- file.path(output_dir, paste0(sampleID, "_database_correctness.png"))
ggsave(output_file, p1, width = 10, height = 6, dpi = 300)
cat("\nSaved visualization to:", output_file, "\n")

cat("\n=== Key Findings ===\n")
for(i in 1:nrow(summary_table)) {
    cat("\n", summary_table$conflict_type_label[i], ":\n", sep = "")
    cat("  Total reads:", summary_table$total[i], "\n")
    cat("  None correct:", summary_table$none_correct[i], 
        sprintf("(%.1f%%)\n", summary_table$none_correct_pct[i]))
    cat("  At least one correct:", summary_table$total[i] - summary_table$none_correct[i], 
        sprintf("(%.1f%%)\n", summary_table$at_least_one_pct[i]))
}

cat("\n=== Visualization complete ===\n")
