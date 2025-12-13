# Generate Supplementary Figure: Classification percentages by taxonomic rank
# Uses kreport approach (before and after Architeuthis filtering)
# Shows how classification % changes across taxonomic ranks and with filtering
# Note: This is a SUPPLEMENTARY figure - main Figure 2 compares across databases (compare_filter.r)

library(tidyverse)
library(ggrepel)
library(ggallin)
library(scales)
library(ggpmisc)
library(rstatix)
library(ggpubr)
library(ggstatsplot)
library(patchwork)
library(RColorBrewer)

# Load classification percentages by rank (from 06_calculate_classification_pct_by_rank.r)
classification_file = "data/NEON_metagenome_classification/summary_files/classification_pct_by_rank_per_sample.csv"
if(!file.exists(classification_file)) {
    # Fallback to old location
    classification_file = "data/classification/analysis_files/classification_pct_by_rank_per_sample.csv"
}

if(!file.exists(classification_file)) {
    stop("❌ MISSING FILE: Classification percentages by rank not found!\n",
         "   Expected: ", classification_file, "\n",
         "   Please run scripts/run_workflow/06_calculate_classification_pct_by_rank.r first.")
}

cat("Loading classification percentages from:", classification_file, "\n")
classification_data <- read_csv(classification_file, show_col_types = FALSE)
cat("Loaded", nrow(classification_data), "records\n")

# Filter to all-domain data only (exclude fungi-specific for main figure)
plot_data <- classification_data %>%
    filter(taxon_group == "all_domain") %>%
    mutate(
        filter_status_pretty = recode(filter_status,
            "before" = "Before filtering",
            "after" = "After filtering"),
        taxonomic_rank_pretty = factor(taxonomic_rank, 
            levels = c("phylum", "class", "order", "family", "genus", "species", "strain"),
            labels = c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"))
    )

# Check if we have data to plot
if(nrow(plot_data) == 0) {
    stop("❌ ERROR: No data available for plotting!\n",
         "   Make sure classification_pct_by_rank_per_sample.csv contains data.")
}

cat("Plotting", nrow(plot_data), "records\n")
cat("  Before filtering:", sum(plot_data$filter_status == "before"), "records\n")
cat("  After filtering:", sum(plot_data$filter_status == "after"), "records\n")

# Create summary statistics for plotting
summary_stats <- plot_data %>%
    group_by(taxonomic_rank_pretty, filter_status_pretty) %>%
    summarize(
        mean_pct = mean(pct_classified, na.rm = TRUE),
        median_pct = median(pct_classified, na.rm = TRUE),
        se_pct = sd(pct_classified, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
    )

# Plotting color palette
filter_colors <- c("Before filtering" = "#E31A1C", "After filtering" = "#1F78B4")
colScale <- scale_color_manual(name = "Filter status", values = filter_colors)
fillScale <- scale_fill_manual(name = "Filter status", values = filter_colors)

# Create main plot: Classification % by rank, comparing before/after filtering
fig2a <- ggplot(plot_data, aes(x = taxonomic_rank_pretty, y = pct_classified, 
                               color = filter_status_pretty, fill = filter_status_pretty)) +
    geom_boxplot(alpha = 0.3, outlier.alpha = 0.5, position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), 
               alpha = 0.4, size = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
                position = position_dodge(width = 0.75)) +
    colScale +
    fillScale +
    labs(
        x = "Taxonomic rank",
        y = "% of total sequencing depth classified\n(cumulative: rank or higher)",
        title = "A. Classification percentage by taxonomic rank"
    ) +
    theme_bw(base_size = 14) +
    theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1))

# Create summary line plot showing mean classification % across ranks
fig2b <- ggplot(summary_stats, aes(x = taxonomic_rank_pretty, y = mean_pct, 
                                   color = filter_status_pretty, group = filter_status_pretty)) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct), 
                  width = 0.2, alpha = 0.6) +
    colScale +
    labs(
        x = "Taxonomic rank",
        y = "Mean % classified (± SE)\n(cumulative: rank or higher)",
        title = "B. Mean classification percentage across ranks"
    ) +
    theme_bw(base_size = 14) +
    theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1))

# Combine plots
fig2 <- fig2a / fig2b + 
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = 'A')

# Save figure as supplementary figure
ggsave("manuscript_figures/fig_s2_classification_by_rank.png", fig2, width = 12, height = 10, units = "in", dpi = 300)

cat("✓ Supplementary figure saved to: manuscript_figures/fig_s2_classification_by_rank.png\n")
cat("   Note: Main Figure 2 compares databases (see compare_filter.r)\n")

# Print summary statistics
cat("\nSummary statistics:\n")
summary_table <- summary_stats %>%
    pivot_wider(names_from = filter_status_pretty, 
               values_from = c(mean_pct, median_pct, se_pct),
               names_sep = "_") %>%
    arrange(taxonomic_rank_pretty)
print(summary_table)
