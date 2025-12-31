library(tidyverse)
library(scales)
library(ggpmisc)
library(patchwork)
library(RColorBrewer)

classification_file <- if(file.exists("data/NEON_metagenome_classification/summary_files/classification_pct_by_rank_per_sample.csv")) {
    "data/NEON_metagenome_classification/summary_files/classification_pct_by_rank_per_sample.csv"
} else {
    "data/classification/analysis_files/classification_pct_by_rank_per_sample.csv"
}

classification_data <- read_csv(classification_file, show_col_types = FALSE)

rank_levels <- c("phylum", "class", "order", "family", "genus", "species", "strain")
rank_labels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
filter_colors <- c("Before filtering" = "#E31A1C", "After filtering" = "#1F78B4")

plot_data <- classification_data %>%
    filter(taxon_group == "all_domain") %>%
    mutate(
        filter_status_pretty = recode(filter_status,
            "before" = "Before filtering",
            "after" = "After filtering"),
        taxonomic_rank_pretty = factor(taxonomic_rank, levels = rank_levels, labels = rank_labels)
    )

summary_stats <- plot_data %>%
    group_by(taxonomic_rank_pretty, filter_status_pretty) %>%
    summarize(
        mean_pct = mean(pct_classified, na.rm = TRUE),
        se_pct = sd(pct_classified, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
    )

colScale <- scale_color_manual(name = "Filter status", values = filter_colors)
fillScale <- scale_fill_manual(name = "Filter status", values = filter_colors)

fig2a <- ggplot(plot_data, aes(x = taxonomic_rank_pretty, y = pct_classified, 
                               color = filter_status_pretty, fill = filter_status_pretty)) +
    geom_boxplot(alpha = 0.3, outlier.alpha = 0.5, position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), 
               alpha = 0.4, size = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
                position = position_dodge(width = 0.75)) +
    colScale + fillScale +
    labs(x = "Taxonomic rank",
         y = "% of total sequencing depth classified\n(cumulative: rank or higher)",
         title = "A. Classification percentage by taxonomic rank") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold")) +
    scale_y_continuous(labels = percent_format(scale = 1))

fig2b <- ggplot(summary_stats, aes(x = taxonomic_rank_pretty, y = mean_pct, 
                                   color = filter_status_pretty, group = filter_status_pretty)) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct), 
                  width = 0.2, alpha = 0.6) +
    colScale +
    labs(x = "Taxonomic rank",
         y = "Mean % classified (± SE)\n(cumulative: rank or higher)",
         title = "B. Mean classification percentage across ranks") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold")) +
    scale_y_continuous(labels = percent_format(scale = 1))

fig2 <- fig2a / fig2b + 
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = 'A')

ggsave("manuscript_figures/fig_s2_classification_by_rank.png", fig2, width = 12, height = 10, units = "in", dpi = 300)

