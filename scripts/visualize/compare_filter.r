library(tidyverse)
library(ggpubr)
library(ggstatsplot)
library(patchwork)
library(RColorBrewer)

# Read classification data for all ranks
# Use species-level data from the updated calculation script
classification_file <- "data/NEON_metagenome_classification/analysis_files/classification_pct_by_rank_per_sample.csv"
if(!file.exists(classification_file)) {
    classification_file <- "data/classification/analysis_files/classification_pct_by_rank_per_sample.csv"
}

cat("Loading classification data from:", classification_file, "\n")
summaries <- read_csv(classification_file, show_col_types = FALSE) %>%
    filter(taxon_group == "all_domain",
           !is.na(pct_classified))

# Define rank order from most general to most specific
rank_order <- c("domain", "phylum", "class", "order", "family", "genus", "species", "strain")
rank_labels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")

# Convert filter_status to metric format to match original structure
all_data <- summaries %>%
    mutate(
        metric = case_when(
            filter_status == "before" ~ "percent_classified",
            filter_status == "after" ~ "percent_passing",
            TRUE ~ NA_character_
        ),
        value = pct_classified,  # Already in percentage form
        Database = recode(db_name,
                         "gtdb_207" = "GTDB r207",
                         "gtdb_207_unfiltered" = "GTDB r207",
                         "soil_microbe_db" = "SoilMicrobeDB",
                         "pluspf" = "PlusPF"),
        pretty_metric = recode(metric,
                              "percent_classified" = "Read classification",
                              "percent_passing" = "Read classification and quality-filtering"),
        pretty_metric = factor(pretty_metric, levels = c("Read classification", "Read classification and quality-filtering")),
        taxonomic_rank = factor(taxonomic_rank, levels = rank_order)) %>%
    filter(db_name %in% c("gtdb_207", "gtdb_207_unfiltered", "pluspf", "soil_microbe_db")) %>%
    filter(!is.na(metric))

# Find common samples across all databases (using any rank to determine commonality)
pluspf_samples <- all_data %>% filter(db_name == "pluspf") %>% distinct(sampleID) %>% pull(sampleID)
smdb_samples <- all_data %>% filter(db_name == "soil_microbe_db") %>% distinct(sampleID) %>% pull(sampleID)
gtdb_unfiltered_samples <- all_data %>% filter(db_name == "gtdb_207_unfiltered") %>% distinct(sampleID) %>% pull(sampleID)
gtdb_filtered_samples <- all_data %>% filter(db_name == "gtdb_207") %>% distinct(sampleID) %>% pull(sampleID)

common_pluspf_smdb <- intersect(pluspf_samples, smdb_samples)
common_all_unfiltered <- intersect(intersect(pluspf_samples, smdb_samples), gtdb_unfiltered_samples)
common_all_filtered <- intersect(intersect(pluspf_samples, smdb_samples), gtdb_filtered_samples)

# Diagnostic: Show sample counts
cat("\nSample counts by database:\n")
cat("  PlusPF:", length(pluspf_samples), "samples\n")
cat("  SoilMicrobeDB:", length(smdb_samples), "samples\n")
cat("  GTDB r207 (unfiltered):", length(gtdb_unfiltered_samples), "samples\n")
cat("  Common to PlusPF + SoilMicrobeDB:", length(common_pluspf_smdb), "samples\n")
cat("  Common to all three:", length(common_all_unfiltered), "samples\n")
cat("  SoilMicrobeDB-only samples:", length(setdiff(smdb_samples, union(pluspf_samples, gtdb_unfiltered_samples))), "samples\n")

# Use unfiltered GTDB for this comparison (only unfiltered should be in the plot)
if(length(common_all_unfiltered) > 0) {
    gtdb_db_to_use <- "gtdb_207_unfiltered"
    common_samples <- common_all_unfiltered
    cat("\nUsing GTDB unfiltered with", length(common_samples), "common samples\n")
    cat("NOTE: Figure 2 only includes samples common to all three databases.\n")
    cat("      This excludes", length(setdiff(smdb_samples, common_samples)), "SoilMicrobeDB-only samples.\n")
} else {
    stop("No common samples found across all databases (PlusPF, SoilMicrobeDB, and GTDB unfiltered)")
}

# Filter to common samples and unfiltered GTDB
plot_data <- all_data %>%
    filter(db_name %in% c("gtdb_207_unfiltered", "pluspf", "soil_microbe_db")) %>%  # Only unfiltered, explicitly exclude gtdb_207
    filter(sampleID %in% common_samples) %>%
    mutate(Database = factor(Database, levels = c("SoilMicrobeDB", 
                                                   "GTDB r207",
                                                   "PlusPF")))

# Colorblind-friendly colors (Dark2 palette)
myColors <- c("SoilMicrobeDB" = "#1B9E77", 
              "GTDB r207" = "#D95F02", 
              "PlusPF" = "#7570B3")
colScale <- scale_colour_manual(name = "Database", values = myColors)

# Create a plot for each rank
plot_list <- list()

for(i in seq_along(rank_order)) {
    rank <- rank_order[i]
    rank_label <- rank_labels[i]
    
    rank_data <- plot_data %>% filter(taxonomic_rank == rank)
    
    if(nrow(rank_data) > 0 && length(unique(rank_data$sampleID)) > 1) {
        # Create plot for this rank
        p <- grouped_ggbetweenstats(rank_data,
                                  Database,
                                  value,
                                  ylab = paste("% reads assigned to", tolower(rank_label)),
                                  xlab = "Database",
                                  ggtheme = theme_bw(base_size = 14),
                                  pairwise.comparisons = TRUE,
                                  grouping.var = pretty_metric,
                                  sample.size.label = FALSE,
                                  ggsignif.args = list(textsize = 3, tip_length = 0.01),
                                  ggplot.component = list(colScale,
                                      theme(axis.title.y.right = element_blank(),
                                            axis.text.y.right = element_blank(),
                                            axis.ticks.y.right = element_blank())),
                                  results.subtitle = FALSE,
                                  centrality.plotting = FALSE) +
            labs(title = rank_label) +
            theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
        
        plot_list[[rank]] <- p
    } else {
        cat("Skipping", rank, "- insufficient data\n")
    }
}

# Combine all plots into a grid
if(length(plot_list) > 0) {
    # Arrange in a grid: 2 columns, multiple rows
    n_cols <- 2
    n_rows <- ceiling(length(plot_list) / n_cols)
    
    fig2_all_ranks <- wrap_plots(plot_list, ncol = n_cols, nrow = n_rows) +
        plot_annotation(tag_levels = 'A')
    
    # Save with appropriate size
    ggsave("manuscript_figures/fig2_all_ranks.png", fig2_all_ranks, 
           width = 16, height = 6 * n_rows, units = "in", dpi = 300)
    cat("\n✅ Saved figure with all ranks to: manuscript_figures/fig2_all_ranks.png\n")
    cat("   Total ranks shown:", length(plot_list), "\n")
    
    # Also create species-only figure (fig2.png)
    if("species" %in% names(plot_list)) {
        fig2_species <- plot_list[["species"]]
        ggsave("manuscript_figures/fig2.png", fig2_species, 
               width = 10, height = 6, units = "in", dpi = 300)
        cat("✅ Saved species-only figure to: manuscript_figures/fig2.png\n")
    } else {
        cat("⚠ Warning: Species plot not found, skipping fig2.png\n")
    }
} else {
    stop("No plots could be created - check data availability")
}

