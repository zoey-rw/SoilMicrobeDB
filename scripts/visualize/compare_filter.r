library(tidyverse)
library(ggpubr)
library(ggstatsplot)
library(patchwork)
library(RColorBrewer)

summaries <- read_csv("data/classification/analysis_files/filter_results_summary.csv", show_col_types = FALSE)

valid_dbs <- c("gtdb_207", "gtdb_207_unfiltered", "pluspf", "soil_microbe_db")

if("samp_name" %in% names(summaries) && "db_name" %in% names(summaries)) {
    summaries <- summaries %>%
        mutate(
            db_name = case_when(
                db_name %in% valid_dbs ~ db_name,
            grepl("gtdb_207_unfiltered", samp_name) ~ "gtdb_207_unfiltered",
            grepl("gtdb_207", samp_name) & !grepl("gtdb_207_unfiltered", samp_name) ~ "gtdb_207",
                grepl("_pluspf$", samp_name) ~ "pluspf",
                grepl("_soil_microbe_db$", samp_name) ~ "soil_microbe_db",
                grepl("_gtdb_207_unfiltered$", samp_name) ~ "gtdb_207_unfiltered",
                grepl("_gtdb_207$", samp_name) & !grepl("_gtdb_207_unfiltered", samp_name) ~ "gtdb_207",
                db_name == "gtdb_207_filtered" ~ "gtdb_207",
            TRUE ~ db_name
            )
        )
}

plot_data <- summaries %>%
    filter(metric %in% c("percent_classified", "percent_passing")) %>%
    mutate(value = value * 100,
           db_name = ifelse(db_name == "gtdb_207_filtered", "gtdb_207", db_name),
           Database = recode(db_name,
                            "gtdb_207" = "GTDB r207 (filtered)",
                            "gtdb_207_unfiltered" = "GTDB r207 (unfiltered)",
                            "soil_microbe_db" = "SoilMicrobeDB",
                            "pluspf" = "PlusPF"),
           pretty_metric = recode(metric,
                                 "percent_classified" = "Read classification",
                                 "percent_passing" = "Read classification and quality-filtering"),
           pretty_metric = factor(pretty_metric, levels = c("Read classification", "Read classification and quality-filtering"))) %>%
    filter(db_name %in% c("gtdb_207", "gtdb_207_unfiltered", "pluspf", "soil_microbe_db"))

pluspf_samples <- plot_data %>% filter(db_name == "pluspf") %>% distinct(sampleID) %>% pull(sampleID)
smdb_samples <- plot_data %>% filter(db_name == "soil_microbe_db") %>% distinct(sampleID) %>% pull(sampleID)
gtdb_unfiltered_samples <- plot_data %>% filter(db_name == "gtdb_207_unfiltered") %>% distinct(sampleID) %>% pull(sampleID)
gtdb_filtered_samples <- plot_data %>% filter(db_name == "gtdb_207") %>% distinct(sampleID) %>% pull(sampleID)

common_pluspf_smdb <- intersect(pluspf_samples, smdb_samples)
common_all_unfiltered <- intersect(intersect(pluspf_samples, smdb_samples), gtdb_unfiltered_samples)
common_all_filtered <- intersect(intersect(pluspf_samples, smdb_samples), gtdb_filtered_samples)

# Use unfiltered GTDB for this comparison (only unfiltered should be in the plot)
if(length(common_all_unfiltered) > 0) {
    gtdb_db_to_use <- "gtdb_207_unfiltered"
    common_samples <- common_all_unfiltered
    cat("Using GTDB unfiltered with", length(common_samples), "common samples\n")
} else {
    stop("No common samples found across all databases (PlusPF, SoilMicrobeDB, and GTDB unfiltered)")
}

# Only include unfiltered GTDB in the plot (explicitly exclude filtered)
plot_data <- plot_data %>%
    filter(db_name %in% c("gtdb_207_unfiltered", "pluspf", "soil_microbe_db")) %>%  # Only unfiltered, explicitly exclude gtdb_207
    filter(sampleID %in% common_samples) %>%
    mutate(Database = factor(Database, levels = c("SoilMicrobeDB", 
                                                   "GTDB r207 (unfiltered)",
                                                   "PlusPF")))

myColors <- brewer.pal(5, "Set1")[1:3]
names(myColors) <- c("SoilMicrobeDB", "GTDB r207 (unfiltered)", "PlusPF")
colScale <- scale_colour_manual(name = "Database", values = myColors)

fig2 <- grouped_ggbetweenstats(plot_data,
                              Database,
                              value,
                              ylab = "% reads",
                              xlab = "Database",
                              ggtheme = theme_bw(base_size = 16),
                              pairwise.comparisons = TRUE,
                              grouping.var = pretty_metric,
                              sample.size.label = FALSE,
                              ggsignif.args = list(textsize = 4, tip_length = 0.01),
                              ggplot.component = list(colScale,
                                  theme(axis.title.y.right = element_blank(),
                                        axis.text.y.right = element_blank(),
                                        axis.ticks.y.right = element_blank())),
                              results.subtitle = FALSE,
                              centrality.plotting = FALSE) &
    plot_annotation(tag_levels = 'A')

ggsave("manuscript_figures/fig2.png", fig2, width = 12, height = 6, units = "in", dpi = 300)

