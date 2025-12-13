
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

# Load pre-processed filter results (extracted by 04_reshape_score_reads.r)
filter_results_file = "data/classification/analysis_files/filter_results_summary.csv"

if(!file.exists(filter_results_file)) {
    stop("❌ MISSING FILE: Filter results summary not found!\n",
         "   Expected: ", filter_results_file, "\n",
         "   Please run scripts/run_workflow/04_reshape_score_reads.r first to extract scoring information.")
}

cat("Loading filter results from:", filter_results_file, "\n")
summaries <- read_csv(filter_results_file, show_col_types = FALSE)
cat("Loaded", nrow(summaries), "records\n")

# Get list of samples that have been evaluated by all databases
sample_count = summaries %>% 
    distinct(samp_name, .keep_all = T) %>%  
    group_by(sampleID) %>% tally
common_samples = sample_count[sample_count$n > 2,]$sampleID
cat("Found", length(common_samples), "samples evaluated by all databases\n")

# Create prettier values for plotting
summaries = summaries %>% mutate(value = value * 100,
    Database = recode(db_name, 
                                                  "gtdb_207" = "GTDB r207",
                                                  "soil_microbe_db" = "SoilMicrobeDB",
                                                  "pluspf" = "PlusPF"),
    pretty_metric = recode(metric,"percent_classified" = "Read classification",
                           "percent_passing"="Read classification and quality-filtering")) %>% filter(metric %in% c("percent_classified","percent_passing"))

summaries$pretty_metric = factor(summaries$pretty_metric, levels = c("Read classification","Read classification and quality-filtering"))

# Plotting color palette
myColors <- brewer.pal(5,"Set1")
names(myColors) <- c("SoilMicrobeDB", "GTDB r207", "PlusPF")
colScale <- scale_colour_manual(name = "Database",values = myColors)

# Filter data for plotting
plot_data = summaries %>% 
    filter(!db_name %in% c("gtdb_207_unfiltered"))

# Use common_samples if available, otherwise use all samples
if(length(common_samples) > 0) {
    plot_data = plot_data %>% filter(sampleID %in% common_samples)
    cat("Plotting", nrow(plot_data), "records from", length(common_samples), "common samples\n")
} else {
    cat("No samples evaluated by all databases. Plotting all available samples.\n")
    cat("Plotting", nrow(plot_data), "records\n")
}

# Check if we have data to plot
if(nrow(plot_data) == 0) {
    stop("❌ ERROR: No data available for plotting!\n",
         "   Make sure filter_results_summary.csv contains data for the selected metrics.")
}

fig2 = grouped_ggbetweenstats(plot_data,
                              
                              #filter(!metric %in% c("seq_depth","n_reads")),
                              Database, 
                              #point.args = list(inherit_aes=T, aes(color=Database)),
                              value, 
                              ylab="% reads",
                              xlab="Database",
                               #type="np",
                              ggtheme = theme_bw(base_size = 16),
                              pairwise.comparisons = TRUE,
                              grouping.var = pretty_metric, 
                              sample.size.label = F,
                              ggsignif.args    = list(textsize = 4, tip_length = 0.01),
                              ggplot.component = list(colScale,
                                  theme(axis.title.y.right = element_blank(), 
                                                            axis.text.y.right = element_blank(), 
                                                            axis.ticks.y.right = element_blank())),
                              results.subtitle=F,
                              centrality.plotting = F) 

fig2 = fig2 & plot_annotation(tag_levels = 'A')

# Save figure
ggsave("manuscript_figures/fig2.png", fig2, width = 12, height = 6, units = "in", dpi = 300)

cat("✓ Figure saved to: manuscript_figures/fig2.png\n")

# Optional: Print summary statistics (commented out View() for non-interactive use)
# summaries %>% 
#     pivot_wider(values_from = "value", names_from = "metric") %>% 
#     distinct() %>% View

