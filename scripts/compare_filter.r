
library(tidyverse)
library(ggrepel)
library(ggallin)
library(scales)
library(ggpmisc)
library(rstatix)
library(ggpubr)
library(ggstatsplot)
library(future.apply)
library(patchwork)
library(RColorBrewer)

# Function to summarize the "scores" output from Architeuthis
# Format is one row per classified read
summarize_filter_scores = function(scores_df, seq_depth_df, max_entropy = .1,
                                   max_multiplicity = 2,
                                   min_consistency = .9) {
    
    scores_df = scores_df %>% mutate(db_name = ifelse(grepl("soil_microbe_db", samp_name), "soil_microbe_db",
                                                      ifelse(grepl("pluspf", samp_name), "pluspf",
                                                             ifelse(grepl("unfiltered", samp_name), "gtdb_207_unfiltered",
                                                                    ifelse(grepl("gtdb", samp_name), "gtdb_207_filtered", NA
                                                                    )))))
    
    # Add column for passing or not
    scores_df = scores_df %>% mutate(pass_filter =
                                         ifelse(multiplicity <= max_multiplicity &
                                                    consistency >= min_consistency &
                                                    entropy <= max_entropy, 1, 0))
    
    # Summarize per-read scores
    score_summary_df = scores_df  %>%
        filter(n_kmers > 0) %>%
        group_by(samp_name) %>%
        reframe(mean_consistency = mean(consistency, na.rm=T),
                mean_multiplicity = mean(multiplicity, na.rm=T),
                mean_entropy = mean(entropy, na.rm=T),
                mean_confidence = mean(confidence,  na.rm=T)) %>% distinct
    
    # Percent passing or not
    filter_summary_df = scores_df %>%
        # filter(n_kmers > 0) %>%
        group_by(db_name, samp_name) %>%
        add_tally(name = "n_total_classified_reads") %>%
        group_by(db_name, samp_name, pass_filter) %>%
        add_tally(name = "n_reads") %>%
        group_by(db_name, samp_name, pass_filter,n_total_classified_reads,n_reads) %>%
        reframe(pct_of_classified_passing = n_reads/n_total_classified_reads) %>% 
        filter(pass_filter == 1) %>% ungroup %>%
        select(-pass_filter) %>% distinct() %>% 
        
        # Parse sample IDs into grouping information
        separate(samp_name,
                 into = c("sampleID","db_name"),
                 sep = "COMP_", remove = F, extra = "merge") %>%
        mutate(sampleID = str_remove(samp_name, paste0("_",db_name))) %>% distinct
    
    # Include sequencing depth to get overall percent
    filter_summary_df = left_join(filter_summary_df, seq_depth_df, by=join_by(sampleID)) %>%
        mutate(percent_classified = n_total_classified_reads / seq_depth,
               percent_passing = n_reads / seq_depth)
    
    # Combine and reshape for output
    out_df = full_join(score_summary_df, filter_summary_df, by=join_by(samp_name)) %>% 
        pivot_longer(cols=-c(db_name, sampleID, samp_name), names_to = "metric") %>% 
        distinct()
    
    out_df
}


seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds") %>% 
    select(-c(db_name, identified_reads))

filter_reads_dir = "data/classification/02_bracken_output"
if(!dir.exists(filter_reads_dir)) {
    stop("❌ MISSING DIRECTORY: Filter scores directory not found!\n",
         "   Expected: ", filter_reads_dir, "\n",
         "   This directory should contain Architeuthis filter score files (*_scores.output).\n",
         "   Please add this directory and the score files to continue.")
}
filter_scores_list <- list.files(filter_reads_dir, 
                                 pattern = "_scores.output", recursive = T, full.names = T) 
if(length(filter_scores_list) == 0) {
    stop("❌ MISSING FILES: No filter score files found!\n",
         "   Expected files matching pattern '*_scores.output' in: ", filter_reads_dir, "\n",
         "   These are output files from Architeuthis filter analysis.\n",
         "   Please add the score files to continue.")
}
filter_scores_list <- filter_scores_list[!grepl("gtdb_207_unfiltered_scores",filter_scores_list)]
# Create list of all samples from each database
filter_scores_list_pluspf <- filter_scores_list[grepl("pluspf_scores", filter_scores_list)]
filter_scores_list_smd <- filter_scores_list[grepl("soil_microbe_db_scores", filter_scores_list)]
filter_scores_list_gtdb <- filter_scores_list[grepl("gtdb_207_scores", filter_scores_list)]
all_filter_scores_list = c(filter_scores_list_pluspf, filter_scores_list_smd, 
                           filter_scores_list_gtdb)

# Names of samples that have been evaluated by all databases
common_samples = substr(basename(all_filter_scores_list), 1, 19) %>%
    table %>% sort %>% stack %>%
    filter(values > 2) %>% select(ind) %>% unlist
# Randomly select 10 samples that have been evaluated by all databases
set.seed(1)
common_samples
common5 = common_samples
common5 = sample(common_samples, 50, replace = F)
common_files <- unique (grep(paste(common5,collapse="|"),
												filter_scores_list, value=TRUE)) 

common_files <- common_files[!grepl("unfiltered", common_files)]


# Run summary function over every input file (returns one row per sample x database)
plan(multisession, workers=18)
Sys.time()
file_names = common_files
summaries_list = future_lapply(file_names, function(x) {
# Read in the (large) files
    df_in = data.table::fread(x)
    df_in$samp_name = sub("_scores.output", "", basename(x))
    out <- summarize_filter_scores(df_in, seq_depth_df)
    return(out) }
)
Sys.time() # takes 5-10 minutes

summaries = data.table::rbindlist(summaries_list)


write_csv(summaries, "data/classification/analysis_files/filter_results_summary.csv")


summaries <- read_csv("data/classification/analysis_files/filter_results_summary.csv")

# Get list of samples that have been evaluated by all databases
# sample_count = summaries %>% 
#     distinct(samp_name, .keep_all = T) %>%  
#     group_by(sampleID) %>% tally
#common_samples = sample_count[sample_count$n>2,]$sampleID

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

fig2 = grouped_ggbetweenstats(summaries %>% 
                                  filter(sampleID %in% common_samples) %>% 
                                  filter(!db_name %in% c("gtdb_207_unfiltered")),
                              
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


summaries %>% 
    pivot_wider(values_from = "value", names_from = "metric") %>% 
    distinct() %>% View

