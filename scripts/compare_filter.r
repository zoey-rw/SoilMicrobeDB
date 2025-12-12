# Compare read-level kmer filtering (on a subset of samples - files are quite large)

library(tidyverse)
library(ggrepel)
library(ggallin)
library(scales)
library(ggpmisc)
library(rstatix)
library(ggpubr)

seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds") %>% 
    select(-c(db_name, identified_reads))

filter_reads_dir = "data/classification/02_bracken_output"
filter_scores_list <- list.files(filter_reads_dir, pattern = "_scores.output", recursive = T, full.names = T)

# Create list of all samples from each database
filter_scores_list_pluspf <- filter_scores_list[grepl("pluspf", filter_scores_list)] 
filter_scores_list_smd <- filter_scores_list[grepl("soil_microbe_db", filter_scores_list)]
filter_scores_list_gtdb <- filter_scores_list[grepl("gtdb_207_scores", filter_scores_list)]

common1 = intersect(substr(basename(filter_scores_list_pluspf), 1, 19), 
                    substr(basename(filter_scores_list_smd), 1, 19))
common2 = intersect(common1, substr(basename(filter_scores_list_gtdb), 1, 19))

# Randomly select 10 samples that have been evaluated by all databases
set.seed(3)
common3 = common2 %>% sample(1)
common_files <- unique (grep(paste(common3,collapse="|"), 
												filter_scores_list, value=TRUE)) 
# Read in the (large) files
filter_scores_files <- common_files %>%
	setNames(., sub("_scores.output", "", basename(.))) %>%
	map(data.table::fread)
filter_scores1 = data.table::rbindlist(filter_scores_files, idcol = "samp_name", fill = T) 

filter_scores1 = filter_scores1 %>%  # default filter values used for Architeuthis
    mutate(pass_filter = ifelse(consistency > .9 &																	
                                    entropy < .1 & multiplicity < 2, 1, 0))

# Summarize by percent passing filter
pass_filter1 = filter_scores1 %>% 
	group_by(samp_name) %>% 
	add_tally(name = "n_classified_reads") %>% 
	group_by(samp_name, n_classified_reads, pass_filter) %>% 
	tally(name = "n_pass_filter") %>% 
	filter(pass_filter==1) %>% 
	mutate(percent_classified_passing = n_pass_filter / n_classified_reads)

# Parse sample IDs into grouping information
pass_filter2 <- pass_filter1 %>% 
	separate(samp_name, 
					 into = c("sampleID","db_name"), 
					 sep = "COMP_", remove = F, extra = "merge") %>% 
	mutate(sampleID = str_remove(samp_name, paste0("_",db_name)))
pass_filter2$siteID = substr(pass_filter2$sampleID, 1, 4)

pass_filter3 = left_join(pass_filter2, seq_depth_df) %>% 
	mutate(percent_classified = n_classified_reads / seq_depth,
	       percent_passing = n_pass_filter / seq_depth)
pass_filter = pass_filter3

# Save intermediate file
saveRDS(pass_filter,"data/classification/analysis_files/filter_summary_20samples.rds")


pass_filter = readRDS("data/classification/analysis_files/filter_summary_20samples.rds")

pass_filter_long = pass_filter %>% 
    pivot_longer(cols = c(percent_passing, percent_classified, percent_classified_passing), 
                                                names_to = "metric") %>% 
    mutate(pretty_metric = recode(metric, "percent_passing" = "% reads classified at\n high quality",
                                      "percent_classified" = "% reads classified",
                           "percent_classified_passing" = "% classified reads retained")) %>% 
    mutate(siteID = substr(sampleID, 1, 4))

pass_filter_long %>% group_by(metric, db_name) %>% 
    summarize(mean = mean(value), median = median(value))

#### Fit Data ####
fit <- pass_filter_long %>% 
    filter(metric != "percent_classified_passing") %>% 
    filter(seq_depth > 1000000) %>% 
    group_by(pretty_metric) %>% 
    anova_test(value ~ db_name) %>% 
    add_significance()

#### Run Tukey ###
tukey <- pass_filter_long %>% 
    filter(metric != "percent_classified_passing") %>% 
    filter(seq_depth > 1000000) %>% 
    group_by(pretty_metric) %>% 
    tukey_hsd(value ~ db_name) %>% 
    add_significance() %>% 
    add_xy_position(step.increase = .06)

ggplot(pass_filter_long %>% 
           filter(metric != "percent_classified_passing") %>% 
           filter(seq_depth > 1000000),
			 aes(y = value, x=reorder(db_name, value), color = pretty_metric)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
	geom_point(
		size=2, alpha=.4, show.legend = F,
		#position=position_jitter(width = .1, height=0.01)) +
		position=position_jitterdodge(jitter.height = .01, 
		                              jitter.width = .1, dodge.width = .75)) +
	theme_bw(base_size = 18) + 
	facet_wrap(~pretty_metric, drop=T, scales="free")  +
	ylab("% reads") + 
    xlab("Database name") +
	guides(color=guide_legend(NULL)) + 
    stat_pvalue_manual(tukey,
                       hide.ns = T) 


# Double-checking calculations: inspect values for one sample 
one_samp = filter_scores1 %>% filter(samp_name %in% c("HARV_001-O-20180710-COMP_gtdb_207",
                                                      "HARV_001-O-20180710-COMP_soil_microbe_db",
                                                      "HARV_001-O-20180710-COMP_pluspf"))


one_samp_pass = pass_filter3 %>% filter(samp_name %in% c("HARV_001-O-20180710-COMP_gtdb_207",
                                                         "HARV_001-O-20180710-COMP_soil_microbe_db",
                                                         "HARV_001-O-20180710-COMP_pluspf"))

ggplot(one_samp_pass) + 
    geom_histogram(aes(x = multiplicity, fill = samp_name)) + 
    facet_grid(rows=vars(samp_name)) + theme_bw()						 			 
ggplot(one_samp) + 
    geom_histogram(aes(x = multiplicity, fill = samp_name)) + 
    facet_grid(rows=vars(samp_name)) + theme_bw()							 			 


