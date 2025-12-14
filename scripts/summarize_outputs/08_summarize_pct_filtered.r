# Calculate the percentage of reads classified at species level
# Input: Species-level merged CSV files from three databases:
#        soil_microbe_db_species_merged.csv, gtdb_207_species_merged.csv, pluspf_species_merged.csv
#        Sequencing depth from seq_depth_df.rds
# Output: pass_filter_summary.csv with percent_classified and percent_passing metrics
#         percent_classified: % of total reads classified by Kraken2
#         percent_passing: % of reads passing Architeuthis quality filter

library(tidyverse)
library(data.table)

species_bracken <- rbindlist(list(
    fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_species_merged.csv"),
    fread("data/classification/taxonomic_rank_summaries/gtdb_207_species_merged.csv"),
    fread("data/classification/taxonomic_rank_summaries/pluspf_species_merged.csv")
))

filter_species <- species_bracken %>% 
    group_by(sample_id) %>% 
    mutate(n_classified_reads = sum(kraken_assigned_reads)) %>% 
    group_by(sample_id, n_classified_reads) %>% 
    summarize(n_pass_filter = sum(new_est_reads), .groups = "drop") %>% 
    separate(sample_id, into = c("sampleID", "db_name"), sep = "COMP_", remove = FALSE, extra = "merge") %>% 
    mutate(sampleID = str_remove(sample_id, paste0("_", db_name)),
           db_name = gsub("_filtered", "", db_name))

pass_filter_species <- left_join(filter_species, 
                                readRDS("data/classification/analysis_files/seq_depth_df.rds") %>% 
                                    select(-c(db_name, sample_id, identified_reads))) %>%
    mutate(percent_classified = n_classified_reads / seq_depth,
           percent_passing = n_pass_filter / seq_depth)

pass_filter_species_long <- pass_filter_species %>% 
    select(sampleID, seq_depth, db_name, percent_passing, percent_classified) %>% 
    pivot_longer(cols = c(percent_passing, percent_classified), names_to = "metric") %>% 
    mutate(pretty_metric = recode(metric, 
                                  "percent_passing" = "% reads classified at\n high quality",
                                  "percent_classified" = "% reads classified"),
           siteID = substr(sampleID, 1, 4))

write.csv(pass_filter_species_long, "data/classification/analysis_files/pass_filter_summary.csv", row.names = FALSE)
