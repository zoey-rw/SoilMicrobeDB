# Calculate the % classified at the species level

library(tidyverse)
library(data.table)
library(ggpubr)
library(ggpmisc)
library(rstatix)
source("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/helper_functions.r")
source("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/comets_shinyapp_example/source.R")


#Read in bracken estimates from 3 databases
species_bracken1 = fread("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_species_merged.csv")
species_bracken2 = fread("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/gtdb_species_merged.csv")
species_bracken3 = fread("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/pluspf_species_merged.csv")


species_bracken = rbindlist(list(species_bracken1, species_bracken2, species_bracken3))


# Summarize by % passing filter
filter_species = species_bracken %>% group_by(sample_id) %>% 
    mutate(n_classified_reads = sum(kraken_assigned_reads)) %>% 
    group_by(sample_id, n_classified_reads) %>% 
    summarize(n_pass_filter = sum(new_est_reads)) %>% 
    separate(sample_id,  
             into = c("sampleID","db_name"), sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(sample_id, paste0("_",db_name))) %>% 
    mutate(db_name = gsub("_filtered","",db_name)) 



# Add in sequencing depth
seq_depth_df <- readRDS("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/seq_depth_df.rds") 

pass_filter_species = left_join(filter_species, #seq_depth_df) %>% 
                                seq_depth_df %>% 
                                    select(-c(db_name, sample_id, identified_reads))) %>% 
    mutate(percent_classified = n_classified_reads / seq_depth,
           percent_passing = n_pass_filter / seq_depth)


# Reshape for plotting
pass_filter_species_long = pass_filter_species %>% 
    select(sampleID,seq_depth,db_name, percent_passing, percent_classified) %>% 
    pivot_longer(cols = c(percent_passing, percent_classified), 
                 names_to = "metric") %>% 
    mutate(pretty_metric = recode(metric, 
                                  "percent_passing" = "% reads classified at\n high quality",
                                  "percent_classified" = "% reads classified")) %>% 
    mutate(siteID = substr(sampleID, 1, 4))

write.csv(pass_filter_species_long, "/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/pass_filter_summary.csv")
