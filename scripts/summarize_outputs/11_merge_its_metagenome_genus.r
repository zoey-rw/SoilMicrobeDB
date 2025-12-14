# Merge ITS amplicon data with metagenome genus-level abundances
# Input: ITS amplicon genus abundances from 03_summarize_its_amplicon.r
#        Metagenome genus abundances from soil_microbe_db_genus_merged_lineage.csv
#        Phylum-level fungal summary for normalization
# Output: genus_ITS_metagenome_comparison.csv comparing ITS and metagenome abundances
#         Matches by genus name, siteID, plotID, dateID, horizon
#         Includes legacy flag for pre-2015 samples

library(tidyverse)
library(data.table)

its_genus <- fread("data/comparison_data/its_amplicon/NEON_ITS_amplicon.csv")
bracken_genus <- fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_genus_merged_lineage.csv", nThread = 8)

phyla_fungi_summary <- read_csv("data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv") %>% 
    rename(fungal_abundance = phyla_fungi_metagenome, metagenome_sample_id = sampleID)

bracken_genus_subset <- bracken_genus %>% 
    mutate(name = tolower(name)) %>%  
    filter(name %in% colnames(its_genus)) %>%
    separate(sample_id, into = c("metagenome_sample_id", "db_name"), sep = "COMP_", remove = FALSE, extra = "merge") %>% 
    mutate(metagenome_sample_id = str_remove(sample_id, paste0("_", db_name))) %>% 
    cbind(microbialForecast:::parseNEONsampleIDs(.$metagenome_sample_id) %>% 
              select(siteID, plotID, dateID, horizon, plot_date)) %>%
    select(-c(taxid_lineage, sample_id, kraken_assigned_reads, added_reads, 
              taxonomy_lvl, new_est_reads, db_name)) %>%
    left_join(phyla_fungi_summary) %>% 
    mutate(metagenome_abundance = fraction_total_reads / fungal_abundance)

its_genus_long <- its_genus %>% 
    cbind(microbialForecast:::parseNEONsampleIDs(.$ITS_sampleID) %>% 
              select(siteID, plotID, dateID, plot_date, horizon)) %>% 
    pivot_longer(cols = -c(siteID, plotID, dateID, ITS_sampleID, plot_date, horizon), 
                 names_to = "name", values_to = "ITS_abundance") %>%
    mutate(name = tolower(name))

df_compare <- left_join(bracken_genus_subset, its_genus_long, 
                       by = join_by(name, siteID, plotID, dateID, plot_date, horizon)) %>%
    mutate(year = as.factor(substr(dateID, 1, 4)),
           legacy = year %in% c("2013", "2014"))

write.csv(df_compare, "data/classification/analysis_files/genus_ITS_metagenome_comparison.csv", row.names = FALSE)
