# Merge ITS amplicon data with metagenome genus-level abundances
# Compares ITS amplicon sequencing with metagenome classification at genus level

library(tidyverse)
library(data.table)

# Read in ITS data from forecasting models
its_genus <- fread("data/comparison_data/its_amplicon/NEON_ITS_amplicon.csv")

# Read in genus-level abundances
bracken_genus <- fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_genus_merged_lineage.csv", nThread = 8)

phyla_fungi_summary <- read_csv("data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv") %>% 
    rename(fungal_abundance = phyla_fungi_metagenome,
           metagenome_sample_id = sampleID)

# Subset metagenome data to taxa in ITS 
bracken_genus_subset <- bracken_genus %>% 
    mutate(name = tolower(name)) %>%  
    filter(name %in% colnames(its_genus)) 

# Parse sample names into grouping information
bracken_genus_subset <- bracken_genus_subset %>% 
    separate(sample_id,  
             into = c("metagenome_sample_id","db_name"), sep = "COMP_", remove = FALSE, extra = "merge") %>% 
    mutate(metagenome_sample_id = str_remove(sample_id, paste0("_",db_name))) %>% 
    cbind(microbialForecast:::parseNEONsampleIDs(.$metagenome_sample_id) %>% 
              select(siteID, plotID, dateID, horizon, plot_date))  %>% 
    select(-c(taxid_lineage, sample_id, kraken_assigned_reads, added_reads, 
              taxonomy_lvl, new_est_reads, db_name))

bracken_genus_subset <- left_join(bracken_genus_subset, phyla_fungi_summary) %>% 
    mutate(metagenome_abundance = fraction_total_reads/fungal_abundance)

# Reshape data for plotting and add sampling info 
its_genus_long <- its_genus %>% 
    cbind(microbialForecast:::parseNEONsampleIDs(.$ITS_sampleID) %>% 
              select(siteID, plotID, dateID, plot_date, horizon)) %>% 
    pivot_longer(cols = -c(siteID, plotID, dateID, ITS_sampleID, plot_date, horizon), 
                 names_to = "name", values_to = "ITS_abundance")
its_genus_long$name <- tolower(its_genus_long$name)

df_compare <- left_join(bracken_genus_subset, its_genus_long, by = 
                           join_by(name, siteID, plotID, dateID, plot_date, horizon)) %>%
    mutate(year = as.factor(substr(dateID, 1, 4)),
           legacy = ifelse(year %in% c("2013","2014"), TRUE, FALSE)) 

cat("✓ Found", length(unique(df_compare$name)), "overlapping genera\n")

# Save comparison file
write.csv(df_compare, "data/classification/analysis_files/genus_ITS_metagenome_comparison.csv", row.names = FALSE)
cat("✅ Saved ITS-metagenome comparison to: data/classification/analysis_files/genus_ITS_metagenome_comparison.csv\n")
