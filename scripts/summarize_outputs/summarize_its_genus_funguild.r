# Summarize genus abundances for comparison with ITS estimates
library(tidyverse)
library(data.table)

# Read in ITS data from forecasting models
its_genus = fread("data/comparison_data/its_amplicon/NEON_ITS_amplicon.csv")

# Read in species-level abundances - not used
bracken_with_lineage=fread("data/classification/taxonomic_rank_summaries/species/soil_microbe_db_filtered_species_merged_lineage.csv", nThread = 8)

# Read in genus-level abundances
bracken_genus =fread("data/classification/taxonomic_rank_summaries/genus/soil_microbe_db_filtered_genus_merged_lineage.csv", nThread = 8)


phyla_fungi_summary <- read_csv("data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv") %>% 
    rename(fungal_abundance = phyla_fungi_metagenome,
           metagenome_sample_id = sampleID)

# Subset metagenome data to taxa in ITS 
bracken_genus_subset = bracken_genus %>% 
    mutate(name = tolower(name)) %>%  
    filter(name %in% colnames(its_genus)) 

unique_genera = bracken_genus %>% select(lineage, name) %>% distinct()
genus_fungi = unique_genera[grepl("Eukary", unique_genera$lineage),]$name

# Parse sample names into grouping information
bracken_genus_subset <- bracken_genus_subset %>% 
    separate(sample_id,  
             into = c("metagenome_sample_id","db_name"), sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(metagenome_sample_id = str_remove(sample_id, paste0("_",db_name))) %>% 
    cbind(microbialForecast:::parseNEONsampleIDs(.$metagenome_sample_id) %>% 
              select(siteID, plotID, dateID, horizon, plot_date))  %>% 
    #rename("metagenome_abundance" = "fraction_total_reads") %>% 
    select(-c(taxid_lineage, sample_id, kraken_assigned_reads, added_reads, 
              taxonomy_lvl, new_est_reads, db_name))

bracken_genus_subset = left_join(bracken_genus_subset, phyla_fungi_summary) %>% 
    mutate(metagenome_abundance = fraction_total_reads/fungal_abundance)

# Reshape data for plotting and add sampling info 
its_genus_long = its_genus %>% 
    cbind(microbialForecast:::parseNEONsampleIDs(.$ITS_sampleID) %>% 
              select(siteID, plotID, dateID, plot_date, horizon)) %>% 
    pivot_longer(cols = -c(siteID, plotID, dateID, ITS_sampleID, plot_date, horizon), 
                 names_to = "name", values_to = "ITS_abundance")
#its_genus_long$name =  str_to_title(its_genus_long$name)
its_genus_long$name = tolower(its_genus_long$name)

df_compare = left_join(bracken_genus_subset, its_genus_long, by = 
                           join_by(name, siteID, plotID, dateID, plot_date, horizon)) %>%
    mutate(year=as.factor(substr(dateID, 1, 4)),
           legacy = ifelse(year %in% c("2013","2014"), T, F)) 

unique(df_compare$name) # 268 overlapping genera

genus_fungi_ITS = unique(its_genus_long$name) %>% str_to_title

# Save!
write.csv(df_compare, "data/classification/analysis_files/genus_ITS_metagenome_comparison.csv")



# Query Funguild using taxon names
results_list_full = lapply(genus_fungi, function(taxon) 
    jsonlite::fromJSON(paste0("https://www.mycoportal.org/funguild/services/api/db_return.php?qDB=funguild_db&qField=taxon&qText=",taxon)))
funguild_results_full = do.call(rbind, results_list_full) %>% 
    as.data.frame() %>% select(taxon, trophicMode, guild, confidenceRanking) 


funguild_results_full = funguild_results_full %>%  mutate(functional_group = 
                                                              ifelse(grepl("Ecto", guild), 
                                                                     "Ectomycorrhizae", 
                                                                     ifelse(grepl("Arbusc", guild), 
                                                                            "Arbuscular mycorrhizae", 
                                                                            ifelse(grepl("Pathogen", guild), 
                                                                                   "Pathogen",
                                                                                   ifelse(grepl("Endophyte", guild), 
                                                                                          "Endophyte",
                                                                                          ifelse(grepl("Sapro", guild), 
                                                                                                 "Other saprotroph",
                                                                                                 NA))))))

funguild_results_full$name = tolower(funguild_results_full$taxon)

write_csv(funguild_results_full, "data/classification/analysis_files/FUNguild_assignments.csv")
