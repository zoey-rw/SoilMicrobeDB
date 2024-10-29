# Compare metagenomic abundances with those from ITS rRNA sequencing

library(tidyverse)
library(data.table)
library(ggpubr)

# Read in ITS data from forecasting models
abundance_ITS = readRDS("/projectnb/dietzelab/zrwerbin/microbialForecasts/data/clean/groupAbundances_ITS_2023.rds")
its_genus = abundance_ITS$genus_fun

# Read in species-level abundances
bracken_with_lineage=fread("/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_species_merged_lineage.csv", nThread = 8)

# Read in genus-level abundances
bracken_genus =fread("/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_genus_merged_lineage.csv", nThread = 8)
bracken_genus_subset = bracken_genus %>% 
    mutate(name = tolower(name)) %>%  
    filter(name %in% colnames(its_genus)) 

# Parse sample names into grouping information
bracken_genus_subset <- bracken_genus_subset %>% 
    separate(sample_id,  
             into = c("compositeSampleID","db_name"), sep = "COMP_", remove = F, extra = "merge") %>% 
     mutate(compositeSampleID = str_remove(sample_id, paste0("_",db_name))) %>% 
    cbind(microbialForecast:::parseNEONsampleIDs(bracken_genus_subset$compositeSampleID) %>% 
              select(siteID, plotID, dateID, horizon, plot_date))

its_genus_long = its_genus %>% 
    pivot_longer(cols = -c(siteID, plotID, dateID, sampleID, dates, plot_date), names_to = "name", values_to = "ITS_abundance") %>% 
    cbind(microbialForecast:::parseNEONsampleIDs(.$sampleID) %>% select(horizon))

df_compare = left_join(bracken_genus_subset, its_genus_long, by = join_by(name, siteID, plotID, dateID, plot_date, horizon))

# Visualize a handful of taxa
ggplot(df_compare %>% filter(name %in% c("cenococcum","glomus","russula","umbelopsis")), 
       aes(x = fraction_total_reads, y = ITS_abundance)) + 
    geom_point(alpha=.5)  +
    stat_smooth(method="lm") +
    scale_y_sqrt() + 
    scale_x_sqrt() + 
    geom_abline(slope=1, intercept = 0, linetype=2) + 
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=6, label.x.npc = .4) + 
    facet_grid(~name, scales="free") + 
    theme_bw(base_size = 18)  +
    xlab("Abundance in metagenome") +
    ylab("Abundance in ITS amplicon")  + theme(axis.text.x = element_text(angle = 270))


