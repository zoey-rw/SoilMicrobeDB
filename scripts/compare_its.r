# Compare metagenomic abundances with those from ITS rRNA sequencing

library(tidyverse)
library(data.table)
library(ggpubr)

# Read in ITS data from forecasting models
abundance_ITS = readRDS("/projectnb/dietzelab/zrwerbin/microbialForecasts/data/clean/groupAbundances_ITS_2023.rds")
its_genus = abundance_ITS$genus_fun

# Read in species-level abundances - not used
bracken_with_lineage=fread("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_species_merged_lineage.csv", nThread = 8)

# Read in genus-level abundances
bracken_genus =fread("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_genus_merged_lineage.csv", nThread = 8)

# Subset metagenome data to taxa in ITS 
bracken_genus_subset = bracken_genus %>% 
    mutate(name = tolower(name)) %>%  
    filter(name %in% colnames(its_genus)) 

# Parse sample names into grouping information
bracken_genus_subset <- bracken_genus_subset %>% 
    separate(sample_id,  
             into = c("metagenome_sample_id","db_name"), sep = "COMP_", remove = F, extra = "merge") %>% 
     mutate(metagenome_sample_id = str_remove(sample_id, paste0("_",db_name))) %>% 
    cbind(microbialForecast:::parseNEONsampleIDs(.$metagenome_sample_id) %>% 
              select(siteID, plotID, dateID, horizon, plot_date))  %>% 
    rename("metagenome_abundance" = "fraction_total_reads") %>% 
    select(-c(taxid_lineage, sample_id, kraken_assigned_reads, added_reads, 
              taxonomy_lvl, new_est_reads, db_name))

# Reshape data for plotting and add sampling horizon 
its_genus_long = its_genus %>% 
    pivot_longer(cols = -c(siteID, plotID, dateID, sampleID, dates, plot_date), names_to = "name", values_to = "ITS_abundance") %>% 
    cbind(microbialForecast:::parseNEONsampleIDs(.$sampleID) %>% 
              select(horizon))
its_genus_long <- its_genus_long %>% rename("ITS_sample_id" = "sampleID")


df_compare = left_join(bracken_genus_subset, its_genus_long, by = 
                           join_by(name, siteID, plotID, dateID, plot_date, horizon))

unique(df_compare$name) # 17 overlapping genera

# Save!
write.csv(df_compare, "/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/genus_ITS_metagenome_comparison.csv")

# Visualize a handful of taxa for manuscript
ggplot(df_compare %>% filter(name %in% c("cenococcum","glomus","russula","umbelopsis")), 
       aes(x = metagenome_abundance, y = ITS_abundance)) + 
    geom_point(alpha=.5)  +
    stat_smooth(method="lm") +
    scale_y_sqrt() + 
    scale_x_sqrt() + 
    geom_abline(slope=1, intercept = 0, linetype=2) + 
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=6, label.x.npc = .4) + 
    facet_grid(~name, scales="free") + 
    theme_bw(base_size = 18)  +
    xlab("Relative abundance in metagenome") +
    ylab("Relative abundance in ITS amplicon")  + 
    theme(axis.text.x = element_text(angle = 270))


ggplot(df_compare,# %>% filter(name %in% c("cenococcum","glomus","russula","umbelopsis")), 
       aes(x = metagenome_abundance, y = ITS_abundance)) + 
    geom_point(alpha=.5)  +
    stat_smooth(method="lm") +
    scale_y_sqrt() + 
    scale_x_sqrt() + 
    geom_abline(slope=1, intercept = 0, linetype=2) + 
    stat_regline_equation(aes(label = ..rr.label..), color=2,
                          show.legend = FALSE, size=7, label.x.npc = .4) + 
    facet_wrap(~name, scales="free") + 
    theme_bw(base_size = 18)  +
    xlab("Relative abundance in metagenome") +
    ylab("Relative abundance in ITS amplicon")  + 
    theme(axis.text.x = element_text(angle = 270))


