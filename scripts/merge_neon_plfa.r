library(neonUtilities)
library(tidyverse)
library(broom)
library(ggpubr)
library(SimplyAgree)

source("scripts/helper_functions.r")

# Read in & combine the metagenome data
domain_bracken <- readRDS("data/classification/taxonomic_rank_summaries/domain/bracken_domain_estimates.rds") 
domain_metagenome = domain_bracken %>% 
    separate(samp_name,  
             into = c("sampleID","db_name"), sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(samp_name, paste0("_",db_name))) %>% 
    select(sampleID, db_name, taxon, percentage) %>% 
    pivot_wider(names_from = "taxon", values_from = "percentage") %>% 
    ungroup() %>% #select(-sample_id) %>% 
    select(-c(Viruses,"Classified at a higher level", Unclassified)) %>% 
    mutate(fungi_domain_metagenome = Eukaryota/(Eukaryota+Bacteria)) %>% 
    filter(db_name=="soil_microbe_db")
phyla_fungi_summary <- read_csv("data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv")
genus_fungi_summary <- read_csv("data/classification/analysis_files/soil_microbe_db_genus_fungi_summary.csv") %>% 
    mutate(genus_fungi_minus_amf_metagenome = 
               genus_fungi_metagenome - genus_amf_metagenome)
metagenome = left_join(domain_metagenome, genus_fungi_summary)
metagenome = left_join(metagenome, phyla_fungi_summary) %>% rename(compositeSampleID=sampleID)

seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds") %>% 
    select(-c(db_name, identified_reads)) %>% rename(compositeSampleID=sampleID)
metagenome = left_join(metagenome, seq_depth_df) %>% 
    filter(seq_depth > 1000000)

# Subset soil data (from source.r)
soil_master_df = soilCores %>% select(compositeSampleID, siteID,sampleID,
                                      plotID, biomassID,geneticSampleID, biome,
                                      nlcdClass, horizon, sampleBottomDepth,sampleTopDepth, 
                                      litterDepth, standingWaterDepth, soilTemp, 
                                      sampleTiming) %>% 
   # distinct(sampleID, .keep_all = T) %>% 
    filter(!is.na(biomassID))

# Read in the processed PLFA data 
plfaData <- readRDS("/projectnb/talbot-lab-data/zrwerbin/SoilBiomassNEON/NEON_microbial_biomass_PLFA.rds")


master_df = left_join(soil_master_df, plfaData) 
master_df = left_join(master_df, metagenome)

write.csv(master_df, "data/classification/analysis_files/plfa_comparison.csv")
