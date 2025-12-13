library(neonUtilities)
library(tidyverse)
library(broom)
library(ggpubr)
library(SimplyAgree)

source("scripts/helper_functions.r")
# Note: source.R from comets_shinyapp_example may not be available locally
# source("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/comets_shinyapp_example/source.R")

# Read in & combine the metagenome data
domain_bracken <- readRDS("data/classification/taxonomic_rank_summaries/domain/bracken_domain_estimates.rds") 
domain_metagenome = domain_bracken %>% 
    separate(samp_name,  
             into = c("compositeSampleID","db_name"), sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(compositeSampleID = str_remove(samp_name, paste0("_",db_name))) %>% 
    select(compositeSampleID, db_name, taxon, percentage) %>% 
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
metagenome = left_join(metagenome, phyla_fungi_summary)

seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds") %>% 
    select(-c(db_name, identified_reads)) %>% 
    rename(compositeSampleID = sampleID)  # sampleID column actually contains compositeSampleIDs
metagenome = left_join(metagenome, seq_depth_df) %>% 
    filter(seq_depth > 1000000)

# Subset soil data (from NEON data RDS file)
# Expected file: data/comparison_data/neon/neon_soil_data_2023.rds
# Or load via: soilData <- readRDS("data/comparison_data/neon/neon_soil_data_2023.rds"); soilCores <- soilData$sls_soilCoreCollection
if(exists("soilCores")) {
soil_master_df = soilCores %>% select(compositeSampleID, siteID,sampleID,
                                      plotID, biomassID,geneticSampleID, biome,
                                      nlcdClass, horizon, sampleBottomDepth,sampleTopDepth, 
                                      litterDepth, standingWaterDepth, soilTemp, 
                                      sampleTiming) %>% 
   # distinct(sampleID, .keep_all = T) %>% 
    filter(!is.na(biomassID))
} else {
    warning("⚠️  MISSING FILE: soilCores not found!")
    warning("   Expected: data/comparison_data/neon/neon_soil_data_2023.rds")
    warning("   Load with: soilData <- readRDS('data/comparison_data/neon/neon_soil_data_2023.rds'); soilCores <- soilData$sls_soilCoreCollection")
    warning("   Continuing with minimal structure - some metadata will be missing")
    # If soilCores is not available, create minimal structure
    soil_master_df = metagenome %>% select(compositeSampleID) %>% distinct()
}

# Read in the processed PLFA data 
plfa_file <- "data/comparison_data/plfa/NEON_microbial_biomass_PLFA.rds"
if(file.exists(plfa_file)) {
    plfaData <- readRDS(plfa_file)
} else {
    stop("❌ MISSING FILE: PLFA data not found!\n",
         "   Expected: ", plfa_file, "\n",
         "   This file contains processed PLFA (phospholipid fatty acid) data from NEON.\n",
         "   Please add this file to continue.")
}


master_df = left_join(soil_master_df, plfaData) 
master_df = left_join(master_df, metagenome)

write.csv(master_df, "data/classification/analysis_files/plfa_comparison.csv")
