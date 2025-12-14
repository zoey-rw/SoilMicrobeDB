# Merge NEON PLFA (phospholipid fatty acid) data with metagenome classification results
# Input: Domain-level Bracken estimates, fungal summaries (genus and phylum)
#        Sequencing depth, PLFA data from data/comparison_data/plfa/NEON_microbial_biomass_PLFA.rds
#        Soil core metadata from helper_functions.r
# Output: plfa_comparison.csv with PLFA and metagenome data merged
#         Filters to samples with seq_depth > 1,000,000 and non-NA biomassID

library(tidyverse)
source("scripts/helper_functions.r")

domain_bracken <- readRDS("data/classification/taxonomic_rank_summaries/domain/bracken_domain_estimates.rds")

domain_metagenome <- domain_bracken %>% 
    separate(samp_name, into = c("compositeSampleID", "db_name"), sep = "COMP_", remove = FALSE, extra = "merge") %>% 
    mutate(compositeSampleID = str_remove(samp_name, paste0("_", db_name))) %>% 
    select(compositeSampleID, db_name, taxon, percentage) %>% 
    pivot_wider(names_from = "taxon", values_from = "percentage") %>% 
    ungroup() %>% 
    select(-c(Viruses, "Classified at a higher level", Unclassified)) %>% 
    mutate(fungi_domain_metagenome = Eukaryota / (Eukaryota + Bacteria)) %>% 
    filter(db_name == "soil_microbe_db")

phyla_fungi_summary <- read_csv("data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv")
genus_fungi_summary <- read_csv("data/classification/analysis_files/soil_microbe_db_genus_fungi_summary.csv") %>% 
    mutate(genus_fungi_minus_amf_metagenome = genus_fungi_metagenome - genus_amf_metagenome)

metagenome <- left_join(domain_metagenome, genus_fungi_summary) %>%
    left_join(phyla_fungi_summary) %>%
    left_join(readRDS("data/classification/analysis_files/seq_depth_df.rds") %>% 
                 select(-c(db_name, identified_reads)) %>% 
                 rename(compositeSampleID = sampleID)) %>%
    filter(seq_depth > 1000000)

source("scripts/helper_functions.r")
if(!exists("soilCores")) soilCores <- load_soilCores()

soil_master_df <- soilCores %>% 
    select(compositeSampleID, siteID, sampleID, plotID, biomassID, geneticSampleID, biome,
           nlcdClass, horizon, sampleBottomDepth, sampleTopDepth, 
           litterDepth, standingWaterDepth, soilTemp, sampleTiming) %>% 
    filter(!is.na(biomassID))

plfa_file <- "data/comparison_data/plfa/NEON_microbial_biomass_PLFA.rds"
if(!file.exists(plfa_file)) {
    stop("PLFA data not found: ", plfa_file)
}

master_df <- left_join(soil_master_df, readRDS(plfa_file)) %>% 
    left_join(metagenome)

write.csv(master_df, "data/classification/analysis_files/plfa_comparison.csv", row.names = FALSE)
