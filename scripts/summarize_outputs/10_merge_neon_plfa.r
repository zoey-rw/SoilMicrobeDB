# Merge NEON PLFA (phospholipid fatty acid) data with metagenome classification results
# Input: Domain-level Bracken estimates, fungal summaries (genus and phylum)
#        Sequencing depth, PLFA data from data/comparison_data/plfa/NEON_microbial_biomass_PLFA.rds
#        Soil core metadata from helper_functions.r
# Output: plfa_comparison.csv with PLFA and metagenome data merged
#         Filters to samples with seq_depth > 1,000,000 and non-NA biomassID

library(tidyverse)
source("scripts/helper_functions.r")

required_files <- c(
    "data/classification/taxonomic_rank_summaries/soil_microbe_db_domain_merged_lineage.csv",
    "data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv",
    "data/classification/analysis_files/soil_microbe_db_genus_fungi_summary.csv",
    "data/classification/analysis_files/seq_depth_df.rds"
)

missing_files <- required_files[!file.exists(required_files)]
if(length(missing_files) > 0) {
    stop("Required files not found:\n  ", paste(missing_files, collapse = "\n  "))
}

domain_bracken <- read_csv("data/classification/taxonomic_rank_summaries/soil_microbe_db_domain_merged_lineage.csv", show_col_types = FALSE)

domain_metagenome <- domain_bracken %>% 
    separate(sample_id, into = c("compositeSampleID", "db_name"), sep = "COMP_", remove = FALSE, extra = "merge") %>% 
    mutate(compositeSampleID = str_remove(sample_id, paste0("_", db_name)),
           db_name = str_remove(db_name, "_domain_filtered"),
           sampleID_for_join = compositeSampleID) %>% 
    select(compositeSampleID, sampleID_for_join, db_name, taxon = name, percentage = fraction_total_reads) %>% 
    pivot_wider(names_from = "taxon", values_from = "percentage") %>% 
    ungroup() %>% 
    mutate(fungi_domain_metagenome = Eukaryota / (Eukaryota + Bacteria)) %>% 
    filter(db_name == "soil_microbe_db")

phyla_fungi_summary <- read_csv("data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv")
genus_fungi_summary <- read_csv("data/classification/analysis_files/soil_microbe_db_genus_fungi_summary.csv") %>% 
    mutate(genus_fungi_minus_amf_metagenome = genus_fungi_metagenome - genus_amf_metagenome)

metagenome <- left_join(domain_metagenome, genus_fungi_summary, 
                       by = c("sampleID_for_join" = "sampleID", "db_name" = "db_name")) %>%
    left_join(phyla_fungi_summary, by = c("sampleID_for_join" = "sampleID")) %>%
    left_join(readRDS("data/classification/analysis_files/seq_depth_df.rds"), 
              by = c("sampleID_for_join" = "sampleID")) %>%
    filter(seq_depth > 1000000) %>%
    select(-sampleID_for_join)

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
