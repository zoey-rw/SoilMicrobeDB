# Calculate proportion of fungal phylum reads not identified to genus level
# Input: Phylum and genus merged lineage CSV files from soil_microbe_db
# Output: Calculates phylum_identified = 1 - (genus_kraken_fungi / phylum_kraken_fungi)
#         Prints histogram and median of phylum_identified
#         This metric indicates how much fungal diversity is unclassified at genus level

library(tidyverse)
library(data.table)

fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota",
                 "Chytridiomycota", "Cryptomycota", "Mucoromycota",
                 "Microsporidia", "Olpidiomycota", "Zoopagomycota")

lineage_phylum <- fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_phylum_merged_lineage.csv")
lineage_genus <- fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_genus_merged_lineage.csv")

lineage_phylum_fungi <- lineage_phylum %>% 
    filter(name %in% fungal_phyla) %>% 
    group_by(sample_id) %>% 
    summarize(phylum_kraken_fungi = sum(kraken_assigned_reads), .groups = "drop") %>%
    separate(sample_id, into = c("sampleID", "db_name"), sep = "COMP_", remove = FALSE, extra = "merge") %>%
    mutate(sampleID = str_remove(sample_id, paste0("_", db_name))) %>%
    select(sampleID, phylum_kraken_fungi)

lineage_genus_fungi <- lineage_genus %>% 
    filter(grepl(paste0(fungal_phyla, collapse = "|"), lineage)) %>%
    group_by(sample_id) %>% 
    summarize(genus_kraken_fungi = sum(kraken_assigned_reads), .groups = "drop") %>%
    separate(sample_id, into = c("sampleID", "db_name"), sep = "COMP_", remove = FALSE, extra = "merge") %>%
    mutate(sampleID = str_remove(sample_id, paste0("_", db_name))) %>%
    select(sampleID, genus_kraken_fungi)

lineage_merged <- merge(lineage_phylum_fungi, lineage_genus_fungi) %>%
    mutate(phylum_identified = 1 - (genus_kraken_fungi / phylum_kraken_fungi))

hist(lineage_merged$phylum_identified)
median(lineage_merged$phylum_identified)
