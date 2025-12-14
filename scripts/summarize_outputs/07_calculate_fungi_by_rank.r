# Calculate fungal abundances at genus and phylum taxonomic ranks
# Input: Merged lineage CSV files from 02_summarize_taxonomic_ranks.r
#        soil_microbe_db_genus_merged_lineage.csv and soil_microbe_db_phylum_merged_lineage.csv
# Output: soil_microbe_db_genus_fungi_summary.csv (includes AMF and total fungi)
#         soil_microbe_db_phyla_fungi_summary.csv (phylum-level fungal abundance)
#         Identifies fungi by lineage matching against known fungal phyla

library(tidyverse)
library(data.table)

fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")

genus_bracken <- fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_genus_merged_lineage.csv") %>%
    mutate(amf = grepl("Glomero", lineage),
           fungi = grepl(paste(fungal_phyla, collapse = "|"), lineage))

amf_summary <- genus_bracken %>%
    group_by(sample_id, amf) %>%
    summarize(genus_amf_metagenome = sum(fraction_total_reads), .groups = "drop") %>%
    filter(amf == TRUE) %>%
    select(sample_id, genus_amf_metagenome)

fungi_summary <- genus_bracken %>%
    group_by(sample_id, fungi) %>%
    summarize(genus_fungi_metagenome = sum(fraction_total_reads), .groups = "drop") %>%
    filter(fungi == TRUE) %>%
    select(sample_id, genus_fungi_metagenome)

genus_fungi_summary <- inner_join(amf_summary, fungi_summary, by = "sample_id") %>%
    separate(sample_id, into = c("sampleID", "db_name"), sep = "COMP_", remove = FALSE, extra = "merge") %>%
    mutate(sampleID = str_remove(sample_id, paste0("_", db_name)),
           db_name = gsub("_genus_filtered", "", db_name)) %>%
    select(sampleID, db_name, genus_amf_metagenome, genus_fungi_metagenome)

dir.create("data/classification/analysis_files", recursive = TRUE, showWarnings = FALSE)
write_csv(genus_fungi_summary, "data/classification/analysis_files/soil_microbe_db_genus_fungi_summary.csv")

phylum_bracken <- fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_phylum_merged_lineage.csv") %>%
    mutate(fungi = grepl(paste(fungal_phyla, collapse = "|"), lineage))

phyla_fungi_summary <- phylum_bracken %>%
    group_by(sample_id, fungi) %>%
    summarize(phyla_fungi_metagenome = sum(fraction_total_reads), .groups = "drop") %>%
    filter(fungi == TRUE) %>%
    separate(sample_id, into = c("sampleID", "db_name"), sep = "COMP_", remove = FALSE, extra = "merge") %>%
    mutate(sampleID = str_remove(sample_id, paste0("_", db_name)),
           db_name = gsub("_genus_filtered", "", db_name)) %>%
    select(sampleID, phyla_fungi_metagenome)

write_csv(phyla_fungi_summary, "data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv")
