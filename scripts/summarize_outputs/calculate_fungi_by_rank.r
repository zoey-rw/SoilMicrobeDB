# Calculate fungal abundances at genus and phylum taxonomic ranks
# Output: Summary files for fungal abundances at different taxonomic levels

library(tidyverse)
library(data.table)

# Define fungal phyla for identification
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")

# ============================================================================
# Genus-level fungal summary
# ============================================================================

genus_bracken <- fread("data/classification/taxonomic_rank_summaries/genus/soil_microbe_db_filtered_genus_merged_lineage.csv")

# Identify fungal and AMF taxa based on lineage
genus_bracken <- genus_bracken %>%
    mutate(
        amf = grepl("Glomero", lineage),
        fungi = grepl(paste(fungal_phyla, collapse = "|"), lineage)
    )

# Calculate AMF abundance by sample
amf_summary <- genus_bracken %>%
    group_by(sample_id, amf) %>%
    summarize(genus_amf_metagenome = sum(fraction_total_reads), .groups = "drop") %>%
    filter(amf == TRUE) %>%
    select(sample_id, genus_amf_metagenome)

# Calculate total fungal abundance by sample
fungi_summary <- genus_bracken %>%
    group_by(sample_id, fungi) %>%
    summarize(genus_fungi_metagenome = sum(fraction_total_reads), .groups = "drop") %>%
    filter(fungi == TRUE) %>%
    select(sample_id, genus_fungi_metagenome)

# Merge AMF and fungal summaries
genus_fungi_summary <- inner_join(amf_summary, fungi_summary, by = "sample_id") %>%
    separate(sample_id, into = c("sampleID", "db_name"), sep = "COMP_", remove = FALSE, extra = "merge") %>%
    mutate(
        sampleID = str_remove(sample_id, paste0("_", db_name)),
        db_name = gsub("_genus_filtered", "", db_name)
    ) %>%
    select(sampleID, db_name, genus_amf_metagenome, genus_fungi_metagenome)

write_csv(genus_fungi_summary, "data/classification/analysis_files/soil_microbe_db_genus_fungi_summary.csv")
cat("✅ Saved genus-level fungal summary to: data/classification/analysis_files/soil_microbe_db_genus_fungi_summary.csv\n")

# ============================================================================
# Phylum-level fungal summary
# ============================================================================

phylum_bracken <- fread("data/classification/taxonomic_rank_summaries/phylum/soil_microbe_db_filtered_phylum_merged_lineage.csv")

# Identify fungal taxa based on lineage
phylum_bracken <- phylum_bracken %>%
    mutate(fungi = grepl(paste(fungal_phyla, collapse = "|"), lineage))

# Calculate fungal abundance by sample at phylum level
phyla_fungi_summary <- phylum_bracken %>%
    group_by(sample_id, fungi) %>%
    summarize(phyla_fungi_metagenome = sum(fraction_total_reads), .groups = "drop") %>%
    filter(fungi == TRUE) %>%
    separate(sample_id, into = c("sampleID", "db_name"), sep = "COMP_", remove = FALSE, extra = "merge") %>%
    mutate(
        sampleID = str_remove(sample_id, paste0("_", db_name)),
        db_name = gsub("_genus_filtered", "", db_name)
    ) %>%
    select(sampleID, phyla_fungi_metagenome)

write_csv(phyla_fungi_summary, "data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv")
cat("✅ Saved phylum-level fungal summary to: data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv\n")
