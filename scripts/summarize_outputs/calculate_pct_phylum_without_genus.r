
# Read in the summary files from all SoilMicrobeDB output files
lineage_phylum = fread("data/classification/taxonomic_rank_summaries/phylum/soil_microbe_db_filtered_phylum_merged_lineage.csv")
lineage_genus = fread("data/classification/taxonomic_rank_summaries/genus/soil_microbe_db_filtered_genus_merged_lineage.csv")

fungal_phyla = c("Ascomycota","Basidiomycota","Blastocladiomycota",
                 "Chytridiomycota","Cryptomycota","Mucoromycota",
                 "Microsporidia","Olpidiomycota","Zoopagomycota")

# Sum the number of phylum-level reads assigned to the fungal phyla listed above
lineage_phylum_fungi = lineage_phylum %>% filter(name %in% fungal_phyla) %>% 
    group_by(sample_id) %>% 
    summarize(phylum_kraken_fungi = sum(kraken_assigned_reads)) %>% 
    separate(sample_id, 
             into = c("sampleID","db_name"), 
             sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(sample_id, paste0("_",db_name))) %>% 
    select(-c(sample_id, db_name))

# Sum the number of genus-level reads assigned to the fungal phyla listed above
lineage_genus_fungi = lineage_genus %>% filter(grepl(paste0(fungal_phyla, collapse="|"), lineage)) %>% group_by(sample_id) %>% 
    summarize(genus_kraken_fungi = sum(kraken_assigned_reads))  %>% 
    separate(sample_id, 
             into = c("sampleID","db_name"), 
             sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(sample_id, paste0("_",db_name))) %>% 
    select(-c(sample_id, db_name))

# Merge together the read counts assigned to fungal phyla vs fungal genera
lineage_merged = merge(lineage_phylum_fungi, lineage_genus_fungi)

# Calculate proportion of genus reads 
lineage_merged$phylum_identified = 1-(lineage_merged$genus_kraken_fungi/lineage_merged$phylum_kraken_fungi)

hist(lineage_merged$phylum_identified )
median(lineage_merged$phylum_identified )
       
