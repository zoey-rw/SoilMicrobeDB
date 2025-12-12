library(tidyverse)
library(data.table)

# Summarize fungal % at phylum or genus level

genus_bracken = fread("data/classification/taxonomic_rank_summaries/genus/soil_microbe_db_filtered_genus_merged_lineage.csv")

fungal_phyla = c("Ascomycota","Basidiomycota","Blastocladiomycota","Chytridiomycota","Cryptomycota","Mucoromycota","Microsporidia","Olpidiomycota","Zoopagomycota")


genus_bracken$amf = ifelse(grepl("Glomero", genus_bracken$lineage), T, F)
genus_bracken$fungi = ifelse(grepl(paste(fungal_phyla, collapse = "|"), genus_bracken$lineage), T, F)
genus_bracken$euk = ifelse(grepl("Eukaryota", genus_bracken$lineage), T, F)

euk_genus = genus_bracken[genus_bracken$euk==T,]$lineage %>% unique
amf_genus = genus_bracken[genus_bracken$amf==T,]$lineage %>% unique

fungi_genus = genus_bracken[genus_bracken$fungi==T,]$lineage %>% unique

amf_summary = genus_bracken %>% 
    group_by(sample_id, amf) %>% 
    summarize(genus_amf_metagenome=sum(fraction_total_reads))  %>% 
    filter(amf==T) %>% select(-amf)

fungi_summary = genus_bracken %>% 
    group_by(sample_id, fungi) %>% 
    summarize(genus_fungi_metagenome=sum(fraction_total_reads)) %>% 
    filter(fungi==T) %>% select(-fungi)

genus_fungi_summary = merge(amf_summary, fungi_summary, by="sample_id")

genus_fungi_summary = genus_fungi_summary %>% 
    separate(sample_id,  
             into = c("sampleID","db_name"), sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(sample_id, paste0("_",db_name))) %>% 
    mutate(db_name = gsub("_genus_filtered","",db_name)) 


write_csv(genus_fungi_summary, "data/classification/analysis_files/soil_microbe_db_genus_fungi_summary.csv")


# Phylum level
phylum_bracken = fread("data/classification/taxonomic_rank_summaries/phylum/soil_microbe_db_filtered_phylum_merged_lineage.csv")


fungal_phyla = c("Ascomycota","Basidiomycota","Blastocladiomycota","Chytridiomycota","Cryptomycota","Mucoromycota","Microsporidia","Olpidiomycota","Zoopagomycota")


phylum_bracken$amf = ifelse(grepl("Glomero", phylum_bracken$lineage), T, F)
phylum_bracken$fungi = ifelse(grepl(paste(fungal_phyla, collapse = "|"), phylum_bracken$lineage), 
                              T, F)
phylum_bracken$euk = ifelse(grepl("Eukaryota", phylum_bracken$lineage), T, F)
euk_phyla = phylum_bracken[phylum_bracken$euk==T,]$lineage %>% unique
fungi_phyla = phylum_bracken[phylum_bracken$fungi==T,]$lineage %>% unique

phyla_fungi_summary = phylum_bracken %>% 
    group_by(sample_id, euk, fungi) %>% 
    summarize(phyla_fungi_metagenome=sum(fraction_total_reads)) %>% 
    filter(fungi==T) %>% ungroup()
phyla_fungi_summary = phyla_fungi_summary %>% 
    separate(sample_id,  
             into = c("sampleID","db_name"), sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(sample_id, paste0("_",db_name))) %>% 
    mutate(db_name = gsub("_genus_filtered","",db_name)) %>% 
    select(sampleID,phyla_fungi_metagenome)


write_csv(phyla_fungi_summary, "data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv")



phylum_to_merge = phylum_bracken %>% 
    separate(sample_id,  
             into = c("compositeSampleID","db_name"), 
             sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(compositeSampleID = 
               str_remove(sample_id, paste0("_",db_name))) %>% 
    mutate(db_name = gsub("_genus_filtered","",db_name)) 

soil_to_merge = soilCores %>% select(compositeSampleID, 
                                     siteID, biome,
                                     nlcdClass, horizon, sampleBottomDepth,sampleTopDepth, litterDepth, standingWaterDepth, soilTemp)  %>% 
    filter(!is.na(compositeSampleID)) %>% 
    distinct(compositeSampleID, .keep_all = T)

phylum_to_plot = left_join(phylum_to_merge, soil_to_merge, relationship =
                               "many-to-many")


phylum_to_plot$siteID = substr(phylum_to_plot$sample_id, 1, 4)

ggplot(phylum_to_plot %>% filter(name %in% fungal_phyla),
       aes(x = name, color=biome, group=biome,
           fill=biome, y = fraction_total_reads)) + 
    #geom_point() + 
    geom_boxplot() +
    facet_wrap(~name, scales = "free") + scale_y_sqrt()


