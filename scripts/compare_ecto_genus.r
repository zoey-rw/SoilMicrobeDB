
# Comparing abundance of ECM genera in PlusPF and Soil Microbe DB

library(tidyverse)
library(data.table)
library(ggpubr)

# Sequencing depth 
seq_depth_df <- readRDS("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/seq_depth_df.rds") %>% 
    select(-c(db_name, identified_reads)) %>% rename("compositeSampleID" = "sampleID")

# Using a few sites where all the samples definitely ran
# Bracken abundances

genus_bracken1 = fread("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/soil_microbe_db_filtered_genus_merged.csv") %>% mutate(db_name="soil_microbe_db")
genus_bracken2 = fread("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/pluspf_filtered_genus_merged.csv") %>% mutate(db_name="pluspf")
#genus_bracken3 = fread("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/gtdb_207_filtered_genus_merged.csv") %>% mutate(db_name="gtdb_207")
genus_bracken = rbindlist(list(genus_bracken1, genus_bracken2)) %>% 
    mutate(compositeSampleID = str_remove(sample_id, paste0("_",db_name,"_genus_filtered"))) 

# Filtere to samples evaluated by both databases
keep_samples = intersect(genus_bracken[genus_bracken$db_name=="pluspf",]$compositeSampleID,
                         genus_bracken[genus_bracken$db_name=="soil_microbe_db",]$compositeSampleID) %>% unique

genus_df = genus_bracken %>% filter(compositeSampleID %in% keep_samples) %>% 
    mutate(siteID = substr(compositeSampleID, 1, 4))

genus_df = left_join(genus_df, seq_depth_df) 


ecto = genus_df %>% filter(name %in% c("Cenococcum", "Gigaspora","Aspergillus",
                                        "Lactifluus","Lactarius","Russula","Xylaria","Trichoderma"))

ecto %>% group_by(db_name,siteID, name) %>% summarize(taxon_mean = mean(fraction_total_reads, na.rm=T))

site_mean = genus_df %>% group_by(db_name,siteID, name) %>% summarize(taxon_mean = mean(fraction_total_reads, na.rm=T))

# Supp figure
ggplot(ecto  %>% 
          # filter(seq_depth > 1000000) %>% 
           filter(db_name %in% c("pluspf","soil_microbe_db")) , 
       aes(x=reorder(name, -fraction_total_reads), y = fraction_total_reads, color=db_name)) + 
    ylab("Abundance") +
    xlab("Genus") +
    geom_point(alpha=.5, position=position_dodge(width = .75)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    scale_y_sqrt() + 
    scale_color_discrete(name = "Database") + 
    stat_compare_means(hide.ns = T, label = "p.signif", #label.y = .1, 
                       label.y.npc = .7, size=7)

