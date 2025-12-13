
# Comparing abundance of ECM genera in PlusPF and Soil Microbe DB

library(tidyverse)
library(data.table)
library(ggpubr)

# Sequencing depth 
# Note: seq_depth_df$sampleID actually contains compositeSampleIDs (ends with "-COMP")
seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds") %>% 
    select(-c(db_name, identified_reads)) %>% 
    rename(compositeSampleID = sampleID)

# Using a few sites where all the samples definitely ran
# Bracken abundances

genus_bracken1 = fread("data/classification/taxonomic_rank_summaries/genus/soil_microbe_db_filtered_genus_merged.csv") %>% mutate(db_name="soil_microbe_db")
genus_bracken2 = fread("data/classification/taxonomic_rank_summaries/genus/pluspf_filtered_genus_merged.csv") %>% mutate(db_name="pluspf")
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



ecto_wide = ecto  %>% 
    # filter(seq_depth > 1000000) %>% 
    filter(db_name %in% c("pluspf","soil_microbe_db")) %>% 
    select(compositeSampleID, fraction_total_reads,db_name, name) %>% 
               pivot_wider(values_from = fraction_total_reads, names_from=db_name)
fig_pluspf_vs_smd_fun = ggplot(ecto_wide) +
    geom_point(aes(x = soil_microbe_db, y = pluspf),
               size=3, alpha=.4
               #position=position_dodge(width = 1)
    ) +
    geom_smooth(aes(x = soil_microbe_db, y = pluspf), 
                method = "lm",
                formula = y ~ x) +
    geom_abline(slope = 1, intercept=0, color=1, linetype=2, linewidth=1.2) +
    facet_wrap(~name,ncol=1, scales = "free")  +
    #facet_grid(rows=vars(taxon), scales = "free")  +
    theme_bw(base_size = 16) +
    ylab("PlusPF rel. abundance") +
    xlab("SoilMicrobeDB rel. abundance") + labs(color = NULL)
