library(tidyverse)
library(data.table)
library(ggpmisc)

source("scripts/helper_functions.r")


# From supplement of Lang et al. 2023
AM_trees = read_csv("reference_data/AM_trees/AM_trees_lang.csv")



bracken_genus_estimates <- readRDS("data/classification/taxonomic_rank_summaries/genus/bracken_genus_estimates.rds")


fungal_phyla = c("Ascomycota","Basidiomycota","Blastocladiomycota","Chytridiomycota","Cryptomycota","Mucoromycota","Microsporidia","Olpidiomycota","Zoopagomycota")


amf_list = c("Diversispora","Archaeospora","Paraglomus",
             "Ambispora","Scutellospora","Rhizophagus",
             "Glomus","Acaulospora","Pacispora","Sacculospora",
             "Redeckera","Corymbiglomus","Gigaspora","Dentiscutalata","Cetraspora","Racocetra","Claroideoglomus","Funneliformis","Septoglomus","Sclerocystis","Geosiphon")

bracken_genus =fread("data/classification/taxonomic_rank_summaries/genus/soil_microbe_db_filtered_genus_merged_lineage.csv", nThread = 8)


bracken_phylum =fread("data/classification/taxonomic_rank_summaries/phylum/soil_microbe_db_filtered_phylum_merged_lineage.csv", nThread = 8)



bracken_phylum$is_mucoro = ifelse(grepl("Mucoro", bracken_phylum$lineage), T, F)
bracken_phylum_estimates <- bracken_phylum %>% filter(is_mucoro)
bracken_phylum_estimates$db_name = "soil_microbe_db"

bracken_phylum_estimates = bracken_phylum_estimates %>% 
    separate(sample_id,  
             into = c("compositeSampleID","db_name"), 
             sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(compositeSampleID = str_remove(sample_id, paste0("_",db_name))) 
bracken_phylum_estimates = 
    cbind(bracken_phylum_estimates,
          microbialForecast:::parseNEONsampleIDs(bracken_phylum_estimates$compositeSampleID) %>% select(siteID, plotID, dateID, horizon, plot_date)) %>% mutate(Site = siteID)

bracken_phylum_estimates$`Plot Number` = substr(bracken_phylum_estimates$plotID, 6, 8) %>% as.numeric(gsub("^0|^00", "", .))

mucoro_counts = bracken_phylum_estimates %>% 
    group_by(siteID, Site, plotID, `Plot Number`, db_name, dateID,
             plot_date, horizon, compositeSampleID) %>% 
        summarize(AMF_percent = sum(fraction_total_reads))



lineage_df =  lapply(unique(bracken_genus$lineage),
                     function(x) {
                         phyloseq::parse_taxonomy_qiime(x) %>% #as.data.frame() %>%
                             t()  %>% as.data.frame()
                     }) %>%
    data.table::rbindlist()
lineage_df = 
    cbind.data.frame(lineage = unique(bracken_genus$lineage), lineage_df)
bracken_genus = left_join(bracken_genus, lineage_df)

bracken_genus$is_fungi = ifelse(bracken_genus$Phylum %in% fungal_phyla, T, F)
bracken_genus$is_amf = ifelse(grepl("Glomero", bracken_genus$lineage), T, F)


amf_counts_rel = bracken_genus %>% group_by(sample_id, is_fungi) %>% 
    mutate(fungi_abun = sum(fraction_total_reads)) %>% ungroup %>% 
    group_by(sample_id,is_fungi, fungi_abun, is_amf) %>% 
    summarize(amf_abun = sum(fraction_total_reads)) %>% ungroup %>% 
    filter(is_amf==T & is_fungi==T) %>% select(-c(is_amf, is_fungi)) %>% 
     mutate(amf_abun_rel = amf_abun/fungi_abun)

amf_counts_rel$db_name = "soil_microbe_db"
amf_counts_rel = amf_counts_rel %>% 
    separate(sample_id,  
             into = c("compositeSampleID","db_name"), 
             sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(compositeSampleID = str_remove(sample_id, paste0("_",db_name))) 
amf_counts_rel = 
    cbind(amf_counts_rel,
          microbialForecast:::parseNEONsampleIDs(amf_counts_rel$compositeSampleID) %>% 
              select(siteID, plotID, dateID, horizon, plot_date)) %>% 
    mutate(Site = siteID)
amf_counts_rel$`Plot Number` = substr(amf_counts_rel$plotID, 6, 8) %>% as.numeric(gsub("^0|^00", "", .))



amf_rel_AM_trees = merge(amf_counts_rel, AM_trees, by = c("Site", "Plot Number"))

# Supplemental figure
p1 = ggplot(amf_rel_AM_trees, 
            aes(x = `Proportion AM basal area`, y = amf_abun_rel)) + 
    geom_jitter(aes(color = Site), alpha=.8, size=2) +
    stat_smooth(method="lm") +
    scale_y_log10() + 
    #scale_x_sqrt() + 
    #    geom_abline(slope=1, intercept = 0, linetype=2) + 
    stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "~~~")),
                 formula = y ~ x, parse = TRUE,
                 show.legend = FALSE, size = 6, label.x.npc = .4)  + 
    theme_bw(base_size = 18)  +
    ylab("Proportion AMF in metagenome")  +
    xlab("Proportion AM tree basal area") 
p1

# Save figure
ggsave("manuscript_figures/fig5b.png", p1, width = 10, height = 8, units = "in", dpi = 300)


######
bracken_genus_estimates <- bracken_genus %>% filter(is_amf)
table(bracken_genus_estimates[bracken_genus_estimates$is_amf,]$name)

bracken_genus_estimates$db_name = "soil_microbe_db"
bracken_genus_estimates = bracken_genus_estimates %>% 
    separate(sample_id,  
             into = c("compositeSampleID","db_name"), 
             sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(compositeSampleID = str_remove(sample_id, paste0("_",db_name))) 
bracken_genus_estimates = 
    cbind(bracken_genus_estimates,
          microbialForecast:::parseNEONsampleIDs(bracken_genus_estimates$compositeSampleID) %>% 
              select(siteID, plotID, dateID, horizon, plot_date)) %>% 
    mutate(Site = siteID)
bracken_genus_estimates$`Plot Number` = substr(bracken_genus_estimates$plotID, 6, 8) %>% as.numeric(gsub("^0|^00", "", .))

table(bracken_genus_estimates[bracken_genus_estimates$is_amf,]$name)

amf_counts = bracken_genus_estimates %>% 
	group_by(siteID, Site, plotID, `Plot Number`, db_name, dateID,
	         plot_date, horizon, compositeSampleID) %>% 
    summarize(AMF_percent = sum(fraction_total_reads))

amf_counts_genus = bracken_genus_estimates %>% 
    group_by(siteID, Site, plotID, `Plot Number`, db_name, dateID,name,
             plot_date, horizon, compositeSampleID) %>% 
    summarize(AMF_percent = sum(fraction_total_reads))


# Struo2 file used for exploration/validation (optional)
# myco_struo = read_tsv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/mycocosm_published_struo.tsv")
# myco_struo[grepl("Ambispora", myco_struo$ncbi_organism_name),]
# myco_struo[grepl("Acaulo", myco_struo$ncbi_organism_name),]
# myco_struo[grepl("Diversispora", myco_struo$ncbi_organism_name),]
# myco_struo[grepl("Archaeospora", myco_struo$ncbi_organism_name),]
# myco_struo[grepl("Paraglomus", myco_struo$ncbi_organism_name),]
# myco_struo[grepl("Rhizophagus", myco_struo$ncbi_organism_name),]


amf_counts_AM_trees = merge(amf_counts, AM_trees)


mucoro_counts_AM_trees = merge(mucoro_counts, AM_trees)


genus_counts_AM_trees = merge(amf_counts_genus, AM_trees)
p1 = ggplot(mucoro_counts_AM_trees, 
            aes(x = `Proportion AM basal area`, y = AMF_percent)) + 
    geom_point(aes(color = Site)) + facet_grid(~horizon) +
    stat_smooth(method="lm") +
    scale_y_log10() + 
    #scale_x_sqrt() + 
#    geom_abline(slope=1, intercept = 0, linetype=2) + 
    stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "~~~")),
                 formula = y ~ x, parse = TRUE,
                 show.legend = FALSE, size = 6, label.x.npc = .4) 


p2 = ggplot(amf_counts_AM_trees %>% filter(dateID > 201501), 
            aes(x = `Proportion AM basal area`, y = AMF_percent)) + 
    geom_point(aes(color = Site)) + #facet_grid(~db_name) +
    stat_smooth(method="lm") +
    scale_y_log10() + 
    #scale_x_sqrt() + 
    #    geom_abline(slope=1, intercept = 0, linetype=2) + 
    stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "~~~")),
                 formula = y ~ x, parse = TRUE,
                 show.legend = FALSE, size = 6, label.x.npc = .4)

ggplot(AM_trees) + geom_point(aes(y = `Proportion AM basal area`, x = Site, color = Site))
