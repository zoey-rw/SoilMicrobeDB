



# From supplement of Lang et al. 2023
AM_trees = read_csv("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/misc_scripts/AM_fungi/AM_trees_lang.csv")


bracken_genus_estimates <- readRDS("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/bracken_genus_estimates.rds")


fungal_phyla = c("Ascomycota","Basidiomycota","Blastocladiomycota","Chytridiomycota","Cryptomycota","Mucoromycota","Microsporidia","Olpidiomycota","Zoopagomycota")

# incomplete list, just for data exploration
emf_list = c("Cenoccoccum","Russula","Umbelopsis","Lactifluus","Lactarius","Boletus","Suillus","Laccaria","Amanita",#"Trichoderma", 
             "Tomentella","Inocybe","Cortinarius","Cadophora","Entoloma","Hebeloma","Tricholoma","Piloderma","Scleroderma","Rhizopogon","Cantharellus","Clavulina","Ceratobasidium","Tulasnella","Hysterangium","Gymnomyces","Serendipita","Thelephora","Elaphomyces","Meliniomyces")

bracken_genus =fread("/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_genus_merged_lineage.csv", nThread = 8)


bracken_phylum =fread("/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_phylum_merged_lineage.csv", nThread = 8)



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
bracken_genus$is_emf = ifelse(bracken_genus$Genus %in% emf_list, T, F)



bracken_genus[bracken_genus$is_emf==T,]$lineage %>% unique

bracken_genus[bracken_genus$is_fungi==T,] %>% View

emf_counts_rel = bracken_genus %>% group_by(sample_id, is_fungi) %>% 
    mutate(fungi_abun = sum(fraction_total_reads)) %>% ungroup %>% 
    group_by(sample_id,is_fungi, fungi_abun, is_emf) %>% 
    summarize(emf_abun = sum(fraction_total_reads)) %>% ungroup %>% 
    filter(is_emf==T & is_fungi==T) %>% select(-c(is_emf, is_fungi)) %>% 
    mutate(emf_abun_rel = emf_abun/fungi_abun)

emf_counts_rel$db_name = "soil_microbe_db"
emf_counts_rel = emf_counts_rel %>% 
    separate(sample_id,  
             into = c("compositeSampleID","db_name"), 
             sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(compositeSampleID = str_remove(sample_id, paste0("_",db_name))) 
emf_counts_rel = 
    cbind(emf_counts_rel,
          microbialForecast:::parseNEONsampleIDs(emf_counts_rel$compositeSampleID) %>% 
              select(siteID, plotID, dateID, horizon, plot_date)) %>% 
    mutate(Site = siteID)
emf_counts_rel$`Plot Number` = substr(emf_counts_rel$plotID, 6, 8) %>% as.numeric(gsub("^0|^00", "", .))



emf_rel_EM_trees = merge(emf_counts_rel, AM_trees)

# Supplemental figure
p2 = ggplot(emf_rel_EM_trees %>% filter(dateID > 201501), 
            aes(x = `Proportion ECM basal area`, y = emf_abun_rel)) + 
    geom_jitter(aes(color = Site), alpha=.8, size=2) + facet_grid(~horizon) +
    stat_smooth(method="lm") +
    scale_y_log10() + 
    #scale_x_sqrt() + 
    #    geom_abline(slope=1, intercept = 0, linetype=2) + 
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=6, label.x.npc = .4)  + 
    theme_bw(base_size = 18)  +
    ylab("Proportion EMF in metagenome")  +
    xlab("Proportion EM tree basal area") 
p2


ggarrange(p1, p2, common.legend = T, nrow = 2)
