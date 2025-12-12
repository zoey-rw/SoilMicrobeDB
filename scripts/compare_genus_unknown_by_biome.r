


genus_bracken1 = fread("data/classification/taxonomic_rank_summaries/genus/soil_microbe_db_filtered_genus_merged.csv")
genus_bracken2 = fread("data/classification/taxonomic_rank_summaries/genus/pluspf_filtered_genus_merged.csv")
genus_bracken3 = fread("data/classification/taxonomic_rank_summaries/genus/gtdb_207_filtered_genus_merged.csv")
genus_bracken = rbindlist(list(genus_bracken1, genus_bracken2, genus_bracken3))


filter_genus = genus_bracken %>% group_by(sample_id) %>% 
    summarize(n_pass_filter = sum(new_est_reads)) %>% 
    separate(sample_id,  
             into = c("sampleID","db_name"), sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(sample_id, paste0("_",db_name))) %>% 
    mutate(db_name = gsub("_genus_filtered","",db_name)) 

#%>% select(-c(db_name, identified_reads))



seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds") #%>% rename(n_classified_reads = identified_reads) #%>% select(-c(db_name, identified_reads))
pass_filter_genus = left_join(filter_genus, #seq_depth_df) %>% 
                              seq_depth_df %>% select(-c(db_name, identified_reads))) %>% 
    mutate(#percent_classified = n_classified_reads / seq_depth,
           percent_passing = n_pass_filter / seq_depth)



pass_filter_genus_long = pass_filter_genus %>% 
    select(sampleID,seq_depth,db_name, percent_passing)  %>%  #, percent_classified) %>% 
    pivot_longer(cols = c(percent_passing), 
                                                       names_to = "metric") %>% 
    mutate(pretty_metric = recode(metric, "percent_passing" = "% reads classified at\n high quality",
                                  "percent_classified" = "% reads classified")) %>% 
    mutate(siteID = substr(sampleID, 1, 4))

pass_filter_genus_long %>% group_by(metric, db_name) %>% 
    summarize(mean = mean(value, na.rm=T), median = median(value, na.rm=T))

pass_filter_genus_long$common =  gsub("-COMP","",pass_filter_genus_long$sampleID)


b <- pass_filter_genus_long  %>% filter(common %in% common2) %>%  
    ggplot(aes(y = value, x=reorder(db_name, -value), color=db_name)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("% metagenomic reads") + 
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    stat_compare_means(show.legend = F)  + 
    scale_y_sqrt() + ggtitle("Classified to the genus level") + ylim(0, .25)
    
pass_filter_genus_long %>% group_by(metric, db_name, siteID) %>% filter(common %in% common2) %>% 
    summarize(mean = mean(value, na.rm=T), median = median(value, na.rm=T))


# categories from https://www.mrlc.gov/data/legends/national-land-cover-database-2011-nlcd2011-legend
soilCore_subset = soilCores %>% select(sampleID = compositeSampleID,plotID,
                                       nlcdClass, horizon, sampleBottomDepth,sampleTopDepth, 
                                       litterDepth, standingWaterDepth, soilTemp, 
                                       sampleTiming, elevation) %>% 
    distinct(sampleID, .keep_all = T) %>% 
    filter(sampleID %in% pass_filter_genus_long$sampleID)  %>% 
    mutate(isForest = ifelse(grepl("Forest", nlcdClass), 
                             "forest habitat", "non forest habitat")) %>% 
    mutate(biome = recode(nlcdClass, "mixedForest" = "Forest",
                                   "evergreenForest" = "Forest",
                          "deciduousForest" = "Forest",
                          "emergentHerbaceousWetlands" = "Wetlands",
                          "woodyWetlands" = "Wetlands",
                          "dwarfScrub" = "Shrubland",
                          "shrubScrub" = "Shrubland",
                          "sedgeHerbaceous" = "Herbaceous",
                          "grasslandHerbaceous" = "Herbaceous",
                          "pastureHay" = "Cultivated",
                          "cultivatedCrops" = "Cultivated"
    ))

pass_filter_genus_soilCore = left_join(pass_filter_genus_long, soilCore_subset) %>% 
    mutate(horizon_pretty = recode(horizon, "M" = "Mineral horizon (30-100cm)",
                                   "O" = "Organic horizon (0-30cm)"))

other_soilData_subset = other_soilData %>% 
    filter(pooled_sampleID %in% pass_filter_long$sampleID) %>% 
    select(sampleID = pooled_sampleID, soilAmmoniumNugPerGram, soilNitrateNitriteNugPerGram, soilInCaClpH, soilMoisture, nitrogenPercent, organicCPercent, CNratio, soilInorganicNugPerGram) %>% 
    distinct(sampleID, .keep_all = T)

pass_filter_genus_soilCore = left_join(pass_filter_genus_soilCore, other_soilData_subset) %>% 
    filter(sampleID != "TALL_004-M-20140708-COMP" & siteID !="GUAN") 




fit <- pass_filter_genus_soilCore %>% 
    filter(db_name =="soil_microbe_db") %>% 
    filter(metric == "percent_passing") %>% 
    filter(seq_depth > 1000000) %>% 
    group_by(db_name, pretty_metric) %>% 
    anova_test(value ~ biome) %>% 
    add_significance()
#### Run Tukey ###
tukey <- pass_filter_genus_soilCore %>% 
    filter(db_name =="soil_microbe_db") %>% 
    filter(metric == "percent_passing") %>% 
    filter(seq_depth > 1000000) %>% 
    group_by(db_name, pretty_metric) %>% 
    tukey_hsd(value ~ biome) %>% 
    add_significance() %>% 
    add_xy_position(step.increase = .01)

tukey$y.position=sqrt(tukey$y.position)
ggplot(pass_filter_genus_soilCore %>% 
           filter(db_name =="soil_microbe_db") %>% 
           filter(metric == "percent_passing") %>% 
           filter(seq_depth > 1000000),
       aes(y = value, x=reorder(biome, -value), color=biome)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    facet_wrap(db_name~pretty_metric, drop=T, scales="free")  +
    ylab("% metagenomic reads") + 
    #xlab("% of fungi in mock community") + 
    #xlab("Database name") +
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    stat_compare_means(show.legend = F)  + 
    scale_y_sqrt() + 
    tagger::tag_facets(tag_levels = "A") + #scale_y_sqrt() +
    stat_pvalue_manual(tukey,y.position = .21,bracket.nudge.y = .34,
                       hide.ns = T)#, y.position = .2, bracket.nudge.y = .32)




ggplot(pass_filter_genus_soilCore %>% 
           #filter(db_name =="soil_microbe_db") %>% 
           filter(seq_depth > 1000000),
       aes(y = value, x=sampleBottomDepth)) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F) +
    
    theme_bw(base_size = 18) + 
    facet_wrap(~pretty_metric, drop=T, scales="free")  +
    ylab("% metagenomic reads") + 
    #xlab("% of fungi in mock community") + 
    #xlab("Database name") +
    xlab(NULL) +
    facet_wrap(~db_name, drop=T, scales="free")  +
    
    guides(color=guide_legend(NULL)) +
    scale_y_sqrt() + geom_smooth() + 	stat_smooth(method = "lm", se= FALSE, colour='green', formula = y~x) 
