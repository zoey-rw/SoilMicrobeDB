library(tidyverse)
library(data.table)

other_soilData <- readRDS("/projectnb/dietzelab/zrwerbin/N-cycle/data/NEON_soilCovariates_allsites.rds")


soilData <- readRDS("/projectnb/dietzelab/zrwerbin/N-cycle/neon_soil_data_2023.rds")
soilCores <- soilData$sls_soilCoreCollection
genomicSamples <- soilData$sls_metagenomicsPooling %>%
	tidyr::separate(genomicsPooledIDList, into=c("first","second","third"),sep="\\|",fill="right") %>%
	dplyr::select(genomicsSampleID,first,second,third)
genSampleExample <- genomicSamples %>%
	tidyr::pivot_longer(cols=c("first","second","third"),values_to = "sampleID") %>%
	dplyr::select(sampleID,genomicsSampleID) %>%
	drop_na()
soilCores$compositeSampleID = genSampleExample[match(soilCores$sampleID,
                                                     genSampleExample$sampleID),]$genomicsSampleID
soilCores <- soilCores %>% filter(!is.na(compositeSampleID))


domain_bracken1=fread("data/classification/taxonomic_rank_summaries/domain/soil_microbe_db_filtered_domain_merged.csv")
domain_bracken = domain_bracken1
# domain_bracken2=fread("/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/pluspf_filtered_domain_merged.csv")
# domain_bracken3=fread("/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/gtdb_207_filtered_domain_merged.csv")
# domain_bracken = rbindlist(list(domain_bracken1, domain_bracken2, domain_bracken3))

filter_domain = domain_bracken %>% group_by(sample_id) %>% 
    summarize(n_pass_filter = sum(new_est_reads)) %>% 
    separate(sample_id,  
             into = c("sampleID","db_name"), sep = "COMP_", remove = F, extra = "merge") %>% 
    
    mutate(sampleID = str_remove(sample_id, paste0("_",db_name))) %>% 
    mutate(db_name = gsub("_domain_filtered","",db_name)) 

seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds") %>% 
    select(-c(db_name, identified_reads)) #%>% 
   # rename(n_classified_reads = identified_reads) 

pass_filter_domain = left_join(filter_domain, seq_depth_df) %>% 
    mutate(#percent_classified = n_classified_reads / seq_depth,
           percent_passing = n_pass_filter / seq_depth)


pass_filter_long = pass_filter_domain %>% pivot_longer(
    cols = percent_passing,
        #c(percent_passing, percent_classified), 
                                                names_to = "metric") %>% 
    mutate(pretty_metric = recode(metric, "percent_passing" = "% reads classified at\n high quality",
                                  "percent_classified" = "% reads classified")) %>% 
    mutate(siteID = substr(sampleID, 1, 4))

pass_filter_long %>% group_by(metric, db_name, siteID) %>% 
    summarize(mean = mean(value), median = median(value))

pass_filter_long %>% group_by(metric, db_name) %>% 
    summarize(mean = mean(value, na.rm=T), median = median(value, na.rm=T))

pass_filter_long$common =  gsub("-COMP","",pass_filter_long$sampleID)
# common2 generated at beginning of "compare_mags.r" script
a <- pass_filter_long  %>% filter(common %in% common2) %>%  
    ggplot(aes(y = value, x=reorder(db_name, -value), color=db_name)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("% metagenomic reads") + 
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    stat_compare_means(show.legend = F, method="t.test")  + 
    scale_y_sqrt()  
    #ggtitle("Classified to the domain level")  + 
    #ylim(0, .25)

ggarrange(a, b, )

soilCore_subset = soilCores %>% select(sampleID = compositeSampleID,plotID,
                     nlcdClass, horizon, sampleBottomDepth,sampleTopDepth, 
                     litterDepth, standingWaterDepth, soilTemp, 
                     sampleTiming, elevation) %>% distinct(sampleID, .keep_all = T) %>% 
    filter(sampleID %in% pass_filter_long$sampleID)  %>% 
    mutate(isForest = ifelse(grepl("Forest", nlcdClass), "forest habitat", "non forest habitat")) %>% 
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

pass_filter_soilCore = left_join(pass_filter_long, soilCore_subset) %>% 
    mutate(horizon_pretty = recode(horizon, "M" = "Mineral horizon (30-100cm)",
                                   "O" = "Organic horizon (0-30cm)"))

other_soilData_subset = other_soilData %>% 
    filter(pooled_sampleID %in% pass_filter_long$sampleID) %>% 
    select(sampleID = pooled_sampleID, soilAmmoniumNugPerGram, soilNitrateNitriteNugPerGram, soilInCaClpH, soilMoisture, nitrogenPercent, organicCPercent, CNratio, soilInorganicNugPerGram) %>% 
    distinct(sampleID, .keep_all = T)

pass_filter_soilCore = left_join(pass_filter_soilCore, other_soilData_subset)
    
ggplot(pass_filter_soilCore %>% 
           filter(db_name =="soil_microbe_db") %>% 
           filter(seq_depth > 1000000),
       aes(y = value, x=reorder(isForest, -value), color=isForest)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
    position=position_jitter(width = .1, height=0.01)) + 

    theme_bw(base_size = 18) + 
    facet_wrap(~pretty_metric, drop=T, scales="free")  +
    ylab("% metagenomic reads") + 
    #xlab("% of fungi in mock community") + 
    #xlab("Database name") +
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    stat_compare_means(show.legend = F)  + 
    scale_y_sqrt() + 
    tagger::tag_facets(tag_levels = "A")


ggplot(pass_filter_soilCore %>% 
           filter(metric != "percent_classified_passing") %>% 
           filter(seq_depth > 1000000),
       aes(y = value, x=reorder(biome, -value), color = pretty_metric)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        #position=position_jitter(width = .1, height=0.01)) +
        position=position_jitterdodge(jitter.height = .01, 
                                      jitter.width = .1, dodge.width = .75)) +
    theme_bw(base_size = 18) + 
    facet_wrap(~pretty_metric, drop=T, scales="free")  +
    ylab("% reads") + 
    #xlab("% of fungi in mock community") + 
    xlab("Database name") +
    guides(color=guide_legend(NULL)) + #scale_y_sqrt() +
    stat_pvalue_manual(tukey,
                       hide.ns = T)


fit <- pass_filter_soilCore %>% 
    filter(metric != "percent_classified_passing") %>% 
    filter(seq_depth > 1000000) %>% 
    group_by(db_name, pretty_metric) %>% 
    anova_test(value ~ biome) %>% 
    add_significance()

#### Run Tukey ###
tukey <- pass_filter_soilCore %>% 
    filter(metric != "percent_classified_passing") %>% 
    filter(seq_depth > 1000000) %>% 
    group_by(db_name, pretty_metric) %>% 
    tukey_hsd(value ~ biome) %>% 
    add_significance() %>% 
    add_xy_position(step.increase = .01)






ggplot(pass_filter_soilCore %>% 
           filter(metric == "percent_passing") %>% 
           filter(seq_depth > 1000000),
       aes(y = value, x=reorder(nlcdClass, -value), color = horizon)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        #position=position_jitter(width = .1, height=0.01)) +
        position=position_jitterdodge(jitter.height = .01, 
                                      jitter.width = .1, dodge.width = .75)) +
    theme_bw(base_size = 18) + 
    facet_wrap(db_name~pretty_metric, drop=T, scales="free")  +
    ylab("% reads") + 
    #xlab("% of fungi in mock community") + 
    xlab("Database name") +
    guides(color=guide_legend(NULL)) 



ggplot(pass_filter_soilCore %>% 
           filter(metric == "percent_classified") %>% 
           filter(seq_depth > 1000000),
       aes(y = value, x=sampleBottomDepth)) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        #position=position_jitter(width = .1, height=0.01)) +
        position=position_jitter()) +
    theme_bw(base_size = 18) + 
    facet_wrap(~pretty_metric, drop=T, scales="free")  +
    ylab("% reads") + 
    #xlab("% of fungi in mock community") + 
    xlab("Database name") +
    guides(color=guide_legend(NULL)) + geom_smooth()



ggplot(pass_filter_soilCore %>% 
           filter(metric == "percent_passing") %>% 
           filter(seq_depth > 1000000),
       aes(y = value, x=soilInCaClpH)) +
    geom_point(
        size=2, alpha=.4, show.legend = F,
        #position=position_jitter(width = .1, height=0.01)) +
        position=position_jitter()) + geom_smooth() + 
    facet_grid(~db_name) + scale_y_sqrt()

genus_bracken = fread("data/classification/taxonomic_rank_summaries/genus/soil_microbe_db_filtered_genus_merged.csv")
