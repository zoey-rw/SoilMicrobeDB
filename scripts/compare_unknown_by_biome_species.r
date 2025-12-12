# Calculate the % classified at the species level across databases, biomes, and seq depth
library(tidyverse)
library(data.table)
library(ggpubr)
library(ggpmisc)
library(rstatix)

#Read in bracken estimates from 3 databases
species_bracken1 = fread("data/classification/taxonomic_rank_summaries/species/soil_microbe_db_filtered_species_merged.csv")
species_bracken2 = fread("data/classification/taxonomic_rank_summaries/species/gtdb_species_merged.csv")
species_bracken3 = fread("data/classification/taxonomic_rank_summaries/species/pluspf_species_merged.csv")


species_bracken = rbindlist(list(species_bracken1, species_bracken2, species_bracken3))

# Summarize by % passing filter
filter_species = species_bracken %>% group_by(sample_id) %>% 
    mutate(n_classified_reads = sum(kraken_assigned_reads)) %>% 
    group_by(sample_id, n_classified_reads) %>% 
    summarize(n_pass_filter = sum(new_est_reads)) %>% 
    separate(sample_id,  
             into = c("sampleID","db_name"), sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(sample_id, paste0("_",db_name))) %>% 
    mutate(db_name = gsub("_filtered","",db_name)) 



# Add in sequencing depth
seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds") 
pass_filter_species = left_join(filter_species, #seq_depth_df) %>% 
                              seq_depth_df %>% 
                                  select(-c(db_name, identified_reads))) %>% 
    mutate(percent_classified = n_classified_reads / seq_depth,
        percent_passing = n_pass_filter / seq_depth)


# Reshape for plotting
pass_filter_species_long = pass_filter_species %>% 
    select(sampleID,seq_depth,db_name, percent_passing, percent_classified) %>% 
    pivot_longer(cols = c(percent_passing, percent_classified), 
                 names_to = "metric") %>% 
    mutate(pretty_metric = recode(metric, 
                                  "percent_passing" = "% reads classified at\n high quality",
                                  "percent_classified" = "% reads classified")) %>% 
    mutate(siteID = substr(sampleID, 1, 4))

# Get list of samples that have been evaluated by all databases
pass_filter_species_long$common =  gsub("-COMP","",pass_filter_species_long$sampleID)
sample_count = pass_filter_species_long %>% group_by(common) %>% tally
common_samples = sample_count[sample_count$n>2,]$common


# Check for significant differences between databases
# Wrap in tryCatch to handle potential errors
db_fit <- tryCatch({
    pass_filter_species_long  %>% filter(common %in% common_samples) %>%  
        filter(metric=="percent_passing") %>% 
        group_by(pretty_metric) %>% 
        anova_test(value ~ db_name) %>% 
        add_significance()
}, error = function(e) {
    warning("ANOVA test failed: ", e$message)
    data.frame()  # Return empty dataframe
})

#### Run Tukey ###
db_tukey <- tryCatch({
    result <- pass_filter_species_long %>% filter(common %in% common_samples) %>%  
        filter(metric=="percent_passing") %>% 
        group_by(pretty_metric) %>% 
        tukey_hsd(value ~ db_name) %>% 
        add_significance()
    
    # Only add xy_position if results exist
    if(nrow(result) > 0) {
        result <- result %>% add_xy_position(step.increase = .01)
        result$y.position=sqrt(result$y.position)
    }
    result
}, error = function(e) {
    warning("Tukey test failed: ", e$message)
    # Create empty dataframe with required columns
    data.frame(pretty_metric = character(0), y.position = numeric(0), 
               group1 = character(0), group2 = character(0))
})

# Visualize classification between databases
fig2 = pass_filter_species_long  %>% 
    filter(common %in% common_samples) %>%  
    filter(metric=="percent_passing") %>% 
    ggplot(aes(y = value, x=reorder(db_name, -value), color=db_name)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("% metagenomic reads classified to species") + 
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    stat_compare_means(show.legend = F)  + 
    scale_y_sqrt() + #scale_y_sqrt() +
    {if(nrow(db_tukey) > 0 && any(db_tukey$pretty_metric == "% reads classified at\n high quality")) {
        stat_pvalue_manual(db_tukey %>%  
                               filter(pretty_metric=="% reads classified at\n high quality"),
                           hide.ns = T, y.position = .3,
                           bracket.nudge.y = .2, step.increase = .05)
    }}#+ ggtitle("Classified to the species level") + ylim(0, .25)
fig2

# Save figure
ggsave("manuscript_figures/fig2.png", fig2, width = 10, height = 8, units = "in", dpi = 300)
cat("✅ Saved figure to: manuscript_figures/fig2.png\n")

# Print stats for % classified and kept after filtering
# Note: Code below may require external soil data file
pass_filter_species_long %>% group_by(metric, db_name) %>% 
    filter(common %in% common_samples) %>% 
    summarize(mean = mean(value, na.rm=T), median = median(value, na.rm=T))



# Now add in biome information for Figure 5 - only SMDB 
soilData_subset = soilCores %>% 
    filter(compositeSampleID %in% pass_filter_species_long$sampleID) %>% 
    select(compositeSampleID, siteID, biome, horizon) %>% 
    distinct(compositeSampleID, .keep_all = T) %>% 
    filter(!is.na(compositeSampleID))

pass_filter_soilCore = left_join(pass_filter_species_long %>%
                                     mutate(compositeSampleID = sampleID), 
                                 soilData_subset)  %>% 
    filter(sampleID != "TALL_004-M-20140708-COMP" ) # this sample got corrupted at some point


# Subset for plotting
species_pass_filter = pass_filter_soilCore %>% 
    filter(db_name =="soil_microbe_db") %>% 
    filter(metric == "percent_passing") %>% 
    filter(seq_depth > 100000) 

species_pass_filter$isWetland = ifelse(species_pass_filter$biome=="Wetlands", T,F)

species_pass_filter$biome = factor(species_pass_filter$biome, levels = c("Shrubland","Herbaceous","Forest", "Cultivated","Wetlands"))

# Check for significant differences between biomes
fit <- species_pass_filter %>% 
    group_by(db_name, pretty_metric) %>% 
    anova_test(value ~ biome) %>% 
    add_significance()
#### Run Tukey ###
tukey <- species_pass_filter %>% 
    group_by(db_name, pretty_metric) %>% 
    tukey_hsd(value ~ biome) %>% 
    add_significance() %>% 
    add_xy_position(step.increase = .01)

tukey$y.position=sqrt(tukey$y.position)


c = ggplot(species_pass_filter,
         #  aes(y = value, x=reorder(biome, -value), color=biome)) +
    aes(y = value, x=biome, color=biome)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .3, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("% reads classified to species at\n high quality") + 
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    stat_compare_means(show.legend = F)  + 
    scale_y_sqrt()  + #scale_y_sqrt() +
    stat_pvalue_manual(tukey,#y.position = .21,
                       #bracket.nudge.y = .34,
                       hide.ns = T, y.position = .34,
                       bracket.nudge.y = .22, step.increase = .01)
c


# Now look at sequencing depth
b <- ggplot(species_pass_filter,
       aes(y = value, x=seq_depth)) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("% reads classified to species at\n high quality") + 
    stat_smooth(method = "loess") +
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    scale_y_sqrt() + 
    #scale_x_log10() + 
    scale_x_continuous(
        trans = "log10",
        breaks=c(1e5,1e6,2.5e6,5e6,1e7,5e7,1e8),
        labels = c("100K","1M","2.5M","5M","10M","50M","100M") )


# Proportion fungi vs sequencing depth
phyla_fungi_summary <- read_csv("data/classification/analysis_files/soil_microbe_db_phyla_fungi_summary.csv")
fungi_proportion <- left_join(species_pass_filter, phyla_fungi_summary)

a <- ggplot(fungi_proportion,
       aes(y = phyla_fungi_metagenome, x=seq_depth)) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("Fungal proportion in \nmetagenome") + 
    #xlab("% of fungi in mock community") + 
    #xlab("Database name") +
    # stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    stat_smooth(method = "loess") +
    
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    scale_y_sqrt() + 
    #tagger::tag_facets(tag_levels = "A")  +
    scale_x_continuous(
        trans = "log10",
        #trans="sqrt",
        breaks=c(1e5,1e6,2.5e6,5e6,1e7,5e7,1e8),labels = c("100K","1M","2.5M","5M","10M","50M","100M") )

left_panel=ggarrange(a, b, nrow=2, labels = c("A","B"))


ggarrange(left_panel, c, nrow=1, labels = c(NA,"C"))
