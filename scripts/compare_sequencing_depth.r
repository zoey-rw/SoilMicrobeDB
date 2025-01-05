

library(neonUtilities)
library(tidyverse)
library(broom)
library(ggpubr)
library(SimplyAgree)

seq_depth_df <- readRDS("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/seq_depth_df.rds") %>% 
    #dplyr::rename(compositeSampleID = sampleID) %>% 
    select(-c(db_name, identified_reads))



bracken_domain_estimates <- readRDS("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/bracken_domain_estimates.rds") 
bracken_fb_ratio = bracken_domain_estimates %>%
    mutate(percentage = percentage*.01) %>%
    pivot_wider(names_from = taxon, values_from = percentage)  %>%
    mutate(#fb_ratio = Eukaryota/Bacteria,
        fb_ratio = Eukaryota/(Bacteria+Archaea),
        f_proportion = Eukaryota/(Bacteria+Archaea+Eukaryota) ) %>% 
    mutate(plotID = substr(samp_name, 1, 8)) %>% 
    filter(f_proportion < .8)
bracken_fb_ratio = left_join(bracken_fb_ratio, seq_depth_df) 




library(ggpmisc)

keep_samples = intersect(bracken_fb_ratio[bracken_fb_ratio$db_name=="pluspf",]$sampleID,
                         bracken_fb_ratio[bracken_fb_ratio$db_name=="soil_microbe_db",]$sampleID)

# Supp figure
ggplot(bracken_fb_ratio %>% 
           filter(db_name %in% c("pluspf","soil_microbe_db") & 
                      sampleID %in% keep_samples) %>% 
           mutate(db_name = recode(db_name, "soil_microbe_db" = "SoilMicrobeDB",
                                   "pluspf" = "PlusPF",
                                   "gtdb_207" = "GTDB r207")), 
       aes(x=seq_depth, y = Eukaryota, color=db_name)) + 
    geom_point(alpha=.5) +
    theme_bw() +
    stat_smooth(method="lm", se= FALSE) +
    stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "*`,`~")), 
                          show.legend = FALSE, size=7)  +
    scale_x_log10() +
    xlab("Sequencing depth (reads)") +
    ylab("Relative abundance of eukaryota") +
    scale_y_log10() + scale_color_discrete(name = "Database")


# Supp figure
ggplot(bracken_fb_ratio %>% 
           filter(seq_depth > 1000000) %>% 
           filter(db_name %in% c("pluspf","soil_microbe_db") & 
                      compositeSampleID %in% keep_samples), 
       aes(x=db_name, y = Eukaryota, color=db_name)) + 
    ylab("Relative abundance of Eukaryota") +
    xlab("Database name") +
    geom_point(alpha=.5, show.legend = F,
               position=position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    theme_bw(base_size = 18) +
    scale_y_log10() + 
    scale_color_discrete(name = "Database") + 
    stat_compare_means(show.legend = F,
                       hide.ns = T, label = "p.signif", #label.y = .1, 
                       label.y.npc = .6, label.x.npc = .5, 
                       size=7)


bracken_fb_ratio %>% filter(db_name %in% c("soil_microbe_db","pluspf") & 
                                compositeSampleID %in% keep_samples) %>% 
    	filter(seq_depth > 1000000) %>%
    #	filter(taxon %in% c("Archaea","Bacteria","Eukaryota")) %>% 
    group_by(db_name) %>% 
    dplyr::summarize( median(Eukaryota)) %>% View






ggplot(bracken_fb_ratio %>% 
           filter(seq_depth > 10000) %>% 
           filter(db_name %in% c("pluspf","soil_microbe_db")),
       aes(y = Eukaryota, x=seq_depth,color=db_name)) +
    geom_point(
        size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("Fungal proportion in \nmetagenome") + 
    xlab("Sample sequencing depth (read count)") +
    stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    #stat_smooth(method = "loess") +
    
    guides(color=guide_legend(NULL)) +
    scale_y_log10() + 
    #tagger::tag_facets(tag_levels = "A")  +
    scale_x_continuous(
        trans = "log10",
        #trans="sqrt",
        breaks=c(1e5,1e6,2.5e6,5e6,1e7,5e7,1e8),labels = c("100K","1M","2.5M","5M","10M","50M","100M") ) +
    
    stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
             # method="pearson", 
            # label.y = 0, 
             size=7)
