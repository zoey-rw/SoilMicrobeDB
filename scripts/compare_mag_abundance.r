library(tidyverse)
library(ggpubr)
library(data.table)
library(ggpmisc)
library(rstatix)

source("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/helper_functions.r")

source("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/comets_shinyapp_example/source.R")
alaska_sites = c("BONA","DEJU","HEAL","TOOL","BARR")
taiga_sites = c("BONA","DEJU","HEAL")
tundra_sites = c("TOOL","BARR")

is_mag_df = read_csv("/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/abundance_MAGs.csv")


mag_tukey <- is_mag_df %>% ungroup %>% 
    #group_by(pretty_metric) %>% 
    tukey_hsd(sum ~ biome) %>% 
    add_significance() %>% 
    add_xy_position(step.increase = .1)

mag_tukey$y.position=sqrt(mag_tukey$y.position * 100)

is_mag_df_tukey = is_mag_df %>% ungroup %>% filter(!is.na(sum)) 
mag_tukey = microbialForecast::tukey(is_mag_df_tukey$custom_biome,is_mag_df_tukey$sum*100,y.offset = 0) 
#mag_tukey$y.position=sqrt(mag_tukey$tot)
mag_tukey$y.position=mag_tukey$tot + 2
mag_tukey$custom_biome = mag_tukey$x
mag_tukey$is_alaska = ifelse(grepl("Alaska", mag_tukey$custom_biome),"Alaska site", "Not Alaska site")

# Figure 6 
fig_6a = ggplot(is_mag_df, aes(y = sum*100, x=reorder(custom_biome, -sum))) +
                    # color=is_alaska)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(aes(color=is_alaska),
        size=2, alpha=.4, show.legend = F,
        #position=position_jitterdodge(jitter.width = .1, jitter.height = .01, dodge.width = .75)) +
    position=position_jitter(width = .1, height = .01)) +
    theme_bw(base_size = 16) + 
    ylab("Rel. abundance of MAGs at species-level") + 
    xlab("Biome") +
    geom_text(data = mag_tukey, 
              aes(x = custom_biome, y = tot + .2, 
                  label = Letters_Tukey), show.legend = F, color = 1, size =5) +
    # stat_pvalue_manual(mag_tukey,#y.position = .21,
    #                    #bracket.nudge.y = .34,
    #                    hide.ns = T, #y.position = 34,
    #                    #bracket.nudge.y = .22, 
    #                    step.increase = .01) +
    #scale_y_sqrt()  +
    #ggtitle("Abundance of uncultured genomes (MAGs)") + 
    scale_color_discrete(name = "") + 
    facet_grid(~is_alaska, scales="free", space="free") + 
    theme(axis.text.x = element_text(angle = 310, vjust = .5, hjust=0))
fig_6a

# fig_6b = ggplot(is_mag_df, aes(y = sum*100, x=reorder(is_alaska, -sum),
#                           color=is_alaska)) +
#     geom_boxplot(outlier.shape = NA, show.legend = F) +
#     geom_point(aes(color=is_alaska),
#                size=2, alpha=.4, show.legend = F,
#                #position=position_jitterdodge(jitter.width = .1, jitter.height = .01, dodge.width = .75)) +
#                position=position_jitter(width = .1, height = .01)) +
#     theme_bw(base_size = 16) + 
#     ylab("Rel. abundance of MAGs at species-level") + 
#     xlab(NULL) +
#     guides(color=guide_legend(NULL)) +
#     stat_compare_means(show.legend = F, method="t.test", hide.ns = T,label.y.npc = .9)  + 
# #                    ))  +     
#     scale_y_sqrt()  +
#     #ggtitle("Abundance of uncultured genomes (MAGs)") + 
#     scale_color_discrete(name = "")

fig_6c = ggplot(is_mag_df, aes(y = sum*100, x=value, color=is_alaska)) +
    geom_point(aes( color=is_alaska),
               size=2, alpha=.4) +
    #           position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 16) + 
    theme(legend.position = "top") +
ylab("Relative abundance of MAGs at species-level") + 
    xlab("% of high-confidence classifications") + 
    guides(color=guide_legend(NULL)) +
    #stat_compare_means(show.legend = F, method="t.test")  + 
    scale_y_sqrt() + 
    stat_smooth(method="lm", se = F, show.legend = F) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=7)  #+
   # facet_wrap(~biome)
fig_6c

# grob1 = ggarrange(fig_6a, fig_6b, nrow=1, #common.legend = T,
#                   widths = c(3,1), labels = c("A","B"))
# ggarrange(grob1,fig_6c, #common.legend = T, 
#           nrow=2, labels = c("","","C"))

ggarrange(fig_6a, fig_6c, common.legend = F,widths = c(2,1),  labels = c("A","B"))

a <-ggplot(is_novel_mag_df, aes(y = sum*100, x=decimalLatitude)) +
    geom_point(aes( color=siteID),
        size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("% metagenomic reads") + 
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    #stat_compare_means(show.legend = F, method="t.test")  + 
    scale_y_sqrt() + 
    ggtitle("Abundances \nassigned ONLY from novel MAGs") +
    stat_smooth(method="lm", show.legend = F) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=7) 


b <- ggplot(is_mag_df, aes(y = sum*100, x=decimalLatitude)) +
    geom_point(aes( color=biome),
               size=2, alpha=.4, #show.legend = F,
               position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("Relative abundance of MAGs at species-level") + 
    xlab("Latitude") +
    guides(color=guide_legend(NULL)) +
    #stat_compare_means(show.legend = F, method="t.test")  + 
    scale_y_sqrt() + #ggtitle("Abundances \nassigned from MAGs") +
stat_smooth(method="lm", show.legend = F) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=7) 

ggarrange(a,b)


a <- ggplot(is_mag_df, aes(y = sum*100, x=reorder(siteID, -sum), 
                            color=biome)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("% metagenomic reads") + 
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    #stat_compare_means(show.legend = F, method="t.test")  + 
    scale_y_sqrt()  + ggtitle("Species-level abundances \nassigned to MAGs") 

b <- ggplot(is_spire_smag_df, aes(y = sum*100, x=reorder(siteID, -sum), 
                                  color=siteID)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("% metagenomic reads") + 
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    #stat_compare_means(show.legend = F, method="t.test")  + 
    scale_y_sqrt() + ggtitle("Species-level abundances \nassigned to SPIRE/SMAG genomes") 

c <- ggplot(is_novel_mag_df, aes(y = sum*100, x=reorder(siteID, -sum), 
                                 color=siteID)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    geom_point(#aes(color=fungal),
        size=2, alpha=.4, show.legend = F,
        position=position_jitter(width = .1, height=0.01)) + 
    theme_bw(base_size = 18) + 
    ylab("% metagenomic reads") + 
    xlab(NULL) +
    guides(color=guide_legend(NULL)) +
    #stat_compare_means(show.legend = F, method="t.test")  + 
    scale_y_sqrt() + ggtitle("Species-level abundances \nassigned ONLY from novel genomes") 
ggarrange(a,b,c)

#lineage_df =  split_taxonomy_ncbi(bracken_with_lineage$lineage)
lineage_df =  lapply(bracken_with_lineage$lineage, 
                     function(x) {
                         phyloseq::parse_taxonomy_qiime(x) %>% #as.data.frame() %>% 
                             t()  %>% as.data.frame()
                         }) %>% 
    data.table::rbindlist()

bracken_with_lineage_full = cbind.data.frame(bracken_with_lineage, lineage_df)
bracken_with_lineage_full$fungi = ifelse(bracken_with_lineage_full$Phylum %in% fungal_phyla, T, F)
bracken_with_lineage_full$pluspf_fungi = ifelse(bracken_with_lineage_full$fungi==T & bracken_with_lineage_full$novel_fungi==F, T, F)



fungal_phyla = c("Ascomycota","Basidiomycota","Blastocladiomycota","Chytridiomycota","Cryptomycota","Mucoromycota","Microsporidia","Olpidiomycota","Zoopagomycota")

fungi_count = phylum_output %>% filter(db_name %in% c("PlusPF","PlusPFP8","soil_microbe_db")) %>%
	group_by(plot_date,db_name, sampleID, seq_depth,sampleID_orig) %>%
	filter(taxon %in% fungal_phyla) %>% summarize(total_fungi=sum(percentage))
