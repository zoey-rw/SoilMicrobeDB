library(tidyverse)
library(ggpubr)
library(data.table)
library(ggpmisc)
library(rstatix)
library(ggstatsplot)


source("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/helper_functions.r")

source("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/comets_shinyapp_example/source.R")

alaska_sites = c("BONA","DEJU","HEAL","TOOL","BARR")
taiga_sites = c("BONA","DEJU","HEAL")
tundra_sites = c("TOOL","BARR")


# Read in sequencing (read) depth
# Generated by scripts/summarize_outputs/calculate_sequencing_depth.r
seq_depth_df <- readRDS("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/seq_depth_df.rds") %>% 
    mutate(sample_id = paste0(sampleID, "_soil_microbe_db_filtered"))

# Read in % of reads passing filter
# Generated by scripts/summarize_outputs/summarize_pct_filtered.r
filter_results = read_csv("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/pass_filter_summary.csv") %>% filter(metric=="percent_classified" & db_name == "soil_microbe_db") %>% select(-1)

# Read in MAG abundances
# Generated by scripts/add_mag_source_info.r
is_mag_df = read_csv("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/abundance_MAGs.csv") %>% select(-1) %>% 
    mutate(percent = sum*100,
           is_alaska = ifelse(grepl("Alaska", custom_biome),
                              "Alaska sites", "Non-Alaska sites"))

is_mag_df <- left_join(is_mag_df, filter_results) 

# Perform Tukey test over all biome groups (separating Alaska vs not)
is_mag_df_tukey = is_mag_df %>% ungroup %>% filter(!is.na(sum)) 
mag_tukey = microbialForecast::tukey(is_mag_df_tukey$custom_biome,is_mag_df_tukey$sum*100,y.offset = 0) %>% 
    mutate(y.position=tot + 2,
           custom_biome = x,
           is_alaska = ifelse(grepl("Alaska", mag_tukey$custom_biome),"Alaska sites", "Non-Alaska sites"),
           biome = str_replace(custom_biome, " Alaska", ""))

    
# Figure 6 
fig_6a = ggplot(is_mag_df, aes(y = percent, x=reorder(biome, -sum))) +
    geom_violin(draw_quantiles = c(.5), show.legend = F) +
    geom_point(aes(color=biome, fill=biome),
        size=2, alpha=.4, show.legend = F,
    position=position_jitter(width = .05, height = .01)) +
    ylab("Rel. abundance of MAGs") + 
    xlab("Biome") +
    geom_text(data = mag_tukey, 
              aes(x = biome, y = tot + 5, 
                  label = Letters_Tukey), show.legend = F, color = 1, size =6) +
    facet_grid(~is_alaska, scales="free_x", space="free") + 
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 310, vjust = .5, hjust=0),
          axis.title = element_text(face = "bold"), 
                                axis.title.y.right = 
              element_text(size = 20), 
          legend.title = element_text(face = "bold"), 
                                plot.title = element_text(size = 22, face = "bold"), 
          strip.text = element_text(face = "bold")) +  
    theme(plot.margin = unit(c(1,25,1,10), "pt"))  +
    scale_color_brewer(palette="Set1")
fig_6a


formula = y ~ poly(x, 3, raw = TRUE) 
fig_6b = ggplot(is_mag_df, 
                aes(y = percent, x=value, color=is_alaska)) +
    geom_point(aes( color=is_alaska),
               size=2, alpha=.4) +
    #           position=position_jitter(width = .1, height=0.01)) + 
    theme(legend.position = "top") +
ylab("Rel. abundance of MAGs") + 
    xlab("% of high-confidence classifications") + 
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 310, vjust = .5, hjust=0),
          axis.title = element_text(face = "bold"), 
          axis.title.y.right = 
              element_text(size = 20), 
          legend.text = element_text(size = 18), 
          plot.title = element_text(size = 22, face = "bold"), 
          strip.text = element_text(face = "bold"),
          legend.position = "right") + #, legend.position.inside = c(.8, .97))  +
    stat_poly_line(formula = formula, se = F, show.legend = F) +
    stat_poly_eq(aes(label =  paste(after_stat(rr.label), "*\", \"*", 
                                    after_stat(p.value.label), "*\".\"",
                                    sep = "")),
                 formula = formula, size = 6,
                 label.x = "right", label.y = "top", show.legend = F) +  
guides(color=guide_legend(NULL)) +  theme(plot.margin = unit(c(10,25,1,10), "pt")) +
    scale_color_brewer(palette="Set1")
fig_6b


fig6 = ggarrange(fig_6a, fig_6b, common.legend = F, heights = c(1,1), nrow=2,
          labels = c("A","B"))

fig6
fig6 %>% 
    ggexport(height = 1000, width = 800, filename = "/projectnb/frpmars/soil_microbe_db/manuscript_figures/fig6.png")


# below - old

# positive relationship MAGs ~ latitude 
a <- ggplot(is_mag_df, aes(y = sum*100, x=decimalLatitude)) +
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
