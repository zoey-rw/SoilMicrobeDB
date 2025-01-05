# Create list of high-abundance taxa, from mapping NEON metagenomes to SoilGenome DB

library(tidyverse)
library(pavian)
library(ggpubr)
source("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")


sampleID_pooled_key <- readRDS("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/NEON_metagenomes/NEON_pooled_sample_key.rds")
pooled_recode_list <- sampleID_pooled_key$genomicsSampleID
names(pooled_recode_list) <- sampleID_pooled_key$sampleID

soil_sample_dir <- "/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/02_bracken_output"

samp_files <- list.files(soil_sample_dir, recursive=T, pattern = "_filtered_kraken_bracken_genuses.kreport", full.names = T)
samp_names <- gsub("_filtered_kraken_bracken_genuses.kreport","",basename(samp_files)) %>% unique()
names(samp_files) <- samp_names



phylum_in = many_files_to_matrix_list(files = samp_files, filter.tax.level = "P", 
                                      include.unclassified = T, percentages = T)[[1]] %>%
    as.data.frame() %>%
    rownames_to_column("taxon") %>%
    pivot_longer(cols = where(is.numeric), names_to="samp_name",values_to="percentage")

phylum_output <- phylum_in %>% 
    separate(samp_name, 
             into = c("sampleID","db_name"), 
             sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(samp_name, paste0("_",db_name)))
phylum_output$siteID = substr(phylum_output$sampleID, 1, 4)

saveRDS(phylum_output,"/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/bracken_phylum_estimates.rds")


phylum_output <- readRDS("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/bracken_phylum_estimates.rds")

# Create list of all samples from each database
common_samples <- phylum_output %>% filter(db_name %in% c("soil_microbe_db","pluspf","gtdb_207")) %>%  distinct(sampleID, db_name) %>% 
    count(sampleID) %>% filter(n > 2) %>% select(sampleID) %>% unlist


phylum_compare = phylum_output %>% mutate(taxon = recode(taxon, 
                                                         "Pseudomonadota" = "Proteobacteria",
                              "Firmicutes" = "Bacillota",
                              #"Acidobacteriota" = "Acidobacteria",
                              #"Actinobacteria" = "Actinomycetota",
                             # "Actinobacteriota" = "Actinomycetota"
                              ))
bracken_phylum_taxon_wide = phylum_output %>% filter(sampleID %in% common_samples) %>%
    pivot_wider(names_from = "taxon", values_from = percentage, id_cols = c(db_name, sampleID)) %>%
    group_by(sampleID, db_name) %>% 
    mutate(Actinomycetota = sum(Actinomycetota, Actinobacteriota),
         #  Acidobacteriota = sum(Acidobacteriota),
           Pseudomonadota = sum(Pseudomonadota, Proteobacteria),
           Bacillota = 	sum(Bacillota, Firmicutes, Firmicutes_A,
                            Firmicutes_B,
                            Firmicutes_C,
                            Firmicutes_D,
                            Firmicutes_E,
                            Firmicutes_F,
                            Firmicutes_G)) %>% select(-c(Firmicutes_A,
                                                         Firmicutes_B,
                                                         Firmicutes_C,
                                                         Firmicutes_D,
                                                         Firmicutes_E,
                                                         Firmicutes_F,
                                                         Firmicutes_G))
bracken_phylum_long = bracken_phylum_taxon_wide %>%
    select(-c(`Classified at a higher level`, "Unclassified")) %>% 
    pivot_longer(3:ncol(.), 
                 names_to = "taxon",values_to = "percentage")  %>%
    pivot_wider(names_from = db_name, values_from = percentage)

phylum_output %>% filter(db_name %in% c("soil_microbe_db","pluspf","gtdb_207")) %>% 
    #	filter(seq_depth > 5000000) %>%
    #	filter(taxon %in% c("Archaea","Bacteria","Eukaryota")) %>% 
    #filter(taxon %in% c("Verrucomicrobiota","Bacillota")) %>% 
    group_by(db_name,taxon) %>% 
    dplyr::summarize( mean(percentage))


pluspf_taxa = bracken_phylum_long %>% ungroup %>%  filter(!is.na(pluspf)) %>% select(taxon) %>% unique %>% unlist

soil_microbe_db_taxa = bracken_phylum_long %>% ungroup %>%  filter(!is.na(soil_microbe_db)) %>% select(taxon) %>% unique %>% unlist

gtdb_taxa = bracken_phylum_long %>% ungroup %>%  filter(!is.na(soil_microbe_db)) %>% select(taxon) %>% unique %>% unlist



bracken_phylum_long %>%
    filter(taxon %in% c("Acidobacteriota","Actinomycetota"))

#gtdb207 has Actinomycetota, so do pluspf and soilmicrobedb
# Actinobacteriota is included but abundance is 0
# Actinomycetota in 

fungal_phyla = c("Ascomycota","Basidiomycota","Blastocladiomycota","Chytridiomycota","Cryptomycota","Mucoromycota","Microsporidia","Olpidiomycota","Zoopagomycota")
archaea_phyla = c("Thermoproteota","Euryarchaeota","Nitrososphaerota")
bracken_phylum_long$Domain = ifelse(bracken_phylum_long$taxon %in% fungal_phyla, "Eukaryota",
                                    ifelse(bracken_phylum_long$taxon %in% archaea_phyla, "Archaea","Bacteria"))

library(RColorBrewer)
myColors <- brewer.pal(3,"Set2")
names(myColors) <- c("Archaea","Eukaryota", "Bacteria")
colScale <- scale_colour_manual(name = "Domain",values = myColors)



also_plot = c("Bacillota",
"Bacteroidota",
"Chytridiomycota",
"Thermoproteota",
#"Cyanobacteriota",
"Euryarchaeota",
#"Myxococcota",
"Nitrososphaerota",
#"Planctomycetota",
"Spirochaetota",
#"Thermodesulfobacteriota",
#"Vulcanimicrobiota",
"Basidiomycota",
"Ascomycota","Mucoromycota",
#"Zoopagomycota",
"Verrucomicrobiota",
#"Nitrospirota",
"Actinomycetota",
"Pseudomonadota",
#"Gemmatimonadota",
"Acidobacteriota")

fig_pluspf_vs_smd_fun = ggplot(bracken_phylum_long %>% 
                                   filter(taxon %in% also_plot & Domain == "Eukaryota")) +
    
                                 #  filter(taxon %in% c("Basidiomycota",
                                                       # "Ascomycota","Mucoromycota",
                                                       # "Zoopagomycota"))) +
    geom_point(aes(x = soil_microbe_db, y = pluspf, color = Domain),
               size=3, alpha=.4,
               #position=position_dodge(width = 1)
    ) +
    geom_smooth(aes(x = soil_microbe_db, y = pluspf, color = Domain), 
                method = "lm",
                formula = y ~ x) +
    geom_abline(slope = 1, intercept=0, color=1, linetype=2, linewidth=1.2) +
    facet_wrap(~taxon,ncol=1, scales = "free")  +
    #facet_grid(rows=vars(taxon), scales = "free")  +
    theme_bw(base_size = 16) +
    ylab("PlusPF rel. abundance") +
    xlab("SoilMicrobeDB rel. abundance") + labs(color = NULL) + colScale

fig_pluspf_vs_smd_arc = ggplot(bracken_phylum_long %>% 
                                   filter(taxon %in% also_plot & Domain == "Archaea")) +
    
    #  filter(taxon %in% c("Basidiomycota",
    # "Ascomycota","Mucoromycota",
    # "Zoopagomycota"))) +
    geom_point(aes(x = soil_microbe_db, y = pluspf, color = Domain),
               size=3, alpha=.4,
               #position=position_dodge(width = 1)
    ) +
    geom_smooth(aes(x = soil_microbe_db, y = pluspf, color = Domain), 
                method = "lm",
                formula = y ~ x) +
    geom_abline(slope = 1, intercept=0, color=1, linetype=2, linewidth=1.2) +
    facet_wrap(~taxon,ncol=1, scales = "free")  +
    #facet_grid(rows=vars(taxon), scales = "free")  +
    theme_bw(base_size = 16) +
    ylab("PlusPF rel. abundance") +
    xlab("SoilMicrobeDB rel. abundance")+ labs(color = NULL) + colScale



fig_pluspf_vs_smd_bac = ggplot(bracken_phylum_long %>%
           filter(taxon %in% also_plot & Domain == "Bacteria")) +
    geom_point(aes(x = soil_microbe_db, y = pluspf, color = Domain),
               size=3, alpha=.4,
               #position=position_dodge(width = 1)
    ) +
    geom_smooth(aes(x = soil_microbe_db, y = pluspf, color = Domain), 
                method = "lm",
                formula = y ~ x) +
    geom_abline(slope = 1, intercept=0, color=1, linetype=2, linewidth=1.2) +
    facet_wrap(~taxon,ncol=1, scales = "free")  +
    #facet_grid(rows=vars(taxon), scales = "free")  +
    theme_bw(base_size = 16) +
    ylab("PlusPF rel. abundance") +
    xlab("SoilMicrobeDB rel. abundance") + labs(color = NULL) + colScale



fig_gtdb_vs_smd_bac = ggplot(bracken_phylum_long %>%
                             filter(taxon %in% also_plot & Domain == "Bacteria")) +    geom_point(aes(x = soil_microbe_db, y = gtdb_207, color = Domain),
               size=3, alpha=.4,
               #position=position_dodge(width = 1)
    ) +
    geom_smooth(aes(x = soil_microbe_db, y = gtdb_207, color = Domain), 
                method = "lm",
                formula = y ~ x) +
    geom_abline(slope = 1, intercept=0, color=1, linetype=2, linewidth=1.2) +
    facet_wrap(~taxon,ncol=1, scales = "free")  +
    #facet_grid(rows=vars(taxon), scales = "free")  +
    theme_bw(base_size = 16) +
    ylab("GTDB r207 rel. abundance") +
    xlab("SoilMicrobeDB rel. abundance")+ labs(color = NULL) + colScale



fig_gtdb_vs_smd_arc = ggplot(bracken_phylum_long %>%
                                 filter(taxon %in% also_plot & Domain == "Archaea")) +    geom_point(aes(x = soil_microbe_db, y = gtdb_207, color = Domain),
                                                                                                      size=3, alpha=.4,
                                                                                                      #position=position_dodge(width = 1)
                                 ) +
    geom_smooth(aes(x = soil_microbe_db, y = gtdb_207, color = Domain), 
                method = "lm",
                formula = y ~ x) +
    geom_abline(slope = 1, intercept=0, color=1, linetype=2, linewidth=1.2) +
    facet_wrap(~taxon,ncol=1, scales = "free")  +
    #facet_grid(rows=vars(taxon), scales = "free")  +
    theme_bw(base_size = 16) +
    ylab("GTDB r207 rel. abundance") +
    xlab("SoilMicrobeDB rel. abundance") + labs(color = NULL) + colScale


library(ggh4x)
library(scales)


fixed_scales = facetted_pos_scales(x = list(
    taxon == "Acidobacteriota" ~ scale_x_continuous(limits = c(0,30)),
    taxon == "Nitrospirota" ~ scale_x_continuous(limits = c(0,20)),
    taxon == "Basidiomycota" ~ scale_x_continuous(limits = c(0,35)),
    taxon == "Actinomycetota" ~ scale_x_continuous(limits = c(0,100)),
    taxon == "Ascomycota" ~ scale_x_continuous(limits = c(0,40)),
    taxon == "Gemmatimonadota" ~ scale_x_continuous(limits = c(0,5)),
    taxon == "Mucoromycota" ~ scale_x_continuous(limits = c(0,20)),
    taxon == "Chytridiomycota" ~ scale_x_continuous(limits = c(0,10)),
    taxon == "Zoopagomycota" ~ scale_x_continuous(limits = c(0,1)),
    taxon == "Pseudomonadota" ~ scale_x_continuous(limits = c(0,75)),
    taxon == "Nitrososphaerota" ~ scale_x_continuous(limits = c(0,3)),
    taxon == "Euryarchaeota" ~ scale_x_continuous(limits = c(0,1.5)),
    taxon == "Thermoproteota" ~ scale_x_continuous(limits = c(0,.1)),
    taxon == "Spirochaetota" ~ scale_x_continuous(#labels = label_number(accuracy = 1), 
                                                  limits = c(0,2.5)),
    taxon == "Bacillota" ~ scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0,6)),
    taxon == "Bacteroidota" ~ scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0,7)),
    taxon == "Verrucomicrobiota" ~ scale_x_continuous(limits = c(0,5))
),
    y = list(taxon == "Acidobacteriota" ~ scale_y_continuous(
limits = c(0,30)),
             taxon == "Nitrospirota" ~ scale_y_continuous(limits = c(0,20)),
             taxon == "Basidiomycota" ~ scale_y_continuous(limits = c(0,35)),
             taxon == "Actinomycetota" ~ scale_y_continuous(limits = c(0,100)),
             taxon == "Ascomycota" ~ scale_y_continuous(limits = c(0,40)),
             taxon == "Gemmatimonadota" ~ scale_y_continuous(limits = c(0,5)),
taxon == "Mucoromycota" ~ scale_y_continuous(limits = c(0,10)),
taxon == "Chytridiomycota" ~ scale_y_continuous(limits = c(0,10)),
taxon == "Zoopagomycota" ~ scale_y_continuous(limits = c(0,1)),
             taxon == "Pseudomonadota" ~ scale_y_continuous(limits = c(0,75)),
             taxon == "Nitrososphaerota" ~ scale_y_continuous(limits = c(0,3)),
taxon == "Euryarchaeota" ~ scale_y_continuous(limits = c(0,1.5)),
taxon == "Thermoproteota" ~ scale_y_continuous(limits = c(0,.1)),
taxon == "Spirochaetota" ~ scale_y_continuous(labels = label_number(accuracy = 1), 
                                                           limits = c(0,2.5)),
taxon == "Bacillota" ~ scale_y_continuous(labels = label_number(accuracy = 1), limits = c(0,6)),
taxon == "Bacteroidota" ~ scale_y_continuous(labels = label_number(accuracy = 1), limits = c(0,7)),
taxon == "Verrucomicrobiota" ~ scale_y_continuous(limits = c(0,5))
    )

)

a <- fig_pluspf_vs_smd_bac + fixed_scales + guides(color = "none") + theme(  panel.margin = unit(0, "lines"),
                                                                             strip.background = element_rect(fill ="white"))

b <- fig_gtdb_vs_smd_bac + fixed_scales  + guides(color = "none")  + theme(panel.margin = unit(0, "lines"),
                                                                           strip.background = element_rect(fill ="white"))

c <- fig_pluspf_vs_smd_arc + fixed_scales + guides(color = "none") + theme(panel.margin = unit(0, "lines"),
                                                                             strip.background = element_rect(fill ="white"))
d <- fig_gtdb_vs_smd_arc + fixed_scales  + guides(color = "none")  + theme(panel.margin = unit(0, "lines"),
                                                                           strip.background = element_rect(fill ="white"))

e <- fig_pluspf_vs_smd_fun + 
    fixed_scales + 
    guides(color = "none")  + 
    theme(panel.margin = unit(0, "lines"),
          strip.background = element_rect(fill ="white"))




plot_for_legend = ggplot(bracken_phylum_long %>%
                                 filter(taxon %in% also_plot)) +    
    geom_point(aes(x = soil_microbe_db, y = gtdb_207, color = Domain))  + theme_bw() + theme(legend.text=element_text(size=16),
                                                                                             legend.title=element_text(size=16))
legend_df <- get_legend(plot_for_legend)


legend_plot <- as_ggplot(legend_df) + theme(panel.margin = unit(0, "lines"), strip.background = element_rect(fill ="white")) 

library(patchwork)


fig_s5 = ((a | b) | 
        ((c | d) + plot_layout(axis_titles = "collect")) / 
        (e | legend_plot) + plot_layout(heights = c(3, 4))) + 
    plot_layout(axis_titles = "collect", widths = c(1,1,2.4)) + 
    #    plot_layout(guides = 'collect') +
    # plot_layout(heights = c(8,1, 2)) +
    plot_annotation(tag_levels = list(c('A',"B","C","D","E")),
                    theme = theme(plot.margin = unit(c(0,0,0,0), 
                                                     'lines')))

fig_s5 %>% 
    ggexport(height = 1000, width = 1200, filename = "/projectnb/frpmars/soil_microbe_db/manuscript_figures/fig_s5.png")
