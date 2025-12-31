library(neonUtilities)
library(tidyverse)
library(data.table)
library(broom)
library(ggpubr)
library(SimplyAgree)
library(ggstatsplot)
library(patchwork)

options(scipen=999)
# Read in processed PLFA data - generated using code in https://github.com/zoey-rw/SoilBiomassNEON
master_df <- fread("data/classification/analysis_files/plfa_comparison.csv")
#master_df <- readRDS("/projectnb/talbot-lab-data/zrwerbin/SoilBiomassNEON/NEON_microbial_biomass_PLFA.rds")


# Read in processed qPCR data - missing lab QC data
qpcr_wide = read_csv("data/comparison_data/qpcr/NEON_qpcr.csv") %>% 
    mutate(fungi = meanCopyNumber_fungi,
           bacarc = `meanCopyNumber_bacteria and archaea`) %>% 
    filter(fungi < 40000 & fungi > 0) %>% 
    filter(bacarc < 5000000000 & bacarc > 0)  %>% 
    filter(!(`meanCopyNumber_bacteria and archaea` < 100 & meanCopyNumber_fungi < 100) ) 

qpcr_to_merge = qpcr_wide  %>% 
    select(-c(sampleID, dataQF.x, dataQF.y)) %>% 
    filter(!is.na(compositeSampleID))


plfa_qpcr_df = left_join(master_df, qpcr_to_merge, 
                         relationship = "many-to-many")

# Average by composite sample to match with metagenome
plfa_qpcr_df_mean = plfa_qpcr_df %>% 
    group_by(compositeSampleID,siteID, biome, horizon, genus_fungi_metagenome) %>% 
    summarize(mean_qpcr_f_proportion = mean(f_proportion, na.rm=T),
              mean_qpcr_f_abun = mean(fungi, na.rm=T),
              mean_plfa_f_proportion = mean(proportion_fungi))


plfa_qpcr_df_long = plfa_qpcr_df %>% 
    rename("qPCR" = f_proportion, "PLFA" = proportion_fungi) %>% 
    pivot_longer(cols=c(qPCR, PLFA), names_to = "metric", values_to = "Fungal abundance")

# Add biome from soilCores if missing (for faceted figures)
if(sum(!is.na(plfa_qpcr_df_long$biome)) == 0) {
    source("scripts/helper_functions.r")
    if(!exists("soilCores")) {
        soilCores <- load_soilCores()
    }
    soilCores_subset = soilCores %>% 
        select(compositeSampleID, nlcdClass) %>% 
        distinct(compositeSampleID, .keep_all = T) %>%
        mutate(biome = recode(nlcdClass, 
                              "mixedForest" = "Forest",
                              "evergreenForest" = "Forest",
                              "deciduousForest" = "Forest",
                              "emergentHerbaceousWetlands" = "Wetlands",
                              "woodyWetlands" = "Wetlands",
                              "dwarfScrub" = "Shrubland",
                              "shrubScrub" = "Shrubland",
                              "sedgeHerbaceous" = "Herbaceous",
                              "grasslandHerbaceous" = "Herbaceous",
                              "pastureHay" = "Cultivated",
                              "cultivatedCrops" = "Cultivated"))
    
    plfa_qpcr_df_long = plfa_qpcr_df_long %>% 
        left_join(soilCores_subset %>% select(compositeSampleID, biome), 
                 by = "compositeSampleID")
}


fig_4ab <- grouped_ggscatterstats(
    plfa_qpcr_df_long,
    x = genus_fungi_metagenome,
    y = `Fungal abundance`,
    xlab = "Fungal relative abundance in metagenome",
    ylab = "Fungal relative abundance",
    results.subtitle = F, , marginal = FALSE,
    grouping.var = metric,
    ggtheme = theme_bw(base_size = 18), 
    #xsidehistogram.args = list( na.rm = TRUE),
    #ysidehistogram.args = list( na.rm = F),
    plotgrid.args = list(widths = c(.5,.5)),
    ggplot.component = list(scale_y_log10(),
                            scale_x_log10(),
                            stat_cor(label.y.npc = .97, size=7, p.accuracy = .0001,
                                     position = position_nudge(x = 0, y = .3)),
                            stat_regline_equation(aes(label = ..rr.label..),
                                                  show.legend = FALSE, size=7, label.y.npc = .9, 
                                                  position = position_nudge(x = 0, y = .3)),
                            theme( plot.margin = unit(c(1.5, 1.5, 1.5, 45), "pt"),
                                   axis.title = element_text(face = "bold")))
) + plot_annotation(tag_levels ="A")
#facet_wrap(~metric) + 
fig_4ab




fig_4ab %>% 
    ggexport(height = 600, width = 1200, filename = "manuscript_figures/fig4.png")


# Figure 4C: qPCR vs PLFA direct comparison
plfa_qpcr_df_filtered = plfa_qpcr_df %>% 
    filter(!is.na(f_proportion) & !is.na(proportion_fungi)) %>%
    filter(f_proportion > 0 & proportion_fungi > 0)

fig_4c <- plfa_qpcr_df_filtered %>%
    ggplot(aes(x = f_proportion, y = proportion_fungi)) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "lm", se = FALSE) +
    stat_cor(label.y.npc = .97, size = 7, p.accuracy = .0001,
             position = position_nudge(x = 0, y = .3)) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size = 7, label.y.npc = .9, 
                          position = position_nudge(x = 0, y = .3)) +
    scale_x_log10() +
    scale_y_log10() +
    xlab("Fungal relative abundance (qPCR)") +
    ylab("Fungal relative abundance (PLFA)") +
    theme_bw(base_size = 18) +
    theme(plot.margin = unit(c(1.5, 1.5, 1.5, 45), "pt"),
          axis.title = element_text(face = "bold"))

fig_4c

fig_4c %>% 
    ggexport(height = 600, width = 600, filename = "manuscript_figures/fig4c.png")


# Figure 4A faceted by biome (PLFA)
plfa_qpcr_df_4a <- plfa_qpcr_df_long %>%
    filter(metric == "PLFA") %>%
    filter(!is.na(`Fungal abundance`) & !is.na(genus_fungi_metagenome) & !is.na(biome)) %>%
    filter(`Fungal abundance` > 0)

if(nrow(plfa_qpcr_df_4a) > 0 && length(unique(plfa_qpcr_df_4a$biome[!is.na(plfa_qpcr_df_4a$biome)])) > 0) {
    fig_4a_biome <- plfa_qpcr_df_4a %>%
        ggplot(aes(x = genus_fungi_metagenome, y = `Fungal abundance`)) +
        geom_point(alpha = 0.5) +
        stat_smooth(method = "lm", se = FALSE) +
        stat_cor(label.y.npc = .97, size = 5, p.accuracy = .0001,
                 position = position_nudge(x = 0, y = .3)) +
        stat_regline_equation(aes(label = ..rr.label..),
                              show.legend = FALSE, size = 5, label.y.npc = .9, 
                              position = position_nudge(x = 0, y = .3)) +
        scale_x_log10() +
        scale_y_log10() +
        xlab("Fungal relative abundance in metagenome") +
        ylab("Fungal relative abundance (PLFA)") +
        facet_wrap(~biome) +
        theme_bw(base_size = 16) +
        theme(plot.margin = unit(c(1.5, 1.5, 1.5, 45), "pt"),
              axis.title = element_text(face = "bold"))
    
    fig_4a_biome
    
    fig_4a_biome %>% 
        ggexport(height = 800, width = 1200, filename = "manuscript_figures/fig4a_biome.png")
}


# Figure 4B faceted by biome (qPCR)
plfa_qpcr_df_4b <- plfa_qpcr_df_long %>%
    filter(metric == "qPCR") %>%
    filter(!is.na(`Fungal abundance`) & !is.na(genus_fungi_metagenome) & !is.na(biome)) %>%
    filter(`Fungal abundance` > 0)

if(nrow(plfa_qpcr_df_4b) > 0 && length(unique(plfa_qpcr_df_4b$biome[!is.na(plfa_qpcr_df_4b$biome)])) > 0) {
    fig_4b_biome <- plfa_qpcr_df_4b %>%
        ggplot(aes(x = genus_fungi_metagenome, y = `Fungal abundance`)) +
        geom_point(alpha = 0.5) +
        stat_smooth(method = "lm", se = FALSE) +
        stat_cor(label.y.npc = .97, size = 5, p.accuracy = .0001,
                 position = position_nudge(x = 0, y = .3)) +
        stat_regline_equation(aes(label = ..rr.label..),
                              show.legend = FALSE, size = 5, label.y.npc = .9, 
                              position = position_nudge(x = 0, y = .3)) +
        scale_x_log10() +
        scale_y_log10() +
        xlab("Fungal relative abundance in metagenome") +
        ylab("Fungal relative abundance (qPCR)") +
        facet_wrap(~biome) +
        theme_bw(base_size = 16) +
        theme(plot.margin = unit(c(1.5, 1.5, 1.5, 45), "pt"),
              axis.title = element_text(face = "bold"))
    
    fig_4b_biome
    
    fig_4b_biome %>% 
        ggexport(height = 800, width = 1200, filename = "manuscript_figures/fig4b_biome.png")
}


## OLD/TESTING


a = plfa_qpcr_df_mean %>% #filter(f_proportion < .4) %>% 
    ggplot(aes(x = genus_fungi_metagenome, y=mean_qpcr_f_proportion,color=biome)) + 
    geom_point(aes(color=biome),
               alpha=.5) +
    
               stat_smooth(method="lm", se=F) +
    # stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=7) + 
    scale_x_log10() + 
    scale_y_log10() + 
    xlab("Fungal abundance in metagenome") +
    ylab("Fungal abundance (qPCR)") +
    #facet_grid(~horizon, scales="free") + 
    theme_bw(base_size = 18) + 
    geom_abline(slope=1, intercept = 0) 




b = plfa_qpcr_df_mean %>% ggplot(aes(x = genus_fungi_metagenome, 
                         y =mean_plfa_f_proportion, color=biome)) + 
    geom_point(aes(color=biome),
        alpha=.5) +
    stat_smooth(method="lm", se=F) +
    # stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=7) + 
    scale_x_log10() + 
    scale_y_log10() + 
    # scale_x_sqrt() + 
    # scale_y_sqrt() + 
    
    xlab("Fungal abundance in metagenome") +
    ylab("Fungal abundance (PLFA)") + theme_bw(base_size = 18) +
    geom_abline(slope=1, intercept = 0) #+ 
   # facet_grid(~biome, scales="free") 

ggarrange(b, a, nrow=1, common.legend = T, labels = c("A","B"))


plfa_qpcr_df_long %>% #filter(f_proportion < .4) %>% 
    ggplot(aes(x = genus_fungi_metagenome, y=`Fungal abundance`,color=biome)) + 
    geom_point(aes(color=biome),
               alpha=.5) +
    
    stat_smooth(method="lm", se=F) +
    # stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=7) + 
    scale_x_log10() + 
    scale_y_log10() + 
    xlab("Fungal abundance in metagenome") +
    ylab("Fungal abundance (qPCR)") +
    facet_grid(metric~biome, scales="free") + 
    theme_bw(base_size = 18) + 
    geom_abline(slope=1, intercept = 0) 

    #scale_y_log10() + 
  #  scale_x_log10() +

ggscatterstats(
    plfa_qpcr_df,
    x = genus_fungi_metagenome,
    y = f_proportion) +
    ggplot2::geom_rug(sides = "b") + 
    scale_x_log10() +
    scale_y_log10() +
    
    stat_smooth(method="lm", se=F) +
    # stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=7) 
    #facet_wrap(~metric)

a1 = plfa_qpcr_df %>% #filter(f_proportion < .4) %>% 
    ggplot(aes(x = genus_fungi_metagenome, y=f_proportion,color=biome)) + 
    geom_point(aes(color=biome),
               alpha=.5) +
    
    stat_smooth(method="lm", se=F) +
    # stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=7) + 
    scale_x_log10() + 
    scale_y_log10() + 
    xlab("Fungal abundance in metagenome") +
    ylab("Fungal abundance (qPCR)") +
    facet_grid(~biome, scales="free") + 
    theme_bw(base_size = 18) + 
    geom_abline(slope=1, intercept = 0) 




b1 = plfa_qpcr_df %>% ggplot(aes(x = genus_fungi_metagenome, 
                                     y =proportion_fungi, color=biome)) + 
    geom_point(aes(color=biome),
               alpha=.5) +
    stat_smooth(method="lm", se=F) +
    # stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size=7) + 
    scale_x_log10() + 
    scale_y_log10() + 
    # scale_x_sqrt() + 
    # scale_y_sqrt() + 
    
    xlab("Fungal abundance in metagenome") +
    ylab("Fungal abundance (PLFA)") + theme_bw(base_size = 18) +
    geom_abline(slope=1, intercept = 0) + facet_grid(~biome, scales="free") 


ggarrange(b1, a1, nrow=2, common.legend = T, labels = c("A","B"))


