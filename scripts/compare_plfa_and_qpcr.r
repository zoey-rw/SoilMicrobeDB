library(neonUtilities)
library(tidyverse)
library(broom)
library(ggpubr)
library(SimplyAgree)
library(ggstatsplot)
library(patchwork)

options(scipen=999)
# Read in processed PLFA data - generated using code in https://github.com/zoey-rw/SoilBiomassNEON
master_df <- fread("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/plfa_comparison.csv")
#master_df <- readRDS("/projectnb/talbot-lab-data/zrwerbin/SoilBiomassNEON/NEON_microbial_biomass_PLFA.rds")


# Read in processed qPCR data - missing lab QC data
qpcr_wide = read_csv("/projectnb/frpmars/soil_microbe_db/ref_data/NEON_qpcr.csv") %>% 
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
    ggexport(height = 600, width = 1200, filename = "/projectnb/frpmars/soil_microbe_db/manuscript_figures/fig4.png")


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
                                     y =proportion_fungi3, color=biome)) + 
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


