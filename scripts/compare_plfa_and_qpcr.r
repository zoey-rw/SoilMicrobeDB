library(neonUtilities)
library(tidyverse)
library(broom)
library(ggpubr)
library(SimplyAgree)


# Read in processed PLFA data - generated using code in https://github.com/zoey-rw/SoilBiomassNEON
master_df <- readRDS("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/plfa_comparison.rds")

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
              mean_plfa_f_proportion = mean(proportion_fungi3))


a = plfa_qpcr_df_mean %>% #filter(f_proportion < .4) %>% 
    ggplot(aes(x = genus_fungi_metagenome, y=mean_qpcr_f_proportion)) + 
    geom_point(aes(color=biome),
               alpha=.5) +
    
               stat_smooth(method="lm") +
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
                         y =mean_plfa_f_proportion)) + 
    geom_point(aes(color=biome),
        alpha=.5) +
    stat_smooth(method="lm") +
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


