library(tidyverse)
library(data.table)
library(ggpubr)

seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds") %>%
    select(-any_of(c("db_name", "identified_reads"))) %>%
    rename(compositeSampleID = sampleID)

genus_bracken <- rbindlist(list(
    fread("data/classification/taxonomic_rank_summaries/genus/soil_microbe_db_filtered_genus_merged.csv") %>%
        mutate(db_name = "soil_microbe_db"),
    fread("data/classification/taxonomic_rank_summaries/genus/pluspf_filtered_genus_merged.csv") %>%
        mutate(db_name = "pluspf")
)) %>%
    mutate(compositeSampleID = str_remove(sample_id, paste0("_", db_name, "_genus_filtered")))

keep_samples <- intersect(
    genus_bracken %>% filter(db_name == "pluspf") %>% pull(compositeSampleID),
    genus_bracken %>% filter(db_name == "soil_microbe_db") %>% pull(compositeSampleID)
) %>% unique()

genus_df <- genus_bracken %>%
    filter(compositeSampleID %in% keep_samples) %>%
    mutate(siteID = substr(compositeSampleID, 1, 4)) %>%
    left_join(seq_depth_df)

ecto_genera <- c("Cenococcum", "Gigaspora", "Aspergillus", "Lactifluus", 
                 "Lactarius", "Russula", "Xylaria", "Trichoderma")

ecto <- genus_df %>%
    filter(name %in% ecto_genera, db_name %in% c("pluspf", "soil_microbe_db"))

ggplot(ecto, aes(x = reorder(name, -fraction_total_reads), y = fraction_total_reads, color = db_name)) +
    geom_point(alpha = 0.5, position = position_dodge(width = 0.75)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_sqrt() +
    scale_color_discrete(name = "Database") +
    stat_compare_means(hide.ns = TRUE, label = "p.signif", label.y.npc = 0.7, size = 7) +
    theme_bw() +
    labs(x = "Genus", y = "Abundance")

ecto_wide <- ecto %>%
    select(compositeSampleID, fraction_total_reads, db_name, name) %>%
    pivot_wider(values_from = fraction_total_reads, names_from = db_name) %>%
    filter(!is.na(soil_microbe_db) & !is.na(pluspf))

fig_pluspf_vs_smd_fun <- ggplot(ecto_wide, aes(x = soil_microbe_db, y = pluspf)) +
    geom_point(size = 3, alpha = 0.4) +
    geom_smooth(method = "lm", formula = y ~ x) +
    geom_abline(slope = 1, intercept = 0, color = 1, linetype = 2, linewidth = 1.2) +
    facet_wrap(~name, ncol = 1, scales = "free") +
    theme_bw(base_size = 16) +
    labs(x = "SoilMicrobeDB rel. abundance", y = "PlusPF rel. abundance", color = NULL)
