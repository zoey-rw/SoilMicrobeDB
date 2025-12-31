library(tidyverse)
library(data.table)
library(ggpubr)

source("scripts/helper_functions.r")

AM_trees <- read_csv("reference_data/AM_trees/AM_trees_lang.csv")

fungal_phyla <- c("Ascomycota","Basidiomycota","Blastocladiomycota","Chytridiomycota",
                  "Cryptomycota","Mucoromycota","Microsporidia","Olpidiomycota","Zoopagomycota")

parse_sample_id <- function(df) {
    df %>%
        separate(sample_id, into = c("compositeSampleID","db_name"), 
                 sep = "COMP_", remove = FALSE, extra = "merge") %>%
        mutate(compositeSampleID = str_remove(sample_id, paste0("_",db_name)))
}

add_neon_metadata <- function(df) {
    neon_meta <- microbialForecast:::parseNEONsampleIDs(df$compositeSampleID) %>%
        select(siteID, plotID, dateID, horizon, plot_date)
    df %>%
        cbind(neon_meta) %>%
        mutate(Site = siteID,
               `Plot Number` = as.numeric(gsub("^0+|^00+", "", substr(plotID, 6, 8))))
}

bracken_genus <- fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_genus_merged_lineage.csv", nThread = 8)
bracken_phylum <- fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_phylum_merged_lineage.csv", nThread = 8)

lineage_df <- lapply(unique(bracken_genus$lineage),
                     function(x) {
                         phyloseq::parse_taxonomy_qiime(x) %>%
                             t() %>% as.data.frame()
                     }) %>%
    data.table::rbindlist() %>%
    cbind.data.frame(lineage = unique(bracken_genus$lineage), .)

bracken_genus <- bracken_genus %>%
    left_join(lineage_df) %>%
    mutate(is_fungi = Phylum %in% fungal_phyla,
           is_amf = grepl("Glomero", lineage))

amf_counts_rel <- bracken_genus %>%
    group_by(sample_id, is_fungi) %>%
    mutate(fungi_abun = sum(fraction_total_reads)) %>%
    ungroup() %>%
    group_by(sample_id, is_fungi, fungi_abun, is_amf) %>%
    summarize(amf_abun = sum(fraction_total_reads), .groups = "drop") %>%
    filter(is_amf & is_fungi) %>%
    mutate(amf_abun_rel = amf_abun/fungi_abun,
           db_name = "soil_microbe_db") %>%
    select(-c(is_amf, is_fungi)) %>%
    parse_sample_id() %>%
    add_neon_metadata()

amf_rel_AM_trees <- merge(amf_counts_rel, AM_trees)

p1 <- ggplot(amf_rel_AM_trees, 
             aes(x = `Proportion AM basal area`, y = amf_abun_rel)) + 
    geom_jitter(aes(color = Site), alpha = 0.8, size = 2) +
    facet_grid(~horizon) +
    stat_smooth(method = "lm") +
    scale_y_log10() +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size = 6, label.x.npc = 0.4) +
    theme_bw(base_size = 18) +
    ylab("Proportion AMF in metagenome") +
    xlab("Proportion AM tree basal area")

bracken_phylum_estimates <- bracken_phylum %>%
    filter(grepl("Mucoro", lineage)) %>%
    mutate(db_name = "soil_microbe_db") %>%
    parse_sample_id() %>%
    add_neon_metadata()

mucoro_counts <- bracken_phylum_estimates %>%
    group_by(siteID, Site, plotID, `Plot Number`, db_name, dateID,
             plot_date, horizon, compositeSampleID) %>%
    summarize(AMF_percent = sum(fraction_total_reads), .groups = "drop")

bracken_genus_estimates <- bracken_genus %>%
    filter(is_amf) %>%
    mutate(db_name = "soil_microbe_db") %>%
    parse_sample_id() %>%
    add_neon_metadata()

amf_counts <- bracken_genus_estimates %>%
    group_by(siteID, Site, plotID, `Plot Number`, db_name, dateID,
             plot_date, horizon, compositeSampleID) %>%
    summarize(AMF_percent = sum(fraction_total_reads), .groups = "drop")

amf_counts_genus <- bracken_genus_estimates %>%
    group_by(siteID, Site, plotID, `Plot Number`, db_name, dateID, name,
             plot_date, horizon, compositeSampleID) %>%
    summarize(AMF_percent = sum(fraction_total_reads), .groups = "drop")

amf_counts_AM_trees <- merge(amf_counts, AM_trees)
mucoro_counts_AM_trees <- merge(mucoro_counts, AM_trees)
genus_counts_AM_trees <- merge(amf_counts_genus, AM_trees)

p1_mucoro <- ggplot(mucoro_counts_AM_trees, 
                    aes(x = `Proportion AM basal area`, y = AMF_percent)) + 
    geom_point(aes(color = Site)) +
    facet_grid(~horizon) +
    stat_smooth(method = "lm") +
    scale_y_log10() +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size = 6, label.x.npc = 0.4)

p2 <- ggplot(amf_counts_AM_trees %>% filter(dateID > 201501), 
             aes(x = `Proportion AM basal area`, y = AMF_percent)) + 
    geom_point(aes(color = Site)) +
    stat_smooth(method = "lm") +
    scale_y_log10() +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size = 6, label.x.npc = 0.4)
