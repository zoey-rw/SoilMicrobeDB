library(tidyverse)
library(data.table)
library(ggpmisc)
library(ggpubr)

source("scripts/helper_functions.r")

AM_trees <- read_csv("reference_data/AM_trees/AM_trees_lang.csv")
bracken_genus <- fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_genus_merged_lineage.csv", nThread = 8)
funguild_results_full <- fread("data/classification/analysis_files/FUNguild_assignments.csv")

fungal_phyla <- c("Ascomycota","Basidiomycota","Blastocladiomycota","Chytridiomycota",
                  "Cryptomycota","Mucoromycota","Microsporidia","Olpidiomycota","Zoopagomycota")

funguild_results_simple <- funguild_results_full %>%
    filter(grepl("Prob", confidenceRanking)) %>%
    select(name, guild) %>%
    mutate(Ectomycorrhizae = grepl("Ecto", guild),
           `Arbuscular mycorrhizae` = grepl("Arbusc", guild)) %>%
    select(-guild)

emf_list <- funguild_results_simple %>% filter(Ectomycorrhizae) %>% pull(name) %>% tolower()
amf_list <- funguild_results_simple %>% filter(`Arbuscular mycorrhizae`) %>% pull(name) %>% tolower()

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
           is_emf = tolower(Genus) %in% emf_list,
           is_amf = tolower(Genus) %in% amf_list)

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

calc_rel_abundance <- function(df, myco_type) {
    type_col <- paste0("is_", myco_type)
    df %>%
        group_by(sample_id, is_fungi) %>%
        mutate(fungi_abun = sum(fraction_total_reads)) %>%
        ungroup() %>%
        group_by(sample_id, is_fungi, fungi_abun, !!sym(type_col)) %>%
        summarize(abun = sum(fraction_total_reads), .groups = "drop") %>%
        filter(!!sym(type_col) & is_fungi) %>%
        mutate(abun_rel = abun/fungi_abun,
               db_name = "soil_microbe_db") %>%
        select(-c(!!sym(type_col), is_fungi)) %>%
        parse_sample_id() %>%
        add_neon_metadata()
}

emf_counts_rel <- calc_rel_abundance(bracken_genus, "emf") %>%
    rename(emf_abun_rel = abun_rel)
amf_counts_rel <- calc_rel_abundance(bracken_genus, "amf") %>%
    rename(amf_abun_rel = abun_rel)

emf_rel_EM_trees <- merge(emf_counts_rel, AM_trees)
amf_rel_AM_trees <- merge(amf_counts_rel, AM_trees)

p2 <- ggplot(emf_rel_EM_trees, aes(x = `Proportion ECM basal area`, y = emf_abun_rel)) +
    geom_jitter(alpha = 0.8, size = 2, width = 0.01, height = 0.01) +
    stat_smooth(method = "lm") +
    scale_y_log10() +
    stat_cor(label.y.npc = 1, size = 7, p.accuracy = 0.0001) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size = 7, label.y.npc = 0.95) +
    theme_bw(base_size = 16) +
    theme(plot.margin = unit(c(3, 25, 25, 45), "pt"),
          axis.title = element_text(face = "bold")) +
    labs(x = "EMF-associated tree basal area", y = "EMF rel.abundance in metagenome")

# Load environmental metadata with nlcdClass for forest type analysis
if(!exists("soilCores")) {
    soilCores <- load_soilCores()
}

# Add forest type classification based on nlcdClass
forest_metadata <- soilCores %>%
    select(compositeSampleID, nlcdClass) %>%
    distinct(compositeSampleID, .keep_all = TRUE) %>%
    mutate(forest_type = case_when(
        nlcdClass == "evergreenForest" ~ "Conifer",
        nlcdClass == "deciduousForest" ~ "Broadleaf",
        nlcdClass == "mixedForest" ~ "Mixed",
        nlcdClass == "woodyWetlands" ~ "Wetlands",
        TRUE ~ NA_character_
    ))

# Merge forest type with EMF data for forest type analysis
emf_rel_EM_trees_with_metadata <- emf_rel_EM_trees %>%
    left_join(forest_metadata, by = "compositeSampleID")

# Check which nlcdClass values exist and how many samples were excluded
nlcdClass_summary <- emf_rel_EM_trees_with_metadata %>%
    count(nlcdClass, forest_type, name = "n_samples") %>%
    arrange(desc(n_samples))

cat("\n=== nlcdClass distribution in EMF data ===\n")
print(nlcdClass_summary)
cat("\nSamples with forest type:", sum(!is.na(emf_rel_EM_trees_with_metadata$forest_type)), "\n")
cat("Samples without forest type:", sum(is.na(emf_rel_EM_trees_with_metadata$forest_type)), "\n")
cat("Excluded nlcdClass categories:", 
    paste(unique(emf_rel_EM_trees_with_metadata$nlcdClass[is.na(emf_rel_EM_trees_with_metadata$forest_type)]), 
          collapse = ", "), "\n\n")

emf_rel_EM_trees_forest <- emf_rel_EM_trees_with_metadata %>%
    filter(!is.na(forest_type))

# Create plot with forest type coloring
p2_forest_type <- ggplot(emf_rel_EM_trees_forest, 
                         aes(x = `Proportion ECM basal area`, y = emf_abun_rel, color = forest_type)) +
    geom_jitter(alpha = 0.8, size = 2, width = 0.01, height = 0.01) +
    stat_smooth(method = "lm", se = TRUE) +
    scale_y_log10() +
    theme_bw(base_size = 16) +
    theme(plot.margin = unit(c(3, 25, 25, 45), "pt"),
          axis.title = element_text(face = "bold"),
          legend.position = "right") +
    labs(x = "EMF-associated tree basal area", 
         y = "EMF rel.abundance in metagenome",
         color = "Forest Type")

# Create individual plots for each forest type
create_forest_plot <- function(data, forest_type_name, is_first = FALSE, is_last = FALSE) {
    data %>%
        filter(forest_type == forest_type_name) %>%
        ggplot(aes(x = `Proportion ECM basal area`, y = emf_abun_rel, color = Site)) +
        geom_jitter(alpha = 0.8, size = 2, width = 0.01, height = 0.01) +
        stat_smooth(method = "lm", color = "black", se = TRUE) +
        scale_y_log10() +
        stat_cor(label.y.npc = 1, size = 5, p.accuracy = 0.0001, color = "black") +
        stat_regline_equation(aes(label = ..rr.label..),
                              show.legend = FALSE, size = 5, label.y.npc = 0.95, color = "black") +
        theme_bw(base_size = 16) +
        theme(plot.margin = unit(c(3, ifelse(is_last, 25, 0), 25, ifelse(is_first, 45, 0)), "pt"),
              axis.title = element_text(face = "bold"),
              plot.title = element_text(face = "bold", size = 14),
              axis.title.x = element_text(margin = margin(t = 0))) +
        labs(x = "", 
             y = ifelse(forest_type_name == "Broadleaf", "EMF rel.abundance in metagenome", ""),
             title = forest_type_name)
}

# Get forest types in order
forest_types <- c("Broadleaf", "Conifer", "Mixed", "Wetlands")
p2_plots <- lapply(seq_along(forest_types), function(i) {
    create_forest_plot(emf_rel_EM_trees_forest, forest_types[i], 
                       is_first = (i == 1), 
                       is_last = (i == length(forest_types)))
})

# Combine plots with A-D labels, with reduced spacing
p2_faceted <- ggarrange(plotlist = p2_plots, 
                        nrow = 1, ncol = 4,
                        labels = c("A", "B", "C", "D"),
                        label.x = 0, label.y = 1,
                        font.label = list(size = 18, face = "bold"),
                        common.legend = TRUE, legend = "right",
                        widths = c(1, 1, 1, 1.15)) %>%
    annotate_figure(bottom = grid::textGrob("EMF-associated tree basal area", 
                                            gp = grid::gpar(fontface = "bold", fontsize = 16)))

ggsave("manuscript_figures/fig5c_forest_type.png", p2_forest_type, width = 12, height = 8, units = "in", dpi = 300)
ggsave("manuscript_figures/fig5c_forest_type_faceted.png", p2_faceted, width = 14, height = 5, units = "in", dpi = 300)

p1 <- ggplot(amf_rel_AM_trees, aes(x = `Proportion AM basal area`, y = amf_abun_rel)) +
    geom_jitter(alpha = 0.8, size = 2, width = 0.01, height = 0.01) +
    stat_smooth(method = "lm") +
    scale_y_log10() +
    stat_cor(label.y.npc = 1, size = 7, p.accuracy = 0.0001) +
    stat_regline_equation(aes(label = ..rr.label..),
                          show.legend = FALSE, size = 7, label.y.npc = 0.95) +
    theme_bw(base_size = 16) +
    theme(plot.margin = unit(c(3, 25, 25, 45), "pt"),
          axis.title = element_text(face = "bold")) +
    labs(x = "AMF-associated tree basal area", y = "AMF rel.abundance in metagenome")

ggsave("manuscript_figures/fig5b.png", p1, width = 10, height = 8, units = "in", dpi = 300)
ggsave("manuscript_figures/fig5c.png", p2, width = 10, height = 8, units = "in", dpi = 300)

if(!exists("fig_5a")) {
    source("scripts/visualize/compare_its.r")
}

if(exists("fig_5a")) {
    fig_5bc <- ggarrange(p1, p2, common.legend = TRUE, nrow = 1, labels = c("B","C"))
    fig5 <- ggarrange(fig_5a, fig_5bc, nrow = 2, labels = c("A","",""), heights = c(2,1))
    ggsave("manuscript_figures/fig5.png", fig5, width = 12, height = 14, units = "in", dpi = 300)
}
