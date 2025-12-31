library(tidyverse)
library(ggpubr)

domain_output <- readRDS("data/classification/taxonomic_rank_summaries/domain/bracken_domain_estimates.rds")

common_samples <- domain_output %>%
    distinct(sampleID, db_name) %>%
    count(sampleID) %>%
    filter(n > 2) %>%
    pull(sampleID)

options(scipen = 999)

fig_s1 <- domain_output %>%
    mutate(db_name = recode(db_name,
                           "soil_microbe_db" = "SoilMicrobeDB",
                           "pluspf" = "PlusPF",
                           "gtdb_207" = "GTDB r207")) %>%
    filter(sampleID %in% common_samples,
           taxon %in% c("Archaea", "Bacteria", "Eukaryota")) %>%
    ggplot(aes(x = db_name, y = percentage, color = db_name)) +
    geom_violin(color = 1, draw_quantiles = 0.5) +
    geom_jitter(size = 3, alpha = 0.2) +
    facet_wrap(~taxon, scale = "free_y") +
    stat_compare_means(method = "t.test",
                      comparisons = list(c("SoilMicrobeDB", "PlusPF"),
                                       c("SoilMicrobeDB", "GTDB r207"),
                                       c("GTDB r207", "PlusPF")),
                      hide.ns = TRUE) +
    theme_minimal(base_size = 20) +
    theme(axis.text.x = element_text(angle = 310, vjust = 1, hjust = 0)) +
    guides(color = "none") +
    labs(x = "Database", y = "Estimated abundances")

ggsave("manuscript_figures/fig_s1.png", fig_s1, width = 14, height = 5, units = "in", dpi = 300)


