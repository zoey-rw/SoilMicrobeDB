# Process NEON ITS amplicon data to genus-level relative abundances
# Input: Phyloseq object from /projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_ITS.rds
#        Contains NEON legacy (pre-2015) and post-2015 ITS amplicon data
# Output: NEON_ITS_amplicon.csv with genus-level relative abundances
#         Filters samples with fewer than 3000 reads

library(phyloseq)
library(tidyverse)
source("scripts/helper_functions.r")

ps_its_paths <- c(
    "data/comparison_data/its_amplicon/phyloseq_ITS.rds",
    "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_ITS.rds"
)

ps_its_file <- ps_its_paths[file.exists(ps_its_paths)][1]
if(is.na(ps_its_file)) {
    stop("phyloseq_ITS.rds not found in:\n  ", paste(ps_its_paths, collapse = "\n  "))
}

ps_its <- readRDS(ps_its_file)

out <- get_tax_level_abun(ps_its, tax_rank_list = "genus", min_seq_depth = 3000)
genus_rel_abun <- out$genus$rel.abundances %>% rownames_to_column("ITS_sampleID")

dir.create("data/comparison_data/its_amplicon", recursive = TRUE, showWarnings = FALSE)
write_csv(genus_rel_abun, "data/comparison_data/its_amplicon/NEON_ITS_amplicon.csv")
