# Process NEON ITS amplicon data to genus-level relative abundances
# Input: Phyloseq object from /projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_ITS.rds
#        Contains NEON legacy (pre-2015) and post-2015 ITS amplicon data
# Output: NEON_ITS_amplicon.csv with genus-level relative abundances
#         Filters samples with fewer than 3000 reads

library(phyloseq)
library(tidyverse)
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/binTaxGroups.r")

ps_its <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_ITS.rds")

out <- get_tax_level_abun(ps_its, tax_rank_list = "genus", min_seq_depth = 3000)
genus_rel_abun <- out$genus$rel.abundances %>% rownames_to_column("ITS_sampleID")

dir.create("data/comparison_data/its_amplicon", recursive = TRUE, showWarnings = FALSE)
write_csv(genus_rel_abun, "data/comparison_data/its_amplicon/NEON_ITS_amplicon.csv")
