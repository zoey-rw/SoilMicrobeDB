
library(phyloseq)
library(tidyverse)

# Function to bin the ASV abundances by taxonomy
#source("https://raw.githubusercontent.com/zoey-rw/neonSoilMicrobes/refs/heads/master/binTaxGroups.r")
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/binTaxGroups.r")

# Read in Phyloseq object with NEON legacy (pre-2015) and post-2015 ITS amplicon data
ps_its <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_ITS.rds")

# Summarize by genus, filter samples with fewer than 3000 reads
out <- get_tax_level_abun(ps_its, 
                          tax_rank_list = "genus", 
                          min_seq_depth = 3000)

genus_rel_abun <- out$genus$rel.abundances %>% rownames_to_column("ITS_sampleID") 
write_csv(genus_rel_abun, "data/comparison_data/its_amplicon/NEON_ITS_amplicon.csv")
