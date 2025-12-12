# Summarize phylum-level abundance from bracken kreport files
# Creates: data/classification/taxonomic_rank_summaries/phylum/bracken_phylum_estimates.rds

library(tidyverse)
library(pavian)
source("scripts/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")

soil_sample_dir <- "data/classification/02_bracken_output"

samp_files <- list.files(soil_sample_dir, recursive=T, pattern = "_filtered_kraken_bracken_genuses.kreport", full.names = T)

# Check if bracken files exist
if(length(samp_files) == 0) {
  stop("No bracken kreport files found in ", soil_sample_dir, 
       ". Please ensure bracken output files are available.")
}

samp_names <- gsub("_filtered_kraken_bracken_genuses.kreport","",basename(samp_files)) %>% unique()
names(samp_files) <- samp_names

phylum_in = many_files_to_matrix_list(files = samp_files, filter.tax.level = "P", 
                                      include.unclassified = T, percentages = T)[[1]] %>%
    as.data.frame() %>%
    rownames_to_column("taxon") %>%
    pivot_longer(cols = where(is.numeric), names_to="samp_name",values_to="percentage")

phylum_output <- phylum_in %>% 
    separate(samp_name, 
             into = c("sampleID","db_name"), 
             sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(samp_name, paste0("_",db_name)))
phylum_output$siteID = substr(phylum_output$sampleID, 1, 4)

saveRDS(phylum_output,"data/classification/taxonomic_rank_summaries/phylum/bracken_phylum_estimates.rds")

cat("✅ Saved phylum estimates to: data/classification/taxonomic_rank_summaries/phylum/bracken_phylum_estimates.rds\n")
