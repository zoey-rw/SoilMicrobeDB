# Summarize domain-level abundance from bracken kreport files
# Creates: data/classification/taxonomic_rank_summaries/domain/bracken_domain_estimates.rds

library(tidyverse)
library(pavian)
source("scripts/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")

soil_sample_dir <- "data/classification/02_bracken_output"

samp_files <- list.files(soil_sample_dir, recursive=T, pattern = "_filtered_kraken_bracken_domains.kreport", full.names = T)

# Check if bracken files exist
if(length(samp_files) == 0) {
  stop("No bracken domain kreport files found in ", soil_sample_dir, 
       ". Please ensure bracken output files are available.")
}

samp_names <- gsub("_filtered_kraken_bracken_domains.kreport","",basename(samp_files)) %>% unique()
names(samp_files) <- samp_names

domain_in = many_files_to_matrix_list(files = samp_files, filter.tax.level = "D", 
                                      include.unclassified = T, percentages = T)[[1]] %>%
    as.data.frame() %>%
    rownames_to_column("taxon") %>%
    pivot_longer(cols = where(is.numeric), names_to="samp_name",values_to="percentage")

domain_output <- domain_in %>% 
    separate(samp_name, 
             into = c("sampleID","db_name"), 
             sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(samp_name, paste0("_",db_name)))
domain_output$siteID = substr(domain_output$sampleID, 1, 4)

saveRDS(domain_output,"data/classification/taxonomic_rank_summaries/domain/bracken_domain_estimates.rds")

cat("✅ Saved domain estimates to: data/classification/taxonomic_rank_summaries/domain/bracken_domain_estimates.rds\n")
