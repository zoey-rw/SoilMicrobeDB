# Summarize taxonomic rank abundances (domain, phylum, genus) from merged CSV files
# Input: Merged CSV files with lineage from 03_add_lineage.sh
#        Expected in data/classification/taxonomic_rank_summaries/
#        Pattern: *{rank}*merged*lineage.csv
# Output: bracken_{rank}_estimates.rds files in data/classification/taxonomic_rank_summaries/{rank}/
#         Each RDS contains taxon, samp_name, sampleID, db_name, percentage, siteID

library(tidyverse)
library(data.table)

ranks <- c("domain", "phylum", "genus")
rank_to_taxonomy_lvl <- c("domain" = "D", "phylum" = "P", "genus" = "G")
rank_to_suffix <- c("domain" = "_domain_filtered", "phylum" = "_phylum_filtered", "genus" = "_genus_filtered")
summary_dir <- "data/classification/taxonomic_rank_summaries"

for(rank in ranks) {
    pattern <- paste0(rank, ".*merged.*lineage\\.csv$")
    csv_files <- list.files(summary_dir, pattern = pattern, full.names = TRUE, recursive = FALSE)
    
    if(length(csv_files) == 0) {
        stop("No ", rank, " CSV files found in ", summary_dir)
    }
    
    rank_output_list <- list()
    taxonomy_lvl_code <- rank_to_taxonomy_lvl[rank]
    suffix_to_remove <- rank_to_suffix[rank]
    
    for(file_path in csv_files) {
        filename <- basename(file_path)
        db_name <- str_remove(filename, paste0("_filtered_", rank, "_merged_lineage\\.csv$"))
        db_name <- str_remove(db_name, paste0("_", rank, "_merged_lineage\\.csv$"))
        db_name <- str_remove(db_name, paste0("_", rank, "_merged\\.csv$"))
        
        rank_data <- fread(file_path, nThread = 8) %>%
            filter(taxonomy_lvl == taxonomy_lvl_code) %>%
            mutate(
                sampleID_temp = str_remove(sample_id, suffix_to_remove),
                taxon = name,
                percentage = fraction_total_reads * 100,
                samp_name = sample_id
            ) %>%
            separate(sampleID_temp, into = c("sampleID", "db_name_from_sample"), 
                    sep = "COMP_", remove = FALSE, extra = "merge") %>%
            mutate(
                sampleID = str_remove(sampleID, "-$"),
                db_name_from_sample = str_remove(db_name_from_sample, suffix_to_remove),
                db_name = ifelse(is.na(db_name_from_sample) | db_name_from_sample == "", 
                               db_name, db_name_from_sample),
                siteID = substr(sampleID, 1, 4)
            ) %>%
            select(taxon, samp_name, sampleID, db_name, percentage, siteID)
        
        rank_output_list[[db_name]] <- rank_data
    }
    
    if(length(rank_output_list) == 0) stop("No ", rank, " data processed")
    
    output_dir <- file.path(summary_dir, rank)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(bind_rows(rank_output_list), file.path(output_dir, paste0("bracken_", rank, "_estimates.rds")))
}
