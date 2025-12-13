# Summarize phylum-level abundance from bracken files
# Uses phylum-specific merged CSV for soil_microbe_db (more accurate)
# Falls back to kreport files for other databases
# Creates: data/classification/taxonomic_rank_summaries/phylum/bracken_phylum_estimates.rds

library(tidyverse)
library(pavian)
library(data.table)
source("scripts/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")

# Check if phylum CSV exists for soil_microbe_db (created from phylum-specific .b2 files)
phylum_csv_file <- "data/classification/taxonomic_rank_summaries/phylum/soil_microbe_db_filtered_phylum_merged_lineage.csv"

phylum_output_list <- list()

# Use CSV file for soil_microbe_db if it exists
if(file.exists(phylum_csv_file)) {
    cat("✓ Using phylum CSV file for soil_microbe_db:", phylum_csv_file, "\n")
    phylum_csv <- fread(phylum_csv_file, nThread = 8)
    
    # Convert CSV format to match expected output format
    phylum_smd <- phylum_csv %>%
        filter(taxonomy_lvl == "P") %>%
        mutate(
            sampleID = str_remove(sample_id, "_phylum_filtered"),
            db_name = "soil_microbe_db",
            taxon = name,
            percentage = fraction_total_reads * 100,
            samp_name = sample_id
        ) %>%
        separate(sampleID, into = c("sampleID", "db_name"), sep = "COMP_", remove = F, extra = "merge") %>%
        mutate(sampleID = str_remove(samp_name, paste0("_", db_name))) %>%
        select(taxon, samp_name, sampleID, db_name, percentage)
    
    phylum_smd$siteID <- substr(phylum_smd$sampleID, 1, 4)
    phylum_output_list[["soil_microbe_db"]] <- phylum_smd
}

# For other databases, use kreport files
soil_sample_dir <- "data/classification/02_bracken_output"
samp_files <- list.files(soil_sample_dir, recursive=T, pattern = "_filtered_kraken_bracken_genuses.kreport", full.names = T)

if(length(samp_files) > 0) {
    # Filter to exclude soil_microbe_db samples (already processed from CSV)
    samp_files <- samp_files[!grepl("soil_microbe_db", samp_files)]
    
    if(length(samp_files) > 0) {
        cat("✓ Using kreport files for other databases\n")
        samp_names <- gsub("_filtered_kraken_bracken_genuses.kreport","",basename(samp_files)) %>% unique()
        names(samp_files) <- samp_names
        
        phylum_in <- many_files_to_matrix_list(files = samp_files, filter.tax.level = "P", 
                                              include.unclassified = T, percentages = T)[[1]] %>%
            as.data.frame() %>%
            rownames_to_column("taxon") %>%
            pivot_longer(cols = where(is.numeric), names_to="samp_name",values_to="percentage")
        
        phylum_other <- phylum_in %>% 
            separate(samp_name, 
                     into = c("sampleID","db_name"), 
                     sep = "COMP_", remove = F, extra = "merge") %>% 
            mutate(sampleID = str_remove(samp_name, paste0("_",db_name)))
        phylum_other$siteID = substr(phylum_other$sampleID, 1, 4)
        
        # Combine with soil_microbe_db data
        phylum_output_list[["other"]] <- phylum_other
    }
}

# Combine all databases
if(length(phylum_output_list) > 0) {
    phylum_output <- bind_rows(phylum_output_list)
} else {
    stop("No phylum data found! Check that either:\n",
         "  1. CSV file exists: ", phylum_csv_file, "\n",
         "  2. OR kreport files exist in: ", soil_sample_dir)
}

saveRDS(phylum_output,"data/classification/taxonomic_rank_summaries/phylum/bracken_phylum_estimates.rds")

cat("✅ Saved phylum estimates to: data/classification/taxonomic_rank_summaries/phylum/bracken_phylum_estimates.rds\n")
cat("   Databases included:", paste(unique(phylum_output$db_name), collapse=", "), "\n")
