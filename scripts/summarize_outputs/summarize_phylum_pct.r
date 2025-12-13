# Summarize phylum-level abundance from bracken files
# Uses phylum-specific merged CSV files from taxonomic_rank_summaries
# Creates: data/classification/taxonomic_rank_summaries/phylum/bracken_phylum_estimates.rds

library(tidyverse)
library(pavian)
library(data.table)
source("scripts/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")

# Find all phylum CSV files directly in taxonomic_rank_summaries (ignore subfolders)
summary_dir <- "data/classification/taxonomic_rank_summaries"
phylum_csv_files <- list.files(summary_dir, pattern = "phylum.*merged.*lineage\\.csv$", 
                                full.names = TRUE, recursive = FALSE)

if(length(phylum_csv_files) == 0) {
    stop("No phylum CSV files found in ", summary_dir, 
         "\nExpected pattern: *phylum*merged*lineage.csv")
}

phylum_output_list <- list()

# Process each phylum CSV file
for(phylum_csv_file in phylum_csv_files) {
    cat("✓ Processing phylum CSV file:", basename(phylum_csv_file), "\n")
    
    # Extract database name from filename
    # Pattern: {db_name}_filtered_phylum_merged_lineage.csv or {db_name}_phylum_merged_lineage.csv
    filename <- basename(phylum_csv_file)
    db_name <- str_remove(filename, "_filtered_phylum_merged_lineage\\.csv$")
    db_name <- str_remove(db_name, "_phylum_merged_lineage\\.csv$")
    db_name <- str_remove(db_name, "_phylum_merged\\.csv$")
    
    phylum_csv <- fread(phylum_csv_file, nThread = 8)
    
    # Convert CSV format to match expected output format
    phylum_data <- phylum_csv %>%
        filter(taxonomy_lvl == "P") %>%
        mutate(
            sampleID_temp = str_remove(sample_id, "_phylum_filtered"),
            taxon = name,
            percentage = fraction_total_reads * 100,
            samp_name = sample_id
        ) %>%
        separate(sampleID_temp, into = c("sampleID", "db_name_from_sample"), sep = "COMP_", remove = F, extra = "merge") %>%
        mutate(
            sampleID = str_remove(sampleID, "-$"),
            db_name_from_sample = str_remove(db_name_from_sample, "_phylum_filtered"),
            db_name = ifelse(is.na(db_name_from_sample) | db_name_from_sample == "", db_name, db_name_from_sample)
        ) %>%
        select(taxon, samp_name, sampleID, db_name, percentage)
    
    phylum_data$siteID <- substr(phylum_data$sampleID, 1, 4)
    phylum_output_list[[db_name]] <- phylum_data
}

# Combine all databases
if(length(phylum_output_list) > 0) {
    phylum_output <- bind_rows(phylum_output_list)
} else {
    stop("No phylum data found in CSV files!")
}

saveRDS(phylum_output,"data/classification/taxonomic_rank_summaries/phylum/bracken_phylum_estimates.rds")

cat("✅ Saved phylum estimates to: data/classification/taxonomic_rank_summaries/phylum/bracken_phylum_estimates.rds\n")
cat("   Databases included:", paste(unique(phylum_output$db_name), collapse=", "), "\n")
