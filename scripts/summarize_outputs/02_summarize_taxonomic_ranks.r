# Summarize taxonomic rank abundances (domain, phylum, genus) from merged CSV files
# Uses merged CSV files with lineage from the classification pipeline
# Creates: bracken_domain_estimates.rds, bracken_phylum_estimates.rds, bracken_genus_estimates.rds

library(tidyverse)
library(data.table)

# Define ranks to process
ranks <- c("domain", "phylum", "genus")

# Map rank names to taxonomy level codes
rank_to_taxonomy_lvl <- c(
    "domain" = "D",
    "phylum" = "P",
    "genus" = "G"
)

# Map rank names to sample_id suffixes to remove
rank_to_suffix <- c(
    "domain" = "_domain_filtered",
    "phylum" = "_phylum_filtered",
    "genus" = "_genus_filtered"
)

summary_dir <- "data/classification/taxonomic_rank_summaries"

# Process each rank
for(rank in ranks) {
    cat("\n=== Processing", rank, "rank ===\n")
    
    # Find all CSV files for this rank
    pattern <- paste0(rank, ".*merged.*lineage\\.csv$")
    csv_files <- list.files(summary_dir, pattern = pattern, 
                            full.names = TRUE, recursive = FALSE)
    
    if(length(csv_files) == 0) {
        stop("❌ MISSING FILES: No ", rank, " CSV files found in ", summary_dir, 
             "\n   Expected pattern: *", rank, "*merged*lineage.csv",
             "\n   Please ensure ", rank, " merged CSV files exist from 03_add_lineage.sh")
    }
    
    cat("Found", length(csv_files), rank, "CSV file(s)\n")
    
    rank_output_list <- list()
    taxonomy_lvl_code <- rank_to_taxonomy_lvl[rank]
    suffix_to_remove <- rank_to_suffix[rank]
    
    # Process each CSV file
    for(file_path in csv_files) {
        cat("  ✓ Processing:", basename(file_path), "\n")
        
        # Extract database name from filename
        filename <- basename(file_path)
        db_name <- str_remove(filename, paste0("_filtered_", rank, "_merged_lineage\\.csv$"))
        db_name <- str_remove(db_name, paste0("_", rank, "_merged_lineage\\.csv$"))
        db_name <- str_remove(db_name, paste0("_", rank, "_merged\\.csv$"))
        
        rank_csv <- fread(file_path, nThread = 8)
        
        # Convert CSV format to match expected output format
        rank_data <- rank_csv %>%
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
                               db_name, db_name_from_sample)
            ) %>%
            select(taxon, samp_name, sampleID, db_name, percentage)
        
        rank_data$siteID <- substr(rank_data$sampleID, 1, 4)
        rank_output_list[[db_name]] <- rank_data
    }
    
    # Combine all databases for this rank
    if(length(rank_output_list) == 0) {
        stop("❌ ERROR: No ", rank, " data processed from CSV files!")
    }
    
    rank_output <- bind_rows(rank_output_list)
    
    # Create output directory if needed
    output_dir <- file.path(summary_dir, rank)
    if(!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Save output
    output_file <- file.path(output_dir, paste0("bracken_", rank, "_estimates.rds"))
    saveRDS(rank_output, output_file)
    
    cat("✅ Saved", rank, "estimates to:", output_file, "\n")
    cat("   Databases included:", paste(unique(rank_output$db_name), collapse=", "), "\n")
    cat("   Total records:", nrow(rank_output), "\n")
}

cat("\n=== Summary complete ===\n")
