library(tidyverse)
library(data.table)

source("scripts/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")
source("scripts/helper_functions.r")

# Set up paths - check both local and HARDDRIVE locations
# Before filtering: original Kraken2 kreport files
kreport_before_local = "data/NEON_metagenome_classification/01_kraken_output"
kreport_before_harddrive = "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/01_kraken_output"

# After filtering: filtered kreport files
kreport_after_local = "data/NEON_metagenome_classification/02_bracken_output"
kreport_after_harddrive = "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/02_bracken_output"

# Find all kreport files (all databases, both before and after filtering)
# Process both _kraken.kreport (before filtering) and _filtered_kraken.kreport (after filtering)
# Sequencing depth is the same regardless of filtering status
samp_files <- character(0)

# Before filtering files
if(dir.exists(kreport_before_local)) {
    samp_files <- c(samp_files, list.files(kreport_before_local, recursive=T, pattern = "_kraken.kreport$", full.names = T))
}
if(dir.exists(kreport_before_harddrive)) {
    before_harddrive <- list.files(kreport_before_harddrive, recursive=T, pattern = "_kraken.kreport$", full.names = T)
    samp_files <- c(samp_files, before_harddrive)
    cat("✓ Found", length(before_harddrive), "kreport files (before filtering) on HARDDRIVE\n")
}

# After filtering files
if(dir.exists(kreport_after_local)) {
    samp_files <- c(samp_files, list.files(kreport_after_local, recursive=T, pattern = "_filtered_kraken.kreport$", full.names = T))
}
if(dir.exists(kreport_after_harddrive)) {
    after_harddrive <- list.files(kreport_after_harddrive, recursive=T, pattern = "_filtered_kraken.kreport$", full.names = T)
    samp_files <- c(samp_files, after_harddrive)
    cat("✓ Found", length(after_harddrive), "kreport files (after filtering) on HARDDRIVE\n")
}

cat("Found", length(samp_files), "total kreport files\n")

if(length(samp_files) == 0) {
    stop("❌ No kreport files found!\n",
         "   Checked: ", kreport_dir_local, "\n",
         "   Checked: ", kreport_dir_harddrive, "\n",
         "   Expected files matching pattern '*_kraken.kreport'")
}

# Function to parse sample ID from kreport filename (same logic as 06_calculate_classification_pct_by_rank.r)
parse_sample_id <- function(filename) {
    is_filtered <- grepl("_filtered_kraken.kreport", filename)
    if(is_filtered) {
        samp_name <- gsub("_filtered_kraken.kreport", "", basename(filename))
    } else {
        samp_name <- gsub("_kraken.kreport", "", basename(filename))
    }
    
    # Parse sampleID and db_name
    parts <- strsplit(samp_name, "COMP_", fixed = TRUE)[[1]]
    if(length(parts) >= 2) {
        sampleID_temp <- parts[1]
        db_name <- parts[2]
        # Remove db_name suffix from sampleID
        sampleID <- sub(paste0("_", db_name), "", samp_name, fixed = TRUE)
        db_name <- gsub("_filtered|_kraken", "", db_name)
    } else {
        # Try alternative parsing for files without COMP_ separator
        if(grepl("soil_microbe_db", samp_name)) {
            sampleID <- gsub("_soil_microbe_db", "", samp_name)
            db_name <- "soil_microbe_db"
        } else {
            sampleID <- samp_name
            db_name <- "unknown"
        }
    }
    
    return(list(sampleID = sampleID, db_name = db_name, samp_name = samp_name))
}

# Process all kreport files
seq_depth = lapply(samp_files, function(report_path){
    tryCatch({
        my_report <- fread_report(report_path) %>% as.data.frame()
        rownames(my_report)[rownames(my_report) == "r_root"] <- "-_root"
        my_report <- my_report[!duplicated(my_report$name),]
        row.names(my_report) <- my_report[["name"]]
        unidentified_reads <- my_report["u_unclassified","cladeReads"]
        identified_reads <- my_report["r_root","cladeReads"]
        number_of_raw_reads = unidentified_reads + identified_reads
        
        # Parse sample ID and database name
        sample_info <- parse_sample_id(report_path)
        
        return(data.frame(
            "sampleID" = sample_info$sampleID,
            "db_name" = sample_info$db_name,
            "number_of_raw_reads" = number_of_raw_reads, 
            "identified_reads" = identified_reads
        ))
    }, error = function(e) {
        cat("    Error processing", basename(report_path), ":", e$message, "\n")
        return(NULL)
    })
})

# Remove NULL results
seq_depth <- Filter(Negate(is.null), seq_depth)

if(length(seq_depth) == 0) {
    stop("❌ No valid kreport files processed!\n",
         "   Check that kreport files are readable and contain data.")
}

seq_depth_out = rbindlist(seq_depth)

# Sequencing depth is sample-specific, not database-specific
seq_depth_df = seq_depth_out %>%
    select(sampleID, seq_depth = number_of_raw_reads) %>%
    distinct(sampleID, .keep_all = TRUE)

cat("✓ Processed", nrow(seq_depth_df), "unique samples\n")
cat("  Total seq_depth records:", nrow(seq_depth_df), "\n")

# Save to both locations for compatibility
output_file_new = "data/NEON_metagenome_classification/seq_depth_df.rds"
output_file_old = "data/classification/analysis_files/seq_depth_df.rds"

# Create output directory if needed
output_dir_new <- dirname(output_file_new)
if(!dir.exists(output_dir_new)) {
    dir.create(output_dir_new, recursive = TRUE)
    cat("Created output directory:", output_dir_new, "\n")
}

saveRDS(seq_depth_df, output_file_new)
cat("✅ Saved sequencing depth to:", output_file_new, "\n")

saveRDS(seq_depth_df, output_file_old)
cat("✅ Saved sequencing depth to:", output_file_old, "\n")
