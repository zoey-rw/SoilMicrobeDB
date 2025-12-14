# Calculate sequencing depth from Kraken2 kreport files
# Input: Kraken2 kreport files from 01_kraken_output and 02_bracken_output directories
#        Checks both local and HARDDRIVE paths
# Output: seq_depth_df.rds with sampleID and seq_depth columns
#         Saved to both data/NEON_metagenome_classification/ and data/classification/analysis_files/

library(tidyverse)
library(data.table)
source("scripts/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")
source("scripts/helper_functions.r")

parse_sample_id <- function(filename) {
    is_filtered <- grepl("_filtered_kraken.kreport", filename)
    samp_name <- ifelse(is_filtered, 
                       gsub("_filtered_kraken.kreport", "", basename(filename)),
                       gsub("_kraken.kreport", "", basename(filename)))
    
    parts <- strsplit(samp_name, "COMP_", fixed = TRUE)[[1]]
    if(length(parts) >= 2) {
        sampleID <- sub(paste0("_", parts[2]), "", samp_name, fixed = TRUE)
        db_name <- gsub("_filtered|_kraken", "", parts[2])
    } else {
        sampleID <- ifelse(grepl("soil_microbe_db", samp_name),
                          gsub("_soil_microbe_db", "", samp_name),
                          samp_name)
        db_name <- ifelse(grepl("soil_microbe_db", samp_name), "soil_microbe_db", "unknown")
    }
    return(list(sampleID = sampleID, db_name = db_name, samp_name = samp_name))
}

kreport_dirs <- c("data/NEON_metagenome_classification/01_kraken_output",
                  "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/01_kraken_output",
                  "data/NEON_metagenome_classification/02_bracken_output",
                  "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/02_bracken_output")

samp_files <- unlist(lapply(kreport_dirs[dir.exists(kreport_dirs)], function(dir) {
    patterns <- ifelse(grepl("01_kraken_output", dir), "_kraken.kreport$", "_filtered_kraken.kreport$")
    list.files(dir, recursive = TRUE, pattern = patterns, full.names = TRUE)
}))

if(length(samp_files) == 0) {
    stop("No kreport files found in: ", paste(kreport_dirs, collapse = ", "))
}

seq_depth <- lapply(samp_files, function(report_path) {
    tryCatch({
        my_report <- fread_report(report_path) %>% as.data.frame()
        rownames(my_report)[rownames(my_report) == "r_root"] <- "-_root"
        my_report <- my_report[!duplicated(my_report$name),]
        row.names(my_report) <- my_report[["name"]]
        
        sample_info <- parse_sample_id(report_path)
        data.frame(
            sampleID = sample_info$sampleID,
            db_name = sample_info$db_name,
            number_of_raw_reads = my_report["u_unclassified","cladeReads"] + my_report["r_root","cladeReads"],
            identified_reads = my_report["r_root","cladeReads"]
        )
    }, error = function(e) NULL)
})

seq_depth_df <- rbindlist(Filter(Negate(is.null), seq_depth)) %>%
    select(sampleID, seq_depth = number_of_raw_reads) %>%
    distinct(sampleID, .keep_all = TRUE)

dir.create(dirname("data/NEON_metagenome_classification/seq_depth_df.rds"), recursive = TRUE, showWarnings = FALSE)
saveRDS(seq_depth_df, "data/NEON_metagenome_classification/seq_depth_df.rds")
saveRDS(seq_depth_df, "data/classification/analysis_files/seq_depth_df.rds")
