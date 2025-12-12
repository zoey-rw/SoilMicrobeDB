library(tidyverse)

source("scripts/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")
source("scripts/helper_functions.r")


unfilter_reads_dir = "data/classification/01_kraken_output"
samp_files <- list.files(unfilter_reads_dir, recursive=T, pattern = "_soil_microbe_db_kraken.kreport", full.names = T)

seq_depth = lapply(samp_files, function(report_path){
	my_report <- fread_report(report_path) %>% as.data.frame()
	rownames(my_report)[rownames(my_report) == "r_root"] <- "-_root"
	my_report <- my_report[!duplicated(my_report$name),]
	row.names(my_report) <- my_report[["name"]]
	unidentified_reads <- my_report["u_unclassified","cladeReads"]
	identified_reads <- my_report["r_root","cladeReads"]
	number_of_raw_reads = unidentified_reads + identified_reads
	return(data.frame("number_of_raw_reads" = number_of_raw_reads, 
	            "identified_reads" = identified_reads))
})


seq_depth_out = rbindlist(seq_depth)

seq_depth_df = cbind.data.frame(samp_files,
                                "sampleID_orig" = gsub("_kraken.kreport","",basename(samp_files)), 
                                seq_depth_out) %>% 
    separate(sampleID_orig, 
             into = c("sampleID","db_name"), 
             sep = "COMP_", remove = F, extra = "merge") %>% 
    mutate(sampleID = str_remove(sampleID_orig, paste0("_",db_name)))
seq_depth_df = seq_depth_df %>% select(sampleID, db_name, seq_depth = number_of_raw_reads, identified_reads)

saveRDS(seq_depth_df,"data/classification/analysis_files/seq_depth_df.rds")
