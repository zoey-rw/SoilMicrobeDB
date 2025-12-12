# Create list of high-abundance taxa, from mapping NEON metagenomes to SoilGenome DB

library(tidyverse)
library(pavian)
library(ggpubr)
source("scripts/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")


sampleID_pooled_key <- readRDS("data/comparison_data/its_amplicon/NEON_pooled_sample_key.rds")
pooled_recode_list <- sampleID_pooled_key$genomicsSampleID
names(pooled_recode_list) <- sampleID_pooled_key$sampleID

soil_sample_dir <- "data/classification/02_bracken_output"

samp_files <- list.files(soil_sample_dir, recursive=T, pattern = "_filtered_kraken_bracken_genuses.kreport", full.names = T)
samp_names <- gsub("_filtered_kraken_bracken_genuses.kreport","",basename(samp_files)) %>% unique()
names(samp_files) <- samp_names

samp_files <- samp_files[grepl("HARV",samp_files)]

genus_in = many_files_to_matrix_list(files = samp_files, filter.tax.level = "G", 
																			include.unclassified = T, percentages = T)[[1]] %>%
	as.data.frame() %>%
	rownames_to_column("taxon") %>%
	pivot_longer(cols = where(is.numeric), names_to="samp_name",values_to="percentage")

genus_output <- genus_in %>% 
	separate(samp_name, 
					 into = c("sampleID","db_name"), 
					 sep = "COMP_", remove = F, extra = "merge") %>% 
	mutate(sampleID = str_remove(samp_name, paste0("_",db_name)))
genus_output$siteID = substr(genus_output$sampleID, 1, 4)

saveRDS(genus_output,"data/classification/taxonomic_rank_summaries/genus/bracken_genus_HARV.rds")



