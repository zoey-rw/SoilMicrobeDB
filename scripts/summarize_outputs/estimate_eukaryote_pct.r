# Create list of high-abundance taxa, from mapping NEON metagenomes to SoilGenome DB

library(tidyverse)
library(pavian)
library(ggpubr)
source("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")


sampleID_pooled_key <- readRDS("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/NEON_metagenomes/NEON_pooled_sample_key.rds")
pooled_recode_list <- sampleID_pooled_key$genomicsSampleID
names(pooled_recode_list) <- sampleID_pooled_key$sampleID

soil_sample_dir <- "/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output"

samp_files <- list.files(soil_sample_dir, recursive=T, pattern = "_filtered_kraken_bracken_domains.kreport", full.names = T)
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


saveRDS(domain_output,"/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/bracken_domain_estimates.rds")


domain_output %>% filter(db_name %in% c("soil_microbe_db","pluspf")) %>% 
	#	filter(seq_depth > 5000000) %>%
#	filter(taxon %in% c("Archaea","Bacteria","Eukaryota")) %>% 
	filter(taxon %in% c("Eukaryota")) %>% 
	group_by(db_name,taxon) %>% 
	dplyr::summarize( mean(percentage))


options(scipen=999)
ggplot(domain_output %>% filter(db_name %in% c("soil_microbe_db","pluspf")) %>% 
			 #	filter(seq_depth > 5000000) %>%
			 	filter(taxon %in% c("Archaea","Bacteria","Eukaryota")),
			 aes(x = db_name, y = percentage, color=db_name)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_jitter( size=3, alpha=.8,
						 #position=position_dodge(width = 1)
	) +
	#geom_smooth(aes(x = seq_depth, y = percentage, color=db_name)) +
	facet_grid(rows=vars(taxon), scale="free_y")  +
	theme_linedraw(base_size = 20) +
	ylab("Estimated abundances") +
	scale_y_sqrt() +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) + 
	stat_compare_means(aes(x = db_name, y = percentage, color=db_name))


ggplot(domain_output %>% filter(db_name %in% c("soil_microbe_db","pluspf")) %>% 
			 	#	filter(seq_depth > 5000000) %>%
			 	filter(taxon %in% c("Eukaryota")),
			 aes(x = db_name, y = percentage, color=db_name)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_jitter( size=3, alpha=.8,
							 #position=position_dodge(width = 1)
	) +
	#geom_smooth(aes(x = seq_depth, y = percentage, color=db_name)) +
	facet_grid(rows=vars(siteID), scale="free_y")  +
	theme_linedraw(base_size = 20) +
	ylab("Estimated abundances") +
	scale_y_sqrt() +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) + 
	stat_compare_means()

ggplot(domain_output %>% filter(db_name %in% c("soil_microbe_db")) %>% 
			 	#	filter(seq_depth > 5000000) %>%
			 	filter(taxon %in% c("Eukaryota")),
			 aes(x = siteID, y = percentage, color=db_name)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_jitter( size=3, alpha=.8,
							 #position=position_dodge(width = 1)
	) +
	#geom_smooth(aes(x = seq_depth, y = percentage, color=db_name)) +
	theme_linedraw(base_size = 20) +
	ylab("Estimated abundances") +
	scale_y_sqrt() +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) + 
	stat_compare_means()


