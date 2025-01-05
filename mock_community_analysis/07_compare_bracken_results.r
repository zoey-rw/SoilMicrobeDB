library(tidyverse)
library(pavian)
library(scales)


source("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/custom_pavian.r")
source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")
source("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/helper_functions.r")


# Get a list of microbes common to PlusPF, to Mycocosm, and to NCBI
prokaryote_to_sample <- readRDS("/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/01_mock_community_input/prokaryote_to_sample.rds")
fungi_to_sample <- readRDS("/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/01_mock_community_input/fungi_to_sample.rds")
fungi_in_any <- readRDS("/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/01_mock_community_input/fungi_any.rds")

genus_bracken_dir = "/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/06_bracken_output/"

sim_reads_dir = "/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/02_simulated_reads"
abundances_list <- list.files(sim_reads_dir, pattern = "abundance.txt", recursive = T, full.names = T)

all_abundance_files <- abundances_list %>%
	setNames(., sub("_abundance.txt", "", basename(.))) %>%
	map(read.table)

mc_input_dir = "/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/01_mock_community_input/"
genome_names_list <- list.files(mc_input_dir, pattern = "genomenames", recursive = T, full.names = T)

all_genome_names_files <- genome_names_list %>%
	setNames(., sub("mc_input_genomenames_|\\.txt", "", basename(.))) %>%
	map(read.table)

genome_names = data.table::rbindlist(all_genome_names_files, idcol = "fungal_proportion", fill = T)  %>% 
	mutate(portal = V1,
				 fungal_proportion = gsub("fungprop_","",fungal_proportion, fixed=T),
				 fungal_proportion = gsub(".txt","",fungal_proportion, fixed=T))

true_abundances_in = data.table::rbindlist(all_abundance_files, idcol = "samp_name", fill = T)   %>% 
	mutate(genome = basename(V1),
				 portal = gsub("_AssemblyScaffolds_Repeatmasked.fasta","",genome, fixed=T), abundance = V2)  %>% 
	mutate(genome_names = 
				 	gsub("_genomic.fna.gz", "", portal, fixed = T))
true_abundances <- true_abundances_in


species_bracken_dir = "/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/05_architeuthis_kraken_output/"


bracken_soil_microbe_db_files <- list.files(species_bracken_dir, pattern = "soil_microbe_db_filtered_kraken_bracken_genuses.kreport", 
																						recursive = T, full.names = T)

names(bracken_soil_microbe_db_files) = basename(bracken_soil_microbe_db_files) %>% gsub("_soil_microbe_db_filtered_kraken_bracken_genuses.kreport","",.)

bracken_mock_filtered_soil_microbe_db = many_files_to_matrix_list(files = bracken_soil_microbe_db_files, 
																									filter.tax.level = "G", 
																									include.unclassified = T, percentages = T)[[1]] %>%
	as.data.frame() %>%
	mutate(db_name="soil_microbe_db") %>% 
	rownames_to_column("taxon") %>%
	pivot_longer(cols = where(is.numeric), names_to="samp_name",values_to="percentage") %>% 
	mutate(genus=tolower(taxon), predicted_abundance = percentage/100)






bracken_gtdb_207_files <- list.files(species_bracken_dir, pattern = "gtdb_207_unfiltered_filtered_kraken_bracken_genuses.kreport", 
																						recursive = T, full.names = T)

names(bracken_gtdb_207_files) = basename(bracken_gtdb_207_files) %>% 
    gsub("_gtdb_207_unfiltered_filtered_kraken_bracken_genuses.kreport","",.)

bracken_mock_filtered_gtdb_207 = many_files_to_matrix_list(files = bracken_gtdb_207_files, 
																																	filter.tax.level = "G", 
																																	include.unclassified = T, percentages = T)[[1]] %>%
	as.data.frame() %>%
	mutate(db_name="gtdb_207") %>% 
	rownames_to_column("taxon") %>%
	pivot_longer(cols = where(is.numeric), names_to="samp_name",values_to="percentage") %>% 
	mutate(genus=tolower(taxon), predicted_abundance = percentage/100)



bracken_pluspf_files <- list.files(species_bracken_dir, pattern = "pluspf_filtered_kraken_bracken_genuses.kreport", 
																		 recursive = T, full.names = T)

names(bracken_pluspf_files) = basename(bracken_pluspf_files) %>% gsub("_pluspf_filtered_kraken_bracken_genuses.kreport","",.)

bracken_mock_filtered_pluspf = many_files_to_matrix_list(files = bracken_pluspf_files, 
																													 filter.tax.level = "G", 
																													 include.unclassified = T, percentages = T)[[1]] %>%
	as.data.frame() %>%
	mutate(db_name="pluspf") %>% 
	rownames_to_column("taxon") %>%
	pivot_longer(cols = where(is.numeric), names_to="samp_name",values_to="percentage") %>% 
	mutate(genus=tolower(taxon), predicted_abundance = percentage/100)




mock_filtered = rbindlist(list(bracken_mock_filtered_gtdb_207, 
															 bracken_mock_filtered_soil_microbe_db, 
															 bracken_mock_filtered_pluspf))




true_abundances$fun_genus = fungi_to_sample[match(true_abundances$genome_names,fungi_to_sample$genome_names),]$genus
true_abundances$fun_genus = tolower(true_abundances$fun_genus)
fungal_mock_genus = true_abundances %>% filter(!is.na(fun_genus)) %>% distinct(samp_name, fun_genus)

fungal_mock_genus$in_fungal_mock_community=T

bacteria_to_sample = prokaryote_to_sample
bacteria_to_sample$portal = basename(bacteria_to_sample$fasta_file_path) %>% gsub(".fna.gz",".fasta",.)
true_abundances$bac_genus = bacteria_to_sample[match(true_abundances$portal,bacteria_to_sample$portal),]$genus
true_abundances$genus = ifelse(!is.na(true_abundances$fun_genus), true_abundances$fun_genus, true_abundances$bac_genus)


mock_genus = true_abundances %>% filter(!is.na(genus)) %>% distinct(samp_name, genus)
mock_genus$in_mock_community=T

true_abundances <- left_join(true_abundances, fungal_mock_genus)
true_abundances <- left_join(true_abundances, mock_genus)

# true_abundances$in_fungal_mock_community = ifelse(tolower(true_abundances$genus) %in% fungal_mock_genus, T, F)
# true_abundances$in_mock_community = ifelse(tolower(true_abundances$genus) %in% mock_genus, T, F)


eval_df = merge(true_abundances, mock_filtered, all=T, by=c("genus","samp_name"), allow.cartesian=TRUE) %>% 
	filter(!is.na(db_name))
eval_df = eval_df %>% 
	separate(samp_name, into = c("fungal_proportion","readdepth"), sep = "readdepth", remove = F) %>% 	
	mutate(fungal_proportion = gsub("fungprop_","",fungal_proportion, fixed=T),
				 fungal_proportion = gsub("_","",fungal_proportion, fixed=T),
				 readdepth = gsub("_","",readdepth, fixed=T))

eval_df$fungal = ifelse(tolower(eval_df$genus) %in% fungi_in_any, "Fungi", "Prokaryote")
eval_df$in_fungal_mock_community = ifelse(is.na(eval_df$in_fungal_mock_community), F, T)
eval_df$in_mock_community = ifelse(is.na(eval_df$in_mock_community), F, T)

# Set true abundance to zero if not in mock community
eval_df[eval_df$in_mock_community==F,]$abundance=0

eval_df = eval_df %>% mutate(db_name = recode(db_name,
																							
														# mycocosm_published = "Mycocosm",
														# mycocosm_all =  "Mycocosm (unpublished)" ,
														gtdb_207= "GTDB 207",
														#with_refsoil = "SoilMicrobeDB",
														soil_microbe_db = "SoilMicrobeDB",
														pluspfp8 = "PlusPFP8",
														pluspf = "PlusPF"),
														readdepth=recode(readdepth,
																						".1M" = 100000,
																						".5M" = 500000,
																						"1M" = 1000000,
																						"5M" = 5000000,
																						"10M"= 10000000,
																						"20M" = 20000000),
														fungal_proportion=as.numeric(fungal_proportion))

write.csv(eval_df, "/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/evaluation_results.csv")




# BELOW THIS LINE IS OLD CODE #


# Estimated at 5% but not in mock community?
"Bremerella" 
"Xylanibacter"

# Estimated at 0%
"Capnocytophaga"
"Paeniclostridium"

ggplot(eval_df %>% filter(samp_name == "fungprop_10_readdepth_5M" & in_mock_community) ,
			 aes(y = predicted_abundance, x=db_name, color=fungal)) +
	#	geom_histogram() +
	geom_point(#aes(color=fungal),
		size=2, alpha=.4,
		#position=position_jitter(width = .1, height=0.01)) +
		position=position_jitterdodge(jitter.height =0, jitter.width = .2, dodge.width = .9)) +
	geom_violin(draw_quantiles = c(.5)) +
	#	geom_smooth(method="lm") +
	theme_bw(base_size = 18) + 	#scale_y_log10() + 
	#facet_grid(~fungal, drop=T, scales="free")  +
	#theme_ridges() + 
	geom_hline(yintercept = .005, linetype=2)  + 
	annotate("text", x = .6, y = .0047, label = "True abundance") + 
	ylab("Predicted abundance") + 
	xlab("Database")  + guides(color=guide_legend(NULL)) + ggtitle("Mock community of 200 species")






ggplot(eval_df %>% filter(in_fungal_mock_community),
			 aes(x = predicted_abundance, y = abundance)) +
	geom_point(aes(color=as.factor(in_fungal_mock_community)), 
						 size=3, alpha=.3, position=position_jitter(width = .0001, height=.0001)) +
	#	geom_smooth(method="lm") +
	theme_bw(base_size = 18) + 	#scale_y_log10() + 
	facet_grid(fungal_proportion+readdepth~db_name, drop=T, scales="free") + 
	ylim(c(0,.01)) +	
	stat_poly_line() +
	stat_poly_eq() + ggtitle("Fungal mock community abundance accuracy (Bracken)") + 
	geom_label_repel(aes(label = genus), max.overlaps = 15) + theme(legend.position="bottom") 


ggplot(eval_df %>% filter(!in_fungal_mock_community),
			 aes(x = predicted_abundance, y = abundance, color=db_name)) +
	geom_boxplot(position=position_dodge(width=1)) +
	geom_point(size=3, alpha=.1, position=position_jitter(width = .0001, height=.0001)) +
	#	geom_smooth(method="lm") +
	theme_bw(base_size = 18) + 	#scale_y_log10() + 
	facet_grid(fungal_proportion~readdepth, drop=T, scales="free") + 
	ylim(c(0,.01)) #+	
#	stat_poly_line() +
#	stat_poly_eq() 



library(ggridges)
ggplot(eval_df %>% filter(readdepth != "10m") %>% filter(in_mock_community),
			 aes(x = predicted_abundance, y=db_name)) +
#	geom_histogram() +
	geom_density_ridges(scale = 1) +
	#geom_point(size=3, alpha=.1, position=position_jitter(width = .0001, height=.0001)) +
	#	geom_smooth(method="lm") +
	theme_bw(base_size = 18) + 	#scale_y_log10() + 
	facet_grid(fungal_proportion~readdepth, drop=T, scales="free")  +
	theme_ridges() + geom_vline(xintercept = .005, linetype=2) + ggtitle("Database abundance estimates")


ggplot(eval_df %>% filter(readdepth != "10m") %>% filter(in_mock_community),
			 aes(y = predicted_abundance, x=db_name)) +
	#	geom_histogram() +
	geom_violin() +
	# geom_point(#aes(color=in_mock_community), 
	# 					 size=2, alpha=.01, 
	# 					 position=position_jitter(width = .1, height=0.01)) +
	#position=position_jitterdodge(jitter.height =0.001, jitter.width = .1, dodge.width = 1)) +
	
	#	geom_smooth(method="lm") +
	theme_bw(base_size = 18) + 	#scale_y_log10() + 
	facet_grid(fungal_proportion~readdepth, drop=T, scales="free")  +
	theme_ridges() + geom_hline(yintercept = .005, linetype=2) 

# 
# eval_df$fungal = ifelse(eval_df$in_fungal_mock_community==T, T, F)
# eval_df$in_mock_community = ifelse(eval_df$in_mock_community==T, T, F)




ggplot(eval_df %>% filter(db_name == "SoilMicrobeDB" & in_fungal_mock_community),
			 aes(y = predicted_abundance, x=as.numeric(readdepth))) +
	#	geom_histogram() +
	#geom_violin() +
	 geom_point(
	 					 size=2, alpha=.1, 
	 					 position=position_jitter(width = .1, height=0)) +
	#position=position_jitterdodge(jitter.height =0.001, jitter.width = .1, dodge.width = 1)) +
	
	#	geom_smooth(method="lm") +
	theme_bw(base_size = 18) + 	#
	#scale_y_log10() + 
	facet_grid(~fungal_proportion, drop=T, scales="free")  +
	scale_y_sqrt() +
	geom_hline(yintercept = .005, linetype=2) 





fig_3a = ggplot(metrics_to_plot %>% filter(metric == "Bias" &
                              fungal_proportion==10),
       aes(y = value, x=db_name, color=fungal)) +
    geom_point(
        size=2, alpha=.4, 
       # position=position_jitter(width = .1, height=0)) +
    position=position_jitterdodge(jitter.height =0, jitter.width = .1, dodge.width = .8)) +
    
    theme_classic(base_size = 18) + 
    facet_wrap(~fungal, scales="free_x") + 
    ylab("Bias (simulated community)") + 
    xlab(NULL) + 
    geom_hline(yintercept = 0, linetype=2)  +
    scale_color_discrete(name="") 
fig_3b = ggplot(metrics_to_plot %>% filter(metric != "Bias" &
                                               fungal_proportion==10),
                aes(y = value, x=db_name, color=fungal)) +
    geom_point(
        size=2, alpha=.4, 
        position=position_jitter(width = .1, height=0)) +
#    position=position_jitterdodge(jitter.height =0, jitter.width = .1, dodge.width = .8)) +
    
    theme_classic(base_size = 18) + 
    facet_wrap(~fungal, scales="free_x") + 
    ylab("RMSE (simulated community)") + 
    xlab("Database") + 
    geom_hline(yintercept = 0, linetype=2)  +
    scale_color_discrete(name="") + 
    stat_compare_means(data = metrics_to_plot %>% filter(metric != "Bias" &
                                                             fungal_proportion==10 & fungal=="Fungi"), inherit.aes = T,
                       comparisons = list(c("PlusPF","SoilMicrobeDB")), method = "t.test", label = "p.signif")  + 
    stat_compare_means(data = metrics_to_plot %>% filter(metric != "Bias" &
                                                             fungal_proportion==10 & fungal !="Fungi"), inherit.aes = T,
                       comparisons = list(  c("PlusPF","SoilMicrobeDB"),
                                            c("GTDB 207","SoilMicrobeDB"),
                                            c("PlusPF","GTDB 207")), method = "t.test", label = "p.signif")
fig_3a
fig_3b


metrics_to_plot %>% filter(metric == "Bias" &
                               fungal_proportion==10),
aes(y = value, x=db_name, color=fungal)) +

ggarrange(fig_3a, fig_3b, nrow = 2, labels = c("A","B"), common.legend = T)


rmse_df <- eval_df %>% 
	group_by(db_name, fungal_proportion, readdepth) %>% 
	rmse(truth = abundance, estimate = predicted_abundance)

adist_df <- eval_df %>% 
	group_by(db_name, fungal_proportion, readdepth, fungal) %>% 
	filter(abundance > 0) %>% 
	summarize(adist = aDist(x = abundance, y = predicted_abundance))


dist_mat_df <- eval_df %>% 
	#group_by(db_name, fungal_proportion, readdepth, fungal) %>% 
	select(abundance, predicted_abundance) #%>% 

vegdist(dist_mat_df, method="robust.aitchison")

fVeg = function(data){
	xveg = data %>% select(abundance, predicted_abundance)
	
	AD = vegdist(xveg, "robust.aitchison")
	AD
#	data %>% mutate(AD = rep(vegdist(xveg, "robust.aitchison"), nrow(xveg)))
} 

xveg = eval_df[1:100,] %>% select(abundance, predicted_abundance)
vegdist(xveg, "robust.aitchison")
fVeg(xveg)

eval_df %>% 
	group_by(db_name, fungal_proportion, readdepth, fungal) %>% 
	group_modify(~fVeg(.x))


aDist(eval_df$abundance, eval_df$predicted_abundance)
mape

library(robCompositions)

aDist
library(vegan)




ggplot(rmse_df %>% filter(readdepth < 1e7),
       aes(y = .estimate, x=db_name, color=db_name)) +
    geom_point(
        size=2, alpha=.5, 
        position=position_jitter(width = .1, height=0)) +
    #position=position_jitterdodge(jitter.height =0.001, jitter.width = .1, dodge.width = 1)) +
    
    #	geom_smooth(method="lm") +
    theme_bw(base_size = 18)  	#
    #facet_grid(db_name, drop=T, scales="free")
