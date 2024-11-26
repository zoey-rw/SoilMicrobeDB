library(tidyverse)

filter_reads_dir = "/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/04_architeuthis_output"
filter_summary_list <- list.files(filter_reads_dir, pattern = "_summary.output", recursive = T, full.names = T)
filter_summary_files <- filter_summary_list %>%
	setNames(., sub("_summary.output", "", basename(.))) %>%
	map(read.table)


filter_reads_dir = "/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/04_architeuthis_output"
filter_scores_list <- list.files(filter_reads_dir, pattern = "_scores.output", recursive = T, full.names = T)
filter_scores_files <- filter_scores_list %>%
	setNames(., sub("_scores.output", "", basename(.))) %>%
	map(data.table::fread)

filter_scores1 = data.table::rbindlist(filter_scores_files, idcol = "samp_name", fill = T) 


filter_scores1 = filter_scores1 %>% mutate(pass_filter = ifelse(consistency > .9 &
																																	entropy < .1 &
																																	multiplicity < 2, 
																																1, 0))
pass_filter = filter_scores1 %>% 
	group_by(samp_name) %>% 
	add_tally(name = "n_reads") %>% 
	group_by(samp_name, n_reads, pass_filter) %>% 
	tally(name = "n_pass_filter") %>% 
	filter(pass_filter==1) %>% 
	mutate(percent_classified_passing = n_pass_filter / n_reads)  %>% 
	separate(samp_name, into = c(NA,"fungal_proportion",NA,"readdepth","db_name"), 
					 sep = "_", remove = F, extra = "merge") %>% 	
	mutate(fungal_proportion = gsub("fungprop_","",fungal_proportion, fixed=T),
				 fungal_proportion = as.numeric(gsub("_","",fungal_proportion, fixed=T)),
				 readdepth = gsub("_","",readdepth, fixed=T),
	readdepth=recode(readdepth,
									 ".1m" = 100000,
									 "1m" = 1000000,
									 "5m" = 5000000,
									 "10m"= 10000000,
									 "20m" = 20000000),
	percent_classified = n_reads / (readdepth/2),
	percent_passing = n_pass_filter / (readdepth/2))
	

pass_filter_long = pass_filter %>% pivot_longer(cols = c(percent_passing, percent_classified), names_to = "metric")

ggplot(pass_filter_long %>% filter(readdepth != 1e7 & metric == "percent_passing"),
			 aes(y = value, x=fungal_proportion, color = db_name)) +
	geom_point(#aes(color=fungal),
		size=2, alpha=.4,
		#position=position_jitter(width = .1, height=0.01)) +
		position=position_jitterdodge(jitter.height = .01, jitter.width = .1, dodge.width = 1)) +
	geom_line() +
	theme_bw(base_size = 18) + 
	#facet_grid(~readdepth, drop=T, scales="free")  +
	ylab("% reads classified and passing filter") + 
	xlab("% of fungi in mock community")  + 
	guides(color=guide_legend(NULL)) 

filter_scores1 %>% 
	group_by(samp_name) %>% summarize(mean_consistency = mean(consistency),
																		mean_confidence = mean(confidence),
																		mean_entropy = mean(entropy),
																		mean_multiplicity = mean(multiplicity))
one_samp = filter_scores1 %>% filter(samp_name =="fungprop_5_readdepth_1m_pluspf")


filter_scores1 %>% 
	group_by(samp_name) %>% summarize(min(consistency),
																		min(confidence))

filter_summary1 = data.table::rbindlist(filter_summary_files, idcol = "samp_name", fill = T) 

filter_summary2 = filter_summary1 %>% 
group_by(samp_name) %>%
	#add_tally(name = "readcount") %>% 
	mutate(readcount =sum (total_reads)) %>% 
	
	group_by(samp_name,readcount, in_lineage) %>%
	#tally(name = "freq") %>% 
	summarize(freq =sum (total_reads)) %>% 

	filter(in_lineage==0) %>% 
	mutate(percent_removed = freq / readcount) 
	
filter_summary3 = filter_summary2 %>% 
	mutate(samp_name = gsub("_summary.output","",samp_name, fixed=T)) %>% 
	separate(samp_name, into = c(NA,"fungal_proportion",NA,"readdepth","db_name"), 
					 sep = "_", remove = F, extra = "merge") %>% 	
	mutate(fungal_proportion = gsub("fungprop_","",fungal_proportion, fixed=T),
				 fungal_proportion = gsub("_","",fungal_proportion, fixed=T),
				 readdepth = gsub("_","",readdepth, fixed=T)) 

write_csv(filter_summary, "/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/overall_filter_results.csv")



ggplot(filter_summary,
			 aes(y = percent_removed, x=db_name, color=fungal_proportion)) +
	#	geom_histogram() +
	geom_point(#aes(color=fungal),
		size=2, alpha=.4,
		#position=position_jitter(width = .1, height=0.01)) +
		position=position_jitterdodge(jitter.height =0, jitter.width = .2, dodge.width = .9)) +
#	geom_violin(draw_quantiles = c(.5)) +
	#	geom_smooth(method="lm") +
	theme_bw(base_size = 18) + 	#scale_y_log10() + 
	#facet_grid(~fungal, drop=T, scales="free")  +
	#theme_ridges() + 
	ylab("") + 
	xlab("Database")  + guides(color=guide_legend(NULL)) #+ ggtitle("Mock community of 200 species")





filter_summary3

filter_summary = read_csv("/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/04_architeuthis_output/fungprop_20_readdepth_1m_pluspfp8_summary.output")


table(filter_summary$in_lineage)
