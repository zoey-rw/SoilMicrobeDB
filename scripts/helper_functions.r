
# Define functions to compute correlations (this could be rewritten as one function with a second argument for variable)
cor_fun_pH <- function(df) cor.test(df$percentage, 
																		df$soilInCaClpH, 
																		method = "spearman", exact = FALSE) %>% 
	tidy()
cor_fun_nitr <- function(df) cor.test(df$percentage, 
																			df$nitrogenPercent, 
																			method = "spearman", exact = FALSE) %>% 
	tidy()
cor_fun_carbon <- function(df) cor.test(df$percentage, 
																				df$organicd13C, 
																				method = "spearman", exact = FALSE) %>% 
	tidy()
cor_fun_temperature <- function(df) cor.test(df$percentage, 
																						 df$soilTemp, 
																						 method = "spearman", exact = FALSE) %>% 
	tidy()
cor_fun_moisture <- function(df) cor.test(df$percentage, 
																					df$soilMoisture, 
																					method = "spearman", exact = FALSE) %>% 
	tidy()


assign_biome_presence = function(test_df) {
	
	# get counts of samples per biome 
	sample_counts = test_df %>% group_by(nlcdClass) %>% tally(name = "biome_count")
	
	# get counts of non-zero abundances per biome
	presence_counts = test_df %>% mutate(present = ifelse(percentage > 0, 1, 0)) %>% 
		group_by(nlcdClass, present) %>% tally(name = "presence_count")
	
	
	presence_counts_expanded = presence_counts %>% ungroup %>%  expand(nlcdClass, present)
	
	# set NA values to zero
	presence_counts_expanded = merge(presence_counts, presence_counts_expanded, all=T) %>% 
		replace_na(list(presence_count = 0))
	
	df_merged = merge(presence_counts_expanded, sample_counts, all=T) %>% 
		filter(present == 1) %>% 
		mutate(prevalence = presence_count/biome_count,
					 present_in_biome = ifelse(prevalence > .02, 1, 0)) %>% 
		select(nlcdClass,biome_count, prevalence, present_in_biome)
	
	df_out = df_merged  %>% 
		select(nlcdClass,biome_count, prevalence, present_in_biome)
	return(df_out)
}

# From testing 


# Combine abundance and metadata into one dataframe 
# merged_df = merge(species_abundances, 
# 									soil_metadata, by.x = "sampleID", by.y = "genomicsSampleID")
# 
# merged_df_nest = merged_df %>% 
# 	select(-c(sampleID)) %>% 
# 	group_by(db_name, taxon) %>% nest
# 
# test_df = merged_df_nest[,3][[1]][[2]]
# 
# 
# # Filter by correlation strength
# data_nest_biome <- merged_df_nest[1:100,] %>% 
# 	mutate(biome_presence = map(data, assign_biome_presence)) 
# 
# df_biome = data_nest_biome %>% 
# 	select(-data) %>% unnest(cols = c(biome_presence)) %>% ungroup 
# 	
# df_biome_wide = df_biome %>% 
# 	select(-c(biome_count, prevalence)) %>% 
# 	pivot_wider(names_from = nlcdClass, values_from = present_in_biome)




