# Tidy various NEON soil data measurements to match with metagenome samples
# Workflow from https://www.neonscience.org/resources/learning-hub/tutorials/neon-data-metagenomics

library(respirometry)
library(neonUtilities)
library(tidyverse)
#devtools:install_github('NEONScience/NEON-geolocation/geoNEON', dependencies=TRUE)
library(geoNEON)

# downloaded with loadByProduct() for all sites, product 'DP1.10086.001'
soilChem = readRDS("/projectnb/dietzelab/zrwerbin/N-cycle/neon_soil_data_2023.rds")

# Get genomic samples
genomicSamples <- soilChem$sls_metagenomicsPooling %>%
	tidyr::separate(genomicsPooledIDList, into=c("first","second","third"),sep="\\|",fill="right") %>%
	dplyr::select(genomicsSampleID,first,second,third) %>% 
	tidyr::pivot_longer(cols=c("first","second","third"),values_to = "sampleID") %>%
	dplyr::select(sampleID,genomicsSampleID) %>%
	drop_na()

# Get soil chemistry (not pH)
chemEx <- soilChem$sls_soilChemistry %>%
	dplyr::select(sampleID,d15N,organicd13C,nitrogenPercent,organicCPercent)

# Get soil moisture
moistureTab <- soilChem$sls_soilMoisture %>%
	dplyr::select(sampleID,soilMoisture)

# Get soil temp
temperatureTab <- soilChem$sls_soilCoreCollection %>%
	dplyr::select(siteID,sampleID,soilTemp)

## now combine the tables 
combinedTab <- left_join(genomicSamples,chemEx, by = "sampleID")# %>% drop_na()
combinedTab <- left_join(combinedTab,temperatureTab, by = "sampleID")# %>% drop_na()
combinedTab <- left_join(combinedTab,moistureTab, by = "sampleID",relationship = "many-to-many")# %>% drop_na()

# Get pH
soilpH_Example <- soilChem$sls_soilpH %>%
	dplyr::filter(sampleID %in% combinedTab$sampleID) %>%
	dplyr::select(sampleID,soilInWaterpH,soilInCaClpH)

combinedTab_pH <- left_join(combinedTab,soilpH_Example, by = "sampleID")

soil_metadata_all <- combinedTab_pH %>%
	group_by(genomicsSampleID) %>%
	{left_join(
		summarize_at(.,vars("d15N","organicd13C","nitrogenPercent",
												"organicCPercent","soilTemp","soilMoisture"), mean),
		summarize_at(.,vars("soilInWaterpH","soilInCaClpH"), mean_pH)
	)} %>% 
		# replace NaN with NA
	 	mutate_all(~ifelse(is.nan(.), NA, .))

# Takes a few minutes to add coordinate data
# coordinate_data <- getLocTOS(soilChem$sls_soilCoreCollection %>%
# 														 	dplyr::filter(sampleID %in% combinedTab$sampleID) , 'sls_soilCoreCollection')
# coordinate_data_only <- coordinate_data %>% select(namedLocation, latitude = adjDecimalLatitude, 
# 																									 longitude = adjDecimalLongitude, 
# 																									 elevation = adjElevation)


coordinate_data <- soilChem$sls_soilCoreCollection %>%
	dplyr::filter(sampleID %in% combinedTab$sampleID) %>%
	dplyr::select(latitude = decimalLatitude, 
								longitude = decimalLongitude,
								elevation,
								sampleTiming,
								nlcdClass,
								# collectDate,
								# sampleTopDepth,
								# sampleBottomDepth,
								sampleID)
genomicSamples_coordinate_data <- left_join(genomicSamples,coordinate_data) %>% select(-sampleID) %>% distinct()
soil_metadata_all <- left_join(soil_metadata_all, genomicSamples_coordinate_data)
write.csv(soil_metadata_all, "/projectnb/frpmars/soil_microbe_db/ref_data/environmental_metadata_NEON.csv")

