# Tidy NEON soil environmental metadata to match with metagenome samples
# Input: NEON soil chemistry data via load_soilChem() helper function
#        Includes: soil chemistry (d15N, organicd13C, nitrogenPercent, organicCPercent),
#                  soil moisture, temperature, pH, coordinates, elevation
#        Uses helper functions: load_soilChem(), load_genSampleExample()
# Output: environmental_metadata_NEON.csv aggregated by genomicsSampleID
#         Averages multiple measurements per genomicsSampleID where applicable

library(tidyverse)
library(neonUtilities)
source("scripts/helper_functions.r")

soilChem <- load_soilChem()
genomicSamples <- load_genSampleExample(soilChem)

chemEx <- soilChem$sls_soilChemistry %>%
    select(sampleID, d15N, organicd13C, nitrogenPercent, organicCPercent)

moistureTab <- soilChem$sls_soilMoisture %>%
    select(sampleID, soilMoisture)

temperatureTab <- soilChem$sls_soilCoreCollection %>%
    select(siteID, sampleID, soilTemp)

combinedTab <- left_join(genomicSamples, chemEx, by = "sampleID") %>%
    left_join(temperatureTab, by = "sampleID") %>%
    left_join(moistureTab, by = "sampleID", relationship = "many-to-many")

soilpH_Example <- soilChem$sls_soilpH %>%
    filter(sampleID %in% combinedTab$sampleID) %>%
    select(sampleID, soilInWaterpH, soilInCaClpH)

soil_metadata_all <- combinedTab %>%
    left_join(soilpH_Example, by = "sampleID") %>%
    group_by(genomicsSampleID) %>%
    {left_join(
        summarize_at(., vars("d15N", "organicd13C", "nitrogenPercent", "organicCPercent", "soilTemp", "soilMoisture"), mean),
        summarize_at(., vars("soilInWaterpH", "soilInCaClpH"), mean_pH)
    )} %>%
    mutate_all(~ifelse(is.nan(.), NA, .))

coordinate_data <- soilChem$sls_soilCoreCollection %>%
    filter(sampleID %in% combinedTab$sampleID) %>%
    select(latitude = decimalLatitude, longitude = decimalLongitude, elevation,
           sampleTiming, nlcdClass, sampleID)

genomicSamples_coordinate_data <- left_join(genomicSamples, coordinate_data) %>%
    select(-sampleID) %>%
    distinct()

soil_metadata_all <- left_join(soil_metadata_all, genomicSamples_coordinate_data)

dir.create("data/environmental_data", recursive = TRUE, showWarnings = FALSE)
write.csv(soil_metadata_all, "data/environmental_data/environmental_metadata_NEON.csv", row.names = FALSE)
