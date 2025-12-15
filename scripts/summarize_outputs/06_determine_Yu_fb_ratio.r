# Match NEON samples to Yu et al. spatial-environmental model fungal-bacterial ratios
# Input: NEON sample coordinates via load_soilChem() helper function
#        Yu et al. fungal-bacterial ratio predictions from /projectnb/dietzelab/zrwerbin/N-cycle/data/ref_data/fun_bac_ratio.csv
#        Uses helper functions: load_soilChem(), load_genSampleExample()
# Output: fun_bac_ratio_NEON.csv with nearest-neighbor matched ratios
#         Uses Haversine distance to find closest Yu et al. prediction point for each NEON sample

library(tidyverse)
library(data.table)
library(geosphere)

expand.grid.df <- function(...) Reduce(function(...) merge(..., by = NULL), list(...))

source("scripts/helper_functions.r")
soilChem <- load_soilChem()
genSampleExample <- load_genSampleExample(soilChem)

sample_lat_lon <- soilChem$sls_soilCoreCollection %>%
    select(sampleID, geneticSampleID, lat = decimalLatitude, lon = decimalLongitude) %>%
    distinct(.keep_all = TRUE)

sample_lat_lon$compositeSampleID <- genSampleExample[match(sample_lat_lon$sampleID,
                                                           genSampleExample$sampleID),]$genomicsSampleID
sample_lat_lon <- sample_lat_lon %>% 
    filter(!is.na(compositeSampleID)) %>% 
    distinct(.keep_all = TRUE)

fun_bac_ratio_paths <- c(
    "data/comparison_data/modeled_fb_ratio/fun_bac_ratio.csv",
    "/projectnb/dietzelab/zrwerbin/N-cycle/data/ref_data/fun_bac_ratio.csv"
)

fun_bac_ratio_file <- fun_bac_ratio_paths[file.exists(fun_bac_ratio_paths)][1]
if(is.na(fun_bac_ratio_file)) {
    stop("fun_bac_ratio.csv not found in:\n  ", paste(fun_bac_ratio_paths, collapse = "\n  "))
}

fun_bac_ratios <- fread(fun_bac_ratio_file) %>%
    mutate(Lat = round(Lat, 3), Lon = round(Lon, 3)) %>%
    select(-V1) %>%
    group_by(Lat, Lon) %>%
    distinct(.keep_all = TRUE)

df_expanded <- expand.grid.df(sample_lat_lon[, c("lat", "lon")],
                              fun_bac_ratios[, c("Lat", "Lon")])

shortest_distance <- left_join(df_expanded, sample_lat_lon) %>%
    mutate(distance = distHaversine(p1 = cbind(lon, lat), p2 = cbind(Lon, Lat))) %>% 
    group_by(geneticSampleID) %>% 
    slice(which.min(distance))

fun_bac_ratios_simple <- fun_bac_ratios %>%
    mutate(fungal_proportion = round(Fun_pro, 3)) %>%
    select(-fungal_proportion) %>%
    distinct(Lat, Lon, .keep_all = TRUE)

fun_bac_proportion <- left_join(shortest_distance, fun_bac_ratios_simple, by = join_by(Lat, Lon)) %>% 
    select(-c(Lat, Lon, distance, lat, lon)) %>%
    ungroup() %>%
    mutate(siteID = substr(geneticSampleID, 1, 4))

dir.create("data/comparison_data/modeled_fb_ratio", recursive = TRUE, showWarnings = FALSE)
write_csv(fun_bac_proportion, "data/comparison_data/modeled_fb_ratio/fun_bac_ratio_NEON.csv")
