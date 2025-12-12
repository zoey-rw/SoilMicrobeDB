
library(tidyverse)
library(neonUtilities)
library(geosphere)
library(data.table)

# Get fungi-bacteria ratios estimated by Yu et al. (spatial-environmental model) 
# By matching lat/lon to NEON sample lat/lon

# Function to get nearest date within list of points
wherenearest <- function(myPoint, allPoints){
    d <- abs(allPoints-myPoint[1])
    index <- which.min(d)
    return( index )
}


# https://stackoverflow.com/questions/21977720/r-finding-closest-neighboring-point-and-number-of-neighbors-within-a-given-rad
# This function creates a big dataframe with every possible combination
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

#fun_bac_ratio = read_csv("/projectnb/dietzelab/zrwerbin/N-cycle/data/ref_data/fun_bac_ratio.csv")

# Read in NEON soil chemistry to link sample IDs, and lat/lon
soilChem <- readRDS("/projectnb/dietzelab/zrwerbin/N-cycle/neon_soil_data_2023.rds")
gen_samp_key = soilChem$sls_soilCoreCollection %>% 
    select(siteID, sampleID, geneticSampleID) %>% 
    distinct(.keep_all=T)
samp_key = merge(samp_key, gen_samp_key, all.x=T)

sample_lat_lon = soilChem$sls_soilCoreCollection %>%
	select(sampleID, geneticSampleID, lat = decimalLatitude, lon = decimalLongitude) %>%
	distinct(.keep_all=T) 

sample_lat_lon$compositeSampleID = genSampleExample[match(sample_lat_lon$sampleID,
																										 genSampleExample$sampleID),]$genomicsSampleID
sample_lat_lon <- sample_lat_lon %>% filter(!is.na(compositeSampleID)) %>% distinct(.keep_all = T)


# Read in ratios of fungal/bacterial abundance from Yu et al. 2022
fun_bac_ratios = data.table::fread("/projectnb/dietzelab/zrwerbin/N-cycle/data/ref_data/fun_bac_ratio.csv")
fun_bac_ratios$Lat = round(fun_bac_ratios$Lat, 3)
fun_bac_ratios$Lon = round(fun_bac_ratios$Lon, 3)
fun_bac_ratios <- fun_bac_ratios %>% select(-V1) %>% group_by(Lat, Lon) %>%  distinct(.keep_all = T)


df_expanded = expand.grid.df(sample_lat_lon[,c("lat","lon")],
														 fun_bac_ratios[,c("Lat","Lon")])


# add sample names back in
df_expanded_samples = left_join(df_expanded, sample_lat_lon)

shortest_distance <- df_expanded_samples %>%
	mutate(distance = distHaversine(p1 = cbind(lon,lat),
																	p2 = cbind(Lon,Lat))) %>% 
	group_by(geneticSampleID) %>% 
	slice(which.min(distance))

fun_bac_ratios_simple = fun_bac_ratios %>%
	mutate(fungal_proportion =  round(Fun_pro, 3)) %>% 
	select(-fungal_proportion) %>% 
	distinct(Lat,Lon, .keep_all = T) #%>% ungroup 

fun_bac_proportion = left_join(shortest_distance, fun_bac_ratios_simple, 
															 by = join_by(Lat, Lon)) %>% 
	select(-c(Lat, Lon, distance))  %>% 
	ungroup %>% 
	select(-c(lat,lon))
fun_bac_proportion$siteID = substr(fun_bac_proportion$geneticSampleID, 1, 4)

write_csv(fun_bac_proportion, "data/comparison_data/modeled_fb_ratio/fun_bac_ratio_NEON.csv")

# Take a look at variation across sites
ggplot(fun_bac_proportion) + geom_point(aes(x = siteID, y = Fun_pro))
