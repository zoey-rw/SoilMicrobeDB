library(neonUtilities)
library(tidyverse)



# Somewhat exploratory - trying to determine the best way to filter and aggregate qPCR data
# without being able to match lab QC results

qpcrData <- loadByProduct(site = "all", dpID = "DP1.10109.001", package = "expanded", check.size = F, #token = neon_token, 
                          include.provisional = T)
## qpcr_in = stackByTable(filepath = "/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/misc_scripts/neon_qPCR/NEON_group-abund-microbe-soil")
#qpcr_in <- qpcrData$mga_soilGroupAbundances

saveRDS(qpcrData, "data/comparison_data/qpcr/NEON_qpcr_full.rds")

qpcrData <- readRDS("data/comparison_data/qpcr/NEON_qpcr_full.rds")




# Aggregate mean/sd, function from: https://stackoverflow.com/questions/9222056/existing-function-to-combine-standard-deviations-in-r
## N: vector of sizes
## M: vector of means
## S: vector of standard deviations
grand.mean <- function(M, N) {weighted.mean(M, N)}
grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                          weighted.mean(M, N)^2)}

qc_info1 = qpcrData$mga_batchResults %>% 
    filter(standardDescription != "genomic DNA") %>%  # couldnt figure out what this meant
    mutate(targetTaxonGroup = ifelse(grepl("Fung",standardDescription),"fun","bac_arc")) %>% 
    select(batchID, targetTaxonGroup, pcrEfficiency) %>% distinct(.keep_all=T)
qc_info2 = qpcrData$mga_soilGroupAbundances %>% 
    filter(is.na(dataQF) | dataQF %in% c("legacyData")) %>% 
    select(dnaSampleID, batchID, nucleicAcidConcentration, meanCqValue, targetTaxonGroup) %>% 
    mutate(targetTaxonGroup = recode(targetTaxonGroup, 
                                     "bacteria and archaea" = "bac_arc",
                                     "bacteria"="bac",
                                     "archaea"="arc","fungi"="fun")) %>% distinct(.keep_all=T) #%>%
qc_info3 = left_join(qc_info2, qc_info1, relationship = "many-to-many") %>% #select(-batchID) %>% 
    distinct(.keep_all=T) 
qc_info =  qc_info3 %>%  pivot_wider(names_from="targetTaxonGroup",values_from = c("nucleicAcidConcentration", "meanCqValue","pcrEfficiency"), values_fun="mean")



qpcrData_qc = qpcrData$mga_soilGroupAbundances %>%
	select(siteID,plotID, collectDate, geneticSampleID, dnaSampleID, targetTaxonGroup, dataQF,
				 #nucleicAcidConcentration, 
				 meanCopyNumber,copyNumberStandardDeviation) %>%
	#filter(dataQF == "LLOQ")
filter(is.na(dataQF) | dataQF %in% c("legacyData")) %>% 
	mutate("copyNumberVariance" = copyNumberStandardDeviation^2)

bac_arc = qpcrData_qc %>% filter(targetTaxonGroup %in% c("bacteria","archaea")) %>% 
	select(dnaSampleID, geneticSampleID, targetTaxonGroup, meanCopyNumber, copyNumberStandardDeviation, dataQF) %>% 
	pivot_wider(names_from="targetTaxonGroup",values_from = c("meanCopyNumber", "copyNumberStandardDeviation"))#

non_bac_arc = qpcrData_qc %>% filter(targetTaxonGroup %in% c("bacteria and archaea")) %>% 
	select(dnaSampleID, geneticSampleID, targetTaxonGroup, meanCopyNumber, copyNumberStandardDeviation, dataQF) %>% 
	pivot_wider(names_from="targetTaxonGroup",values_from = c("meanCopyNumber", "copyNumberStandardDeviation"))#

fung_only = qpcrData_qc %>% filter(targetTaxonGroup %in% c("fungi")) %>% 
	select(dnaSampleID, geneticSampleID, targetTaxonGroup, meanCopyNumber, copyNumberStandardDeviation, dataQF) %>% 
	pivot_wider(names_from="targetTaxonGroup",values_from = c("meanCopyNumber", "copyNumberStandardDeviation"))#

# 
# bac_arc %>% mutate(combined_sd = grand.sd(S = c(copyNumberStandardDeviation_bacteria, copyNumberStandardDeviation_archaea), 
# 																				M =c(meanCopyNumber_bacteria, meanCopyNumber_archaea), 
# 																				N = c(3, 3)))

combined_sd = list()
for (i in 1:nrow(bac_arc)){
	x = bac_arc[i,]
	combined_sd[i] = grand.sd(#S = c(x["copyNumberStandardDeviation_bacteria"],x["copyNumberStandardDeviation_archaea"]), 
		
		#			M = c(x["copyNumberStandardDeviation_bacteria"],x["copyNumberStandardDeviation_archaea"]), 
		S =c(x$copyNumberStandardDeviation_bacteria[i], x$copyNumberStandardDeviation_archaea[i]), 
		M =c(x$meanCopyNumber_bacteria[i], x$meanCopyNumber_archaea[i]), 
		N = c(3, 3))
}
bac_arc$`copyNumberStandardDeviation_bacteria and archaea` = unlist(combined_sd)
bac_arc$`meanCopyNumber_bacteria and archaea` = bac_arc$meanCopyNumber_bacteria + bac_arc$meanCopyNumber_archaea
bac_arc$copyNumberStandardDeviation_bacteria=NULL
bac_arc$copyNumberStandardDeviation_archaea=NULL
bac_arc$meanCopyNumber_bacteria=NULL
bac_arc$meanCopyNumber_archaea=NULL

non_fungi_qpcrData = data.table::rbindlist(list(non_bac_arc,  
																			bac_arc), use.names = T)

qpcrData_out = merge(fung_only, non_fungi_qpcrData, all=T, by = c("geneticSampleID","dnaSampleID")) %>% 
	mutate(legacy = ifelse(grepl("legacy",dataQF.x)|grepl("legacy",dataQF.y), T, F)) %>% 
	mutate(fb_ratio = `meanCopyNumber_fungi`/`meanCopyNumber_bacteria and archaea`,
				 f_proportion = `meanCopyNumber_fungi`/(`meanCopyNumber_fungi` + `meanCopyNumber_bacteria and archaea`),
				 siteID = substr(dnaSampleID, 1, 4)) 
qpcrData_out[qpcrData_out$legacy,]$siteID %>% table
qpcrData_out[!qpcrData_out$legacy,]$siteID %>% table


qpcrData_out <- left_join(qpcrData_out, qc_info)

ggplot(qpcrData_out %>% filter(siteID %in% c("CPER", "DSNY", "HARV", "OSBS", "STER", "TALL", "WOOD"))) + 
	geom_histogram(aes(x = `meanCopyNumber_bacteria and archaea`)) + 
	facet_grid(legacy~siteID, scales = "free") + scale_x_log10()
ggplot(qpcrData_out %>% filter(siteID %in% c("CPER", "DSNY", "HARV", "OSBS", "STER", "TALL", "WOOD"))) + 
	geom_histogram(aes(x = `meanCopyNumber_fungi`)) + 
	facet_grid(legacy~siteID, scales = "free") + scale_x_log10()
ggplot(qpcrData_out %>% filter(siteID %in% c("CPER", "DSNY", "HARV", "OSBS", "STER", "TALL", "WOOD"))) + 
	geom_histogram(aes(x = f_proportion)) + 
	facet_grid(legacy~siteID, scales = "free") + scale_x_log10()
ggplot(qpcrData_out) + geom_histogram(aes(x = fb_ratio)) + 
	facet_grid(rows=vars(legacy)) + scale_x_log10()

qpcrData_out <- qpcrData_out %>% 
	mutate(test_fun = ifelse(legacy==T, meanCopyNumber_fungi/1000, meanCopyNumber_fungi))
ggplot(qpcrData_out %>% filter(siteID %in% c("CPER", "DSNY", "HARV", "OSBS", "STER", "TALL", "WOOD"))) + 
	geom_histogram(aes(x = test_fun)) + 
	facet_grid(legacy~siteID, scales = "free") + scale_x_log10()


qpcrData_out <- qpcrData_out %>% 
	mutate(test_bac = ifelse(legacy==T, `meanCopyNumber_bacteria and archaea`/1000, `meanCopyNumber_bacteria and archaea`))

ggplot(qpcrData_out %>% filter(siteID %in% c("CPER", "DSNY", "HARV", "OSBS", "STER", "TALL", "WOOD"))) + 
	geom_histogram(aes(x = test_bac)) + 
	facet_grid(legacy~siteID, scales = "free") + scale_x_log10()


ggplot(qpcrData_out %>% filter(siteID %in% c("CPER", "DSNY", "HARV", "OSBS", "STER", "TALL", "WOOD"))) + 
    geom_point(aes(x = test_bac, y = meanCqValue_bac_arc)) + 
    facet_grid(~siteID, scales = "free") + scale_x_log10()





# Load soilCores using helper function
if(!exists("soilCores")) {
    source("scripts/helper_functions.r")
    soilCores <- load_soilCores()
}


qpcrData_out$sampleID = soilCores[match(qpcrData_out$geneticSampleID,
																		 soilCores$geneticSampleID),]$sampleID
qpcrData_out$compositeSampleID = genSampleExample[match(qpcrData_out$sampleID,
																										 genSampleExample$sampleID),]$genomicsSampleID

write_csv(qpcrData_out, "data/comparison_data/qpcr/NEON_qpcr.csv")
