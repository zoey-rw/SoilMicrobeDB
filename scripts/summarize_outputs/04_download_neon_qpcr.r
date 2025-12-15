# Download and process NEON qPCR data (DP1.10109.001)
# Downloads only if NEON_qpcr_full.rds doesn't exist
# Input: NEON API via neonUtilities::loadByProduct (all sites, expanded package, includes provisional)
#        Uses helper functions: load_soilCores(), load_genSampleExample()
# Output: NEON_qpcr.csv with fungal-bacterial ratios and qPCR copy numbers
#         Combines bacteria and archaea measurements, calculates combined standard deviations
#         Links to sampleID and compositeSampleID via soilCores and genSampleExample

library(neonUtilities)
library(tidyverse)
library(data.table)

qpcr_full_file <- "data/comparison_data/qpcr/NEON_qpcr_full.rds"
qpcr_output_file <- "data/comparison_data/qpcr/NEON_qpcr.csv"

if(!file.exists(qpcr_full_file)) {
    dir.create(dirname(qpcr_full_file), recursive = TRUE, showWarnings = FALSE)
    qpcrData <- loadByProduct(site = "all", dpID = "DP1.10109.001", package = "expanded", 
                             check.size = FALSE, include.provisional = TRUE)
    saveRDS(qpcrData, qpcr_full_file)
}

qpcrData <- readRDS(qpcr_full_file)

grand.mean <- function(M, N) weighted.mean(M, N)
grand.sd <- function(S, M, N) sqrt(weighted.mean(S^2 + M^2, N) - weighted.mean(M, N)^2)

qc_info1 <- qpcrData$mga_batchResults %>% 
    filter(standardDescription != "genomic DNA") %>% 
    mutate(targetTaxonGroup = ifelse(grepl("Fung", standardDescription), "fun", "bac_arc")) %>% 
    select(batchID, targetTaxonGroup, pcrEfficiency) %>% 
    distinct(.keep_all = TRUE)

qc_info2 <- qpcrData$mga_soilGroupAbundances %>% 
    filter(is.na(dataQF) | dataQF %in% c("legacyData")) %>% 
    select(dnaSampleID, batchID, nucleicAcidConcentration, meanCqValue, targetTaxonGroup) %>% 
    mutate(targetTaxonGroup = recode(targetTaxonGroup, 
                                     "bacteria and archaea" = "bac_arc",
                                     "bacteria" = "bac",
                                     "archaea" = "arc",
                                     "fungi" = "fun")) %>% 
    distinct(.keep_all = TRUE)

qc_info <- left_join(qc_info2, qc_info1, relationship = "many-to-many") %>% 
    distinct(.keep_all = TRUE) %>% 
    pivot_wider(names_from = "targetTaxonGroup", 
               values_from = c("nucleicAcidConcentration", "meanCqValue", "pcrEfficiency"), 
               values_fn = mean)

qpcrData_qc <- qpcrData$mga_soilGroupAbundances %>%
    select(siteID, plotID, collectDate, geneticSampleID, dnaSampleID, targetTaxonGroup, dataQF,
           meanCopyNumber, copyNumberStandardDeviation) %>%
    filter(is.na(dataQF) | dataQF %in% c("legacyData")) %>% 
    mutate(copyNumberVariance = copyNumberStandardDeviation^2)

bac_arc <- qpcrData_qc %>% 
    filter(targetTaxonGroup %in% c("bacteria", "archaea")) %>% 
    select(dnaSampleID, geneticSampleID, targetTaxonGroup, meanCopyNumber, copyNumberStandardDeviation, dataQF) %>% 
    pivot_wider(names_from = "targetTaxonGroup", 
               values_from = c("meanCopyNumber", "copyNumberStandardDeviation"))

non_bac_arc <- qpcrData_qc %>% 
    filter(targetTaxonGroup == "bacteria and archaea") %>% 
    select(dnaSampleID, geneticSampleID, targetTaxonGroup, meanCopyNumber, copyNumberStandardDeviation, dataQF) %>% 
    pivot_wider(names_from = "targetTaxonGroup", 
               values_from = c("meanCopyNumber", "copyNumberStandardDeviation"))

fung_only <- qpcrData_qc %>% 
    filter(targetTaxonGroup == "fungi") %>% 
    select(dnaSampleID, geneticSampleID, targetTaxonGroup, meanCopyNumber, copyNumberStandardDeviation, dataQF) %>% 
    pivot_wider(names_from = "targetTaxonGroup", 
               values_from = c("meanCopyNumber", "copyNumberStandardDeviation"))

combined_sd <- lapply(1:nrow(bac_arc), function(i) {
    grand.sd(S = c(bac_arc$copyNumberStandardDeviation_bacteria[i], 
                   bac_arc$copyNumberStandardDeviation_archaea[i]),
              M = c(bac_arc$meanCopyNumber_bacteria[i], 
                    bac_arc$meanCopyNumber_archaea[i]),
              N = c(3, 3))
})

bac_arc$`copyNumberStandardDeviation_bacteria and archaea` <- unlist(combined_sd)
bac_arc$`meanCopyNumber_bacteria and archaea` <- bac_arc$meanCopyNumber_bacteria + bac_arc$meanCopyNumber_archaea
bac_arc <- bac_arc %>% select(-copyNumberStandardDeviation_bacteria, -copyNumberStandardDeviation_archaea,
                               -meanCopyNumber_bacteria, -meanCopyNumber_archaea)

non_fungi_qpcrData <- rbindlist(list(non_bac_arc, bac_arc), use.names = TRUE)

qpcrData_out <- merge(fung_only, non_fungi_qpcrData, all = TRUE, by = c("geneticSampleID", "dnaSampleID")) %>% 
    mutate(legacy = ifelse(grepl("legacy", dataQF.x) | grepl("legacy", dataQF.y), TRUE, FALSE),
           fb_ratio = `meanCopyNumber_fungi` / `meanCopyNumber_bacteria and archaea`,
           f_proportion = `meanCopyNumber_fungi` / (`meanCopyNumber_fungi` + `meanCopyNumber_bacteria and archaea`),
           siteID = substr(dnaSampleID, 1, 4),
           test_fun = ifelse(legacy == TRUE, meanCopyNumber_fungi / 1000, meanCopyNumber_fungi),
           test_bac = ifelse(legacy == TRUE, `meanCopyNumber_bacteria and archaea` / 1000, 
                            `meanCopyNumber_bacteria and archaea`)) %>%
    left_join(qc_info)

source("scripts/helper_functions.r")
if(!exists("soilCores")) soilCores <- load_soilCores()
genSampleExample <- load_genSampleExample()

qpcrData_out$sampleID <- soilCores[match(qpcrData_out$geneticSampleID, soilCores$geneticSampleID),]$sampleID
qpcrData_out$compositeSampleID <- genSampleExample[match(qpcrData_out$sampleID, genSampleExample$sampleID),]$genomicsSampleID

dir.create(dirname(qpcr_output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(qpcrData_out, qpcr_output_file)
