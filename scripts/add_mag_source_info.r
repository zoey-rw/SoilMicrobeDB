# Add MAG source information to species-level abundances
# This relies on large data files associated with each database - cannot be run 

library(tidyverse)
library(ggpubr)
library(data.table)
library(ggpmisc)
library(rstatix)
source("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/helper_functions.r")
source("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/comets_shinyapp_example/source.R")


# Read in species-level abundances
bracken_with_lineage=fread("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_species_merged_lineage.csv", nThread = 8)

# Read in sequencing depth 
seq_depth_df <- readRDS("/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/seq_depth_df.rds") 

# Now read in all the Kraken and study-specific files
seqid_map_gtdb207 = fread("/projectnb/microbiome/ref_db/GTDB_207_kraken2/seqid2taxid.map")
seqid_map_gtdb207_filtered = fread("/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/seqid2taxid.map")

gtdb_struo = read_tsv("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/gtdb_207_filtered_struo.tsv")

# Label GTDB genomes as MAGs or not
gtdb_struo$is_MAG2 = 
    gtdb_207_metadata[match(gtdb_struo$ncbi_species_taxid,gtdb_207_metadata$ncbi_taxid),]$is_MAG
gtdb_struo$is_MAG = 
    gtdb_207_metadata[match(gtdb_struo$accession,gtdb_207_metadata$ncbi_genbank_assembly_accession),]$is_MAG

# Other struo files
spire_struo = read_tsv("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/spire_MAGs_struo.tsv") %>% 
    mutate(is_MAG=T)
SMAG_struo = read_tsv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/SMAG_struo.tsv") %>% 
    mutate(is_MAG=T)
nayfach_struo = read_tsv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/nayfach_MAGs_struo.tsv") %>% 
    mutate(source="GEM catalog") %>% mutate(is_MAG=T)
JGI_struo = read_tsv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/ncbi_struo.tsv")
myco_struo = read_tsv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/mycocosm_published_struo.tsv")

# Struo file for entire database
soil_microbe_db_struo = read_tsv("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/soil_genome_db_struo.tsv")

# Reference for PlusPF taxa
tax_lineage_key <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/merged_db_taxonomy.rds")

pluspf_tax_lineage_key = tax_lineage_key[tax_lineage_key$db_name=="PlusPF",]
pluspf_taxids = pluspf_tax_lineage_key$kraken_tax_id %>% unique

# Count for Table 1
soil_microbe_db_struo %>% group_by(kingdom, source) %>% tally()

ncbi_taxids = JGI_struo$ncbi_species_taxid %>% unique
gtdb_taxids = gtdb_struo$ncbi_species_taxid %>% unique

# Create list of all MAGs
MAG_ids = unique(c(spire_struo$ncbi_species_taxid, 
                   SMAG_struo$ncbi_species_taxid, 
                   nayfach_struo$ncbi_species_taxid,
                   gtdb_struo[which(gtdb_struo$is_MAG),]$ncbi_species_taxid))

# Create lists of MAGs unique to each source 
spire_unique <- setdiff(spire_struo$ncbi_species_taxid, gtdb_taxids) %>% unique 
smag_unique <-  setdiff(SMAG_struo$ncbi_species_taxid, gtdb_taxids) %>% unique 
nayfach_unique <-  setdiff(nayfach_struo$ncbi_species_taxid, gtdb_taxids) %>% unique 
myco_unique <-  setdiff(myco_struo$ncbi_species_taxid, pluspf_taxids) %>% unique 

spire_or_smag <- unique(c(spire_struo$ncbi_species_taxid, SMAG_struo$ncbi_species_taxid))

fungi_unique <-  soil_microbe_db_struo[which(soil_microbe_db_struo$kingdom=="Fungi",),]$ncbi_species_taxid

# Use lists to label abundance data
soil_microbe_db_struo$is_MAG = ifelse(soil_microbe_db_struo$ncbi_species_taxid %in% MAG_ids, T, F)

write.csv(soil_microbe_db_struo %>% select(-c(fasta_file_path, tax_id, tax_name)), "/projectnb/frpmars/soil_microbe_db/data/soil_microbe_db_genome_table.csv")

# Sanity check - get counts
table(soil_microbe_db_struo$is_MAG)/nrow(soil_microbe_db_struo)

# Use specific MAG lists to assign categories to abundances
bracken_with_lineage$smag_unique = ifelse(bracken_with_lineage$taxonomy_id %in% smag_unique, T, F)
bracken_with_lineage$spire_unique = ifelse(bracken_with_lineage$taxonomy_id %in% spire_unique, T, F)
bracken_with_lineage$nayfach_unique = ifelse(bracken_with_lineage$taxonomy_id %in% nayfach_unique, T, F)
bracken_with_lineage$is_fungi = ifelse(bracken_with_lineage$taxonomy_id %in% fungi_unique, T, F)
bracken_with_lineage$myco_unique = ifelse(bracken_with_lineage$taxonomy_id %in% myco_unique, T, F)
bracken_with_lineage$is_MAG = ifelse(bracken_with_lineage$taxonomy_id %in% MAG_ids, T, F)
bracken_with_lineage$spire_or_smag = ifelse(bracken_with_lineage$taxonomy_id %in% spire_or_smag, T, F)
bracken_with_lineage = bracken_with_lineage %>% mutate(novel_mag = ifelse(smag_unique == T | 
                                                                              spire_unique == T | 
                                                                              nayfach_unique == T, T, F),
                                                       novel_fungi = ifelse(myco_unique == T, T, F))

# Label source using Struo input TSV
bracken_with_lineage$source = soil_microbe_db_struo[match(bracken_with_lineage$taxonomy_id,
                                                          soil_microbe_db_struo$ncbi_species_taxid),]$source

# Save abundances for visualization portal
species_abundances = bracken_with_lineage %>% select(sample_id, name, taxonomy_id, percentage = fraction_total_reads,
                                                     lineage, source, is_MAG, taxid_lineage)
write.csv(species_abundances, "/projectnb/frpmars/soil_microbe_db/shiny_app/soil_microbe_db_abundances.csv")



# Now add in biome information for summary stats - only for SMDB 
soilData_subset = soilCores %>% 
    mutate(sample_id = paste0(compositeSampleID, "_soil_microbe_db_filtered")) %>% 
    filter(sample_id %in% bracken_with_lineage$sample_id) %>% 
    select(sample_id, compositeSampleID, siteID, biome, horizon, nlcdClass, elevation, decimalLatitude) %>% 
    distinct(compositeSampleID, .keep_all = T) %>% 
    filter(!is.na(compositeSampleID))

bracken_with_lineage_biome = left_join(bracken_with_lineage,
                                       soilData_subset) 

bracken_with_lineage_biome$is_alaska = ifelse(bracken_with_lineage_biome$siteID %in% alaska_sites, 
                                              "Alaska site", "Not Alaska site")
bracken_with_lineage_biome = bracken_with_lineage_biome %>% 
    mutate(custom_biome = 
               ifelse(siteID %in% alaska_sites,
                      paste0(biome, " Alaska"), biome))


weird = c("OSBS_005-O-20200706-COMP_soil_microbe_db_filtered",
          "UNDE_043-M-20180724-COMP_soil_microbe_db_filtered",
          "BLAN_033-M-20170714-COMP_soil_microbe_db_filtered",
          "TALL_004-M-20140708-COMP_soil_microbe_db_filtered")

# These files only have a couple lines - corrupted
bracken_with_lineage_biome <- bracken_with_lineage_biome %>% filter(!sample_id %in% weird)

# Calculate alpha diversity (# of unique species per sample, for MAGs and non-MAGs)
species_alpha_diversity = bracken_with_lineage_biome %>% 
    group_by(custom_biome, biome,decimalLatitude, nlcdClass, siteID, is_alaska,sample_id, is_MAG) %>% 
    distinct(taxonomy_id, .keep_all = T) %>% tally(name = "alpha_div")

# Summarize for plotting
is_novel_mag_df = bracken_with_lineage_biome %>%
    group_by(custom_biome, biome,decimalLatitude, nlcdClass, siteID, is_alaska,sample_id, novel_mag) %>% 
    summarize(sum = sum(fraction_total_reads)) %>% filter(novel_mag==T)
is_fungi_df = bracken_with_lineage_biome %>% 
    group_by(custom_biome, biome, decimalLatitude,nlcdClass, siteID, is_alaska,sample_id, novel_fungi) %>% 
    summarize(sum = sum(fraction_total_reads)) %>% filter(novel_fungi==T)
is_mag_df = bracken_with_lineage_biome %>% 
    group_by(custom_biome, biome, decimalLatitude,nlcdClass, siteID, is_alaska,sample_id, is_MAG) %>%
    summarize(sum = sum(fraction_total_reads))  %>% filter(is_MAG==T)


# Run beginning of /projectnb/frpmars/soil_microbe_db/scripts/compare_unknown_by_biome_species.r
is_mag_df = left_join(is_mag_df,pass_filter_species_long %>% filter(metric == "percent_classified"))

seq_depth_df$percent_classified2 = seq_depth_df$identified_reads/seq_depth_df$seq_depth
is_mag_df = left_join(is_mag_df,seq_depth_df)  %>% ungroup

write.csv(is_mag_df, "/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files/abundance_MAGs.csv")
