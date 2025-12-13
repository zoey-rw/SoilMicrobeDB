# Add MAG source information to species-level abundances
# This relies on large data files associated with each database - cannot be run unless all files are present

library(tidyverse)
library(ggpubr)
library(data.table)
library(ggpmisc)
library(rstatix)
source("scripts/helper_functions.r")

# Define file paths
struo_dir <- "data/genome_database/struo2_input_tables"
bracken_file <- "data/classification/taxonomic_rank_summaries/species/soil_microbe_db_filtered_species_merged_lineage.csv"
seq_depth_file <- "data/classification/analysis_files/seq_depth_df.rds"
pass_filter_file <- "data/classification/analysis_files/pass_filter_summary.csv"

# Taxonomy file may be on external hard drive (large file) - check both locations
tax_lineage_file_local <- "data/genome_database/merged_db_taxonomy.rds"
tax_lineage_file_external <- "/Volumes/HARDDRIVE/genome_databases/merged_db_taxonomy.rds"

tax_lineage_file <- NULL
if(file.exists(tax_lineage_file_local)) {
    tax_lineage_file <- tax_lineage_file_local
    cat("✓ Using taxonomy file (local):", tax_lineage_file_local, "\n")
} else if(file.exists(tax_lineage_file_external)) {
    tax_lineage_file <- tax_lineage_file_external
    cat("✓ Using taxonomy file (external):", tax_lineage_file_external, "\n")
} else {
    stop("❌ MISSING FILE: merged_db_taxonomy.rds not found!\n",
         "   Checked: ", tax_lineage_file_local, "\n",
         "   Checked: ", tax_lineage_file_external, "\n",
         "   Please ensure the file exists in one of these locations.")
}

# Load soilCores using helper function
if(!exists("soilCores")) {
    soilCores <- load_soilCores()
}

# gtdb_207_metadata comes from helper_functions.r (already sourced above)
# alaska_sites needs to be defined - define it here if not already loaded from source.R
if(!exists("alaska_sites")) {
    alaska_sites <- c("BONA","DEJU","HEAL","TOOL","BARR")
}

# Check for required data objects
if(!exists("gtdb_207_metadata")) {
    stop("❌ MISSING DATA OBJECT: gtdb_207_metadata\n",
         "   This should be loaded from scripts/helper_functions.r\n",
         "   Please ensure helper_functions.r is sourced correctly.")
}

# Read in species-level abundances
bracken_with_lineage <- fread(bracken_file, nThread = 8)

# Read in sequencing depth 
seq_depth_df <- readRDS(seq_depth_file)

# Read in all Struo files
gtdb_struo <- read_tsv(file.path(struo_dir, "gtdb_207_filtered_struo.tsv"))
spire_struo <- read_tsv(file.path(struo_dir, "spire_MAGs_struo.tsv")) %>% mutate(is_MAG=T)
SMAG_struo <- read_tsv(file.path(struo_dir, "SMAG_struo.tsv")) %>% mutate(is_MAG=T)
nayfach_struo <- read_tsv(file.path(struo_dir, "nayfach_MAGs_struo.tsv")) %>% 
    mutate(source="GEM catalog") %>% mutate(is_MAG=T)
JGI_struo <- read_tsv(file.path(struo_dir, "ncbi_struo.tsv"))
myco_struo <- read_tsv(file.path(struo_dir, "mycocosm_published_struo.tsv"))
soil_microbe_db_struo <- read_tsv(file.path(struo_dir, "soil_genome_db_struo.tsv"))

# Label GTDB genomes as MAGs or not
gtdb_struo$is_MAG2 <- 
    gtdb_207_metadata[match(gtdb_struo$ncbi_species_taxid, gtdb_207_metadata$ncbi_taxid),]$is_MAG
gtdb_struo$is_MAG <- 
    gtdb_207_metadata[match(gtdb_struo$accession, gtdb_207_metadata$ncbi_genbank_assembly_accession),]$is_MAG

# Reference for PlusPF taxa
tax_lineage_key <- readRDS(tax_lineage_file)
pluspf_tax_lineage_key <- tax_lineage_key[tax_lineage_key$db_name=="PlusPF",]
pluspf_taxids <- pluspf_tax_lineage_key$kraken_tax_id %>% unique

# Count for Table 1
table1_counts <- soil_microbe_db_struo %>% group_by(kingdom, source) %>% tally()
write.csv(table1_counts, "manuscript_figures/table1_genome_counts.csv", row.names = FALSE)

ncbi_taxids <- JGI_struo$ncbi_species_taxid %>% unique
gtdb_taxids <- gtdb_struo$ncbi_species_taxid %>% unique

# Create list of all MAGs
MAG_ids <- unique(c(spire_struo$ncbi_species_taxid, 
                   SMAG_struo$ncbi_species_taxid, 
                   nayfach_struo$ncbi_species_taxid,
                   gtdb_struo[which(gtdb_struo$is_MAG),]$ncbi_species_taxid))

# Create lists of MAGs unique to each source 
spire_unique <- setdiff(spire_struo$ncbi_species_taxid, gtdb_taxids) %>% unique 
smag_unique <- setdiff(SMAG_struo$ncbi_species_taxid, gtdb_taxids) %>% unique 
nayfach_unique <- setdiff(nayfach_struo$ncbi_species_taxid, gtdb_taxids) %>% unique 
myco_unique <- setdiff(myco_struo$ncbi_species_taxid, pluspf_taxids) %>% unique 

spire_or_smag <- unique(c(spire_struo$ncbi_species_taxid, SMAG_struo$ncbi_species_taxid))

fungi_unique <- soil_microbe_db_struo[which(soil_microbe_db_struo$kingdom=="Fungi",),]$ncbi_species_taxid

# Use lists to label abundance data
soil_microbe_db_struo$is_MAG <- ifelse(soil_microbe_db_struo$ncbi_species_taxid %in% MAG_ids, T, F)

write.csv(soil_microbe_db_struo %>% select(-c(fasta_file_path, tax_id, tax_name)), 
          "data/genome_database/soil_microbe_db_genome_table.csv")

# Sanity check - get counts
table(soil_microbe_db_struo$is_MAG)/nrow(soil_microbe_db_struo)

# Use specific MAG lists to assign categories to abundances
bracken_with_lineage$smag_unique <- ifelse(bracken_with_lineage$taxonomy_id %in% smag_unique, T, F)
bracken_with_lineage$spire_unique <- ifelse(bracken_with_lineage$taxonomy_id %in% spire_unique, T, F)
bracken_with_lineage$nayfach_unique <- ifelse(bracken_with_lineage$taxonomy_id %in% nayfach_unique, T, F)
bracken_with_lineage$is_fungi <- ifelse(bracken_with_lineage$taxonomy_id %in% fungi_unique, T, F)
bracken_with_lineage$myco_unique <- ifelse(bracken_with_lineage$taxonomy_id %in% myco_unique, T, F)
bracken_with_lineage$is_MAG <- ifelse(bracken_with_lineage$taxonomy_id %in% MAG_ids, T, F)
bracken_with_lineage$spire_or_smag <- ifelse(bracken_with_lineage$taxonomy_id %in% spire_or_smag, T, F)
bracken_with_lineage <- bracken_with_lineage %>% mutate(novel_mag = ifelse(smag_unique == T | 
                                                                              spire_unique == T | 
                                                                              nayfach_unique == T, T, F),
                                                       novel_fungi = ifelse(myco_unique == T, T, F))

# Label source using Struo input TSV
bracken_with_lineage$source <- soil_microbe_db_struo[match(bracken_with_lineage$taxonomy_id,
                                                          soil_microbe_db_struo$ncbi_species_taxid),]$source

# Save abundances for visualization portal
species_abundances <- bracken_with_lineage %>% select(sample_id, name, taxonomy_id, percentage = fraction_total_reads,
                                                     lineage, source, is_MAG, taxid_lineage)
# Save abundances for visualization portal (external path - commented out for local use)
# write.csv(species_abundances, "/projectnb/frpmars/soil_microbe_db/shiny_app/soil_microbe_db_abundances.csv")

# Now add in biome information for summary stats - only for SMDB 
soilData_subset <- soilCores %>% 
    mutate(sample_id = paste0(compositeSampleID, "_soil_microbe_db_filtered")) %>% 
    filter(sample_id %in% bracken_with_lineage$sample_id) %>% 
    select(sample_id, compositeSampleID, siteID, biome, horizon, nlcdClass, elevation, decimalLatitude) %>% 
    distinct(sample_id, .keep_all = T) %>% 
    filter(!is.na(compositeSampleID))

bracken_with_lineage_biome <- left_join(bracken_with_lineage, soilData_subset, by = "sample_id") 

bracken_with_lineage_biome$is_alaska <- ifelse(bracken_with_lineage_biome$siteID %in% alaska_sites, 
                                              "Alaska site", "Not Alaska site")
bracken_with_lineage_biome <- bracken_with_lineage_biome %>% 
    mutate(custom_biome = 
               ifelse(siteID %in% alaska_sites,
                      paste0(biome, " Alaska"), biome))

weird <- c("OSBS_005-O-20200706-COMP_soil_microbe_db_filtered",
          "UNDE_043-M-20180724-COMP_soil_microbe_db_filtered",
          "BLAN_033-M-20170714-COMP_soil_microbe_db_filtered",
          "TALL_004-M-20140708-COMP_soil_microbe_db_filtered")

# These files only have a couple lines - corrupted
bracken_with_lineage_biome <- bracken_with_lineage_biome %>% filter(!sample_id %in% weird)

# Calculate alpha diversity (# of unique species per sample, for MAGs and non-MAGs)
species_alpha_diversity <- bracken_with_lineage_biome %>% 
    group_by(custom_biome, biome, decimalLatitude, nlcdClass, siteID, is_alaska, sample_id, is_MAG) %>% 
    distinct(taxonomy_id, .keep_all = T) %>% tally(name = "alpha_div")

# Summarize for plotting
is_novel_mag_df <- bracken_with_lineage_biome %>%
    group_by(custom_biome, biome, decimalLatitude, nlcdClass, siteID, is_alaska, sample_id, novel_mag) %>% 
    summarize(sum = sum(fraction_total_reads)) %>% filter(novel_mag==T)
is_fungi_df <- bracken_with_lineage_biome %>% 
    group_by(custom_biome, biome, decimalLatitude, nlcdClass, siteID, is_alaska, sample_id, novel_fungi) %>% 
    summarize(sum = sum(fraction_total_reads)) %>% filter(novel_fungi==T)
is_mag_df <- bracken_with_lineage_biome %>% 
    group_by(custom_biome, biome, decimalLatitude, nlcdClass, siteID, is_alaska, sample_id, is_MAG) %>%
    summarize(sum = sum(fraction_total_reads)) %>% filter(is_MAG==T)

# Join percent_classified from pass_filter_summary
is_mag_df <- left_join(is_mag_df, 
                       read_csv(pass_filter_file) %>% 
                           filter(metric == "percent_classified"))

seq_depth_df$percent_classified2 <- seq_depth_df$identified_reads/seq_depth_df$seq_depth
is_mag_df <- left_join(is_mag_df, seq_depth_df) %>% ungroup

write.csv(is_mag_df, "data/classification/analysis_files/abundance_MAGs.csv")

cat("Completed add_mag_source_info.r\n")
cat("   Saved: manuscript_figures/table1_genome_counts.csv\n")
cat("   Saved: data/genome_database/soil_microbe_db_genome_table.csv\n")
cat("   Saved: data/classification/analysis_files/abundance_MAGs.csv\n")
