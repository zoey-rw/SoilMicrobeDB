# Add MAG source information to species-level abundances
# This relies on large data files associated with each database - cannot be run 

library(tidyverse)
library(ggpubr)
library(data.table)
library(ggpmisc)
library(rstatix)
source("scripts/helper_functions.r")
# Note: source.R from comets_shinyapp_example may not be available locally - commented out
# source("/projectnb/talbot-lab-data/zrwerbin/soil_microbe_GEMs/comets_shinyapp_example/source.R")


# Read in species-level abundances
bracken_with_lineage=fread("data/classification/taxonomic_rank_summaries/species/soil_microbe_db_filtered_species_merged_lineage.csv", nThread = 8)

# Read in sequencing depth 
seq_depth_df <- readRDS("data/classification/analysis_files/seq_depth_df.rds") 

# Now read in all the Kraken and study-specific files
# Note: These are external files from cluster - Struo files contain genome metadata
struo_dir <- "data/genome_database/struo_files"
if(!dir.exists(struo_dir)) {
    warning("⚠️  MISSING DIRECTORY: Struo files directory not found!")
    warning("   Expected: ", struo_dir)
    warning("   Struo files contain genome metadata needed for MAG source analysis.")
    warning("   Continuing without Struo files - some analyses may be incomplete.")
    # Create empty data frames to prevent errors
    gtdb_struo <- data.frame()
    spire_struo <- data.frame()
    SMAG_struo <- data.frame()
    nayfach_struo <- data.frame()
    JGI_struo <- data.frame()
    myco_struo <- data.frame()
    soil_microbe_db_struo <- data.frame()
} else {
    # Try to load Struo files if directory exists
    gtdb_struo_file <- file.path(struo_dir, "gtdb_207_filtered_struo.tsv")
    if(file.exists(gtdb_struo_file)) {
        gtdb_struo = read_tsv(gtdb_struo_file)
    } else {
        warning("⚠️  MISSING FILE: ", gtdb_struo_file)
        gtdb_struo <- data.frame()
    }
    
    spire_struo_file <- file.path(struo_dir, "spire_MAGs_struo.tsv")
    if(file.exists(spire_struo_file)) {
        spire_struo = read_tsv(spire_struo_file) %>% mutate(is_MAG=T)
    } else {
        warning("⚠️  MISSING FILE: ", spire_struo_file)
        spire_struo <- data.frame()
    }
    
    SMAG_struo_file <- file.path(struo_dir, "SMAG_struo.tsv")
    if(file.exists(SMAG_struo_file)) {
        SMAG_struo = read_tsv(SMAG_struo_file) %>% mutate(is_MAG=T)
    } else {
        warning("⚠️  MISSING FILE: ", SMAG_struo_file)
        SMAG_struo <- data.frame()
    }
    
    nayfach_struo_file <- file.path(struo_dir, "nayfach_MAGs_struo.tsv")
    if(file.exists(nayfach_struo_file)) {
        nayfach_struo = read_tsv(nayfach_struo_file) %>% mutate(source="GEM catalog") %>% mutate(is_MAG=T)
    } else {
        warning("⚠️  MISSING FILE: ", nayfach_struo_file)
        nayfach_struo <- data.frame()
    }
    
    JGI_struo_file <- file.path(struo_dir, "ncbi_struo.tsv")
    if(file.exists(JGI_struo_file)) {
        JGI_struo = read_tsv(JGI_struo_file)
    } else {
        warning("⚠️  MISSING FILE: ", JGI_struo_file)
        JGI_struo <- data.frame()
    }
    
    myco_struo_file <- file.path(struo_dir, "mycocosm_published_struo.tsv")
    if(file.exists(myco_struo_file)) {
        myco_struo = read_tsv(myco_struo_file)
    } else {
        warning("⚠️  MISSING FILE: ", myco_struo_file)
        myco_struo <- data.frame()
    }
    
    soil_microbe_db_struo_file <- file.path(struo_dir, "soil_genome_db_struo.tsv")
    if(file.exists(soil_microbe_db_struo_file)) {
        soil_microbe_db_struo = read_tsv(soil_microbe_db_struo_file)
    } else {
        warning("⚠️  MISSING FILE: ", soil_microbe_db_struo_file)
        soil_microbe_db_struo <- data.frame()
    }
}

# External seqid map files (from cluster) - commented out as they're not needed for local analysis
# seqid_map_gtdb207 = fread("/projectnb/microbiome/ref_db/GTDB_207_kraken2/seqid2taxid.map")
# seqid_map_gtdb207_filtered = fread("/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/seqid2taxid.map")

# Label GTDB genomes as MAGs or not (only if gtdb_struo was loaded)
if(nrow(gtdb_struo) > 0 && exists("gtdb_207_metadata")) {
    gtdb_struo$is_MAG2 = 
        gtdb_207_metadata[match(gtdb_struo$ncbi_species_taxid,gtdb_207_metadata$ncbi_taxid),]$is_MAG
    gtdb_struo$is_MAG = 
        gtdb_207_metadata[match(gtdb_struo$accession,gtdb_207_metadata$ncbi_genbank_assembly_accession),]$is_MAG
}

# Reference for PlusPF taxa
tax_lineage_file <- "data/genome_database/merged_db_taxonomy.rds"
if(file.exists(tax_lineage_file)) {
    tax_lineage_key <- readRDS(tax_lineage_file)
    pluspf_tax_lineage_key = tax_lineage_key[tax_lineage_key$db_name=="PlusPF",]
    pluspf_taxids = pluspf_tax_lineage_key$kraken_tax_id %>% unique
} else {
    warning("⚠️  MISSING FILE: Taxonomy lineage key not found!")
    warning("   Expected: ", tax_lineage_file)
    warning("   Continuing without PlusPF taxonomy reference")
    tax_lineage_key <- data.frame()
    pluspf_taxids <- numeric(0)
}

# Count for Table 1 (only if Struo files were loaded)
if(nrow(soil_microbe_db_struo) > 0) {
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
    
    write.csv(soil_microbe_db_struo %>% select(-c(fasta_file_path, tax_id, tax_name)), "data/genome_database/soil_microbe_db_genome_table.csv")
    
    # Sanity check - get counts
    table(soil_microbe_db_struo$is_MAG)/nrow(soil_microbe_db_struo)
} else {
    # If Struo files are missing, create empty variables to prevent errors
    ncbi_taxids <- numeric(0)
    gtdb_taxids <- numeric(0)
    MAG_ids <- numeric(0)
    spire_unique <- numeric(0)
    smag_unique <- numeric(0)
    nayfach_unique <- numeric(0)
    myco_unique <- numeric(0)
    spire_or_smag <- numeric(0)
    fungi_unique <- numeric(0)
}

# Use specific MAG lists to assign categories to abundances (only if Struo files were loaded)
if(nrow(soil_microbe_db_struo) > 0) {
    bracken_with_lineage$smag_unique = ifelse(bracken_with_lineage$taxonomy_id %in% smag_unique, T, F)
    bracken_with_lineage$spire_unique = ifelse(bracken_with_lineage$taxonomy_id %in% spire_unique, T, F)
    bracken_with_lineage$nayfach_unique = ifelse(bracken_with_lineage$taxonomy_id %in% nayfach_unique, T, F)
    bracken_with_lineage$is_fungi = ifelse(bracken_with_lineage$taxonomy_id %in% fungi_unique, T, F)
    bracken_with_lineage$myco_unique = ifelse(bracken_with_lineage$taxonomy_id %in% myco_unique, T, F)
    bracken_with_lineage$is_MAG = ifelse(bracken_with_lineage$taxonomy_id %in% MAG_ids, T, F)
    bracken_with_lineage$spire_or_smag = ifelse(bracken_with_lineage$taxonomy_id %in% spire_or_smag, T, F)
} else {
    # Set defaults when Struo files are missing
    bracken_with_lineage$smag_unique = F
    bracken_with_lineage$spire_unique = F
    bracken_with_lineage$nayfach_unique = F
    bracken_with_lineage$is_fungi = F
    bracken_with_lineage$myco_unique = F
    bracken_with_lineage$is_MAG = F
    bracken_with_lineage$spire_or_smag = F
}
bracken_with_lineage = bracken_with_lineage %>% mutate(novel_mag = ifelse(smag_unique == T | 
                                                                              spire_unique == T | 
                                                                              nayfach_unique == T, T, F),
                                                       novel_fungi = ifelse(myco_unique == T, T, F))

# Label source using Struo input TSV (only if Struo files were loaded)
if(nrow(soil_microbe_db_struo) > 0) {
    bracken_with_lineage$source = soil_microbe_db_struo[match(bracken_with_lineage$taxonomy_id,
                                                              soil_microbe_db_struo$ncbi_species_taxid),]$source
} else {
    bracken_with_lineage$source = NA
}

# Save abundances for visualization portal
species_abundances = bracken_with_lineage %>% select(sample_id, name, taxonomy_id, percentage = fraction_total_reads,
                                                     lineage, source, is_MAG, taxid_lineage)
# Save abundances for visualization portal (external path - commented out for local use)
# write.csv(species_abundances, "/projectnb/frpmars/soil_microbe_db/shiny_app/soil_microbe_db_abundances.csv")



# Now add in biome information for summary stats - only for SMDB 
# Expected file: data/comparison_data/neon/neon_soil_data_2023.rds
# Or load via: soilData <- readRDS("data/comparison_data/neon/neon_soil_data_2023.rds"); soilCores <- soilData$sls_soilCoreCollection
if(!exists("soilCores")) {
    stop("❌ MISSING DATA: soilCores not found!\n",
         "   Expected: data/comparison_data/neon/neon_soil_data_2023.rds\n",
         "   Load with: soilData <- readRDS('data/comparison_data/neon/neon_soil_data_2023.rds'); soilCores <- soilData$sls_soilCoreCollection\n",
         "   This is required for adding site metadata to MAG source analysis.")
}
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

write.csv(is_mag_df, "data/classification/analysis_files/abundance_MAGs.csv")
