# Add MAG (Metagenome-Assembled Genome) source information to species-level abundances
# Input: Species-level merged lineage CSV from soil_microbe_db
#        Struo2 input TSV files for each genome database (GTDB, SPIRE, SMAG, Nayfach, JGI, Mycocosm)
#        Merged taxonomy RDS (may be on external HARDDRIVE)
#        Requires gtdb_207_metadata from helper_functions.r
# Output: abundance_MAGs.csv with MAG abundance summaries by biome
#         table1_genome_counts.csv with genome counts by kingdom and source
#         soil_microbe_db_genome_table.csv with genome metadata

library(tidyverse)
library(data.table)
source("scripts/helper_functions.r")

struo_dir <- "data/genome_database/struo2_input_tables"
bracken_file <- "data/classification/taxonomic_rank_summaries/soil_microbe_db_species_merged_lineage.csv"
seq_depth_file <- "data/classification/analysis_files/seq_depth_df.rds"
pass_filter_file <- "data/classification/analysis_files/pass_filter_summary.csv"

tax_lineage_file <- ifelse(file.exists("data/genome_database/merged_db_taxonomy.rds"),
                          "data/genome_database/merged_db_taxonomy.rds",
                          "/Volumes/HARDDRIVE/genome_databases/merged_db_taxonomy.rds")

if(!file.exists(tax_lineage_file)) {
    stop("merged_db_taxonomy.rds not found")
}

if(!exists("soilCores")) soilCores <- load_soilCores()
if(!exists("alaska_sites")) alaska_sites <- c("BONA", "DEJU", "HEAL", "TOOL", "BARR")
if(!exists("gtdb_207_metadata")) {
    stop("gtdb_207_metadata not found - ensure helper_functions.r is sourced correctly")
}

bracken_with_lineage <- fread(bracken_file, nThread = 8)
seq_depth_df <- readRDS(seq_depth_file)

gtdb_struo <- read_tsv(file.path(struo_dir, "gtdb_207_filtered_struo.tsv"))
spire_struo <- read_tsv(file.path(struo_dir, "spire_MAGs_struo.tsv")) %>% mutate(is_MAG = TRUE)
SMAG_struo <- read_tsv(file.path(struo_dir, "SMAG_struo.tsv")) %>% mutate(is_MAG = TRUE)
nayfach_struo <- read_tsv(file.path(struo_dir, "nayfach_MAGs_struo.tsv")) %>% 
    mutate(source = "GEM catalog", is_MAG = TRUE)
JGI_struo <- read_tsv(file.path(struo_dir, "ncbi_struo.tsv"))
myco_struo <- read_tsv(file.path(struo_dir, "mycocosm_published_struo.tsv"))
soil_microbe_db_struo <- read_tsv(file.path(struo_dir, "soil_genome_db_struo.tsv"))

gtdb_struo$is_MAG2 <- gtdb_207_metadata[match(gtdb_struo$ncbi_species_taxid, gtdb_207_metadata$ncbi_taxid),]$is_MAG
gtdb_struo$is_MAG <- gtdb_207_metadata[match(gtdb_struo$accession, gtdb_207_metadata$ncbi_genbank_assembly_accession),]$is_MAG

tax_lineage_key <- readRDS(tax_lineage_file)
pluspf_taxids <- tax_lineage_key[tax_lineage_key$db_name == "PlusPF",]$kraken_tax_id %>% unique

table1_counts <- soil_microbe_db_struo %>% group_by(kingdom, source) %>% tally()
dir.create("manuscript_figures", recursive = TRUE, showWarnings = FALSE)
write.csv(table1_counts, "manuscript_figures/table1_genome_counts.csv", row.names = FALSE)

ncbi_taxids <- JGI_struo$ncbi_species_taxid %>% unique
gtdb_taxids <- gtdb_struo$ncbi_species_taxid %>% unique

MAG_ids <- unique(c(spire_struo$ncbi_species_taxid, 
                   SMAG_struo$ncbi_species_taxid, 
                   nayfach_struo$ncbi_species_taxid,
                   gtdb_struo[which(gtdb_struo$is_MAG),]$ncbi_species_taxid))

spire_unique <- setdiff(spire_struo$ncbi_species_taxid, gtdb_taxids) %>% unique 
smag_unique <- setdiff(SMAG_struo$ncbi_species_taxid, gtdb_taxids) %>% unique 
nayfach_unique <- setdiff(nayfach_struo$ncbi_species_taxid, gtdb_taxids) %>% unique 
myco_unique <- setdiff(myco_struo$ncbi_species_taxid, pluspf_taxids) %>% unique 
spire_or_smag <- unique(c(spire_struo$ncbi_species_taxid, SMAG_struo$ncbi_species_taxid))
fungi_unique <- soil_microbe_db_struo[which(soil_microbe_db_struo$kingdom == "Fungi",),]$ncbi_species_taxid

soil_microbe_db_struo$is_MAG <- ifelse(soil_microbe_db_struo$ncbi_species_taxid %in% MAG_ids, TRUE, FALSE)

dir.create("data/genome_database", recursive = TRUE, showWarnings = FALSE)
write.csv(soil_microbe_db_struo %>% select(-c(fasta_file_path, tax_id, tax_name)), 
          "data/genome_database/soil_microbe_db_genome_table.csv", row.names = FALSE)

bracken_with_lineage$smag_unique <- bracken_with_lineage$taxonomy_id %in% smag_unique
bracken_with_lineage$spire_unique <- bracken_with_lineage$taxonomy_id %in% spire_unique
bracken_with_lineage$nayfach_unique <- bracken_with_lineage$taxonomy_id %in% nayfach_unique
bracken_with_lineage$is_fungi <- bracken_with_lineage$taxonomy_id %in% fungi_unique
bracken_with_lineage$myco_unique <- bracken_with_lineage$taxonomy_id %in% myco_unique
bracken_with_lineage$is_MAG <- bracken_with_lineage$taxonomy_id %in% MAG_ids
bracken_with_lineage$spire_or_smag <- bracken_with_lineage$taxonomy_id %in% spire_or_smag
bracken_with_lineage <- bracken_with_lineage %>% 
    mutate(novel_mag = smag_unique | spire_unique | nayfach_unique,
           novel_fungi = myco_unique)

bracken_with_lineage$source <- soil_microbe_db_struo[match(bracken_with_lineage$taxonomy_id,
                                                          soil_microbe_db_struo$ncbi_species_taxid),]$source

soilData_subset <- soilCores %>% 
    mutate(sample_id = paste0(compositeSampleID, "_soil_microbe_db_filtered")) %>% 
    filter(sample_id %in% bracken_with_lineage$sample_id) %>% 
    select(sample_id, compositeSampleID, siteID, biome, horizon, nlcdClass, elevation, decimalLatitude) %>% 
    distinct(sample_id, .keep_all = TRUE) %>% 
    filter(!is.na(compositeSampleID))

bracken_with_lineage_biome <- left_join(bracken_with_lineage, soilData_subset, by = "sample_id") %>%
    mutate(is_alaska = ifelse(siteID %in% alaska_sites, "Alaska site", "Not Alaska site"),
           custom_biome = ifelse(siteID %in% alaska_sites, paste0(biome, " Alaska"), biome)) %>%
    filter(!sample_id %in% c("OSBS_005-O-20200706-COMP_soil_microbe_db_filtered",
                             "UNDE_043-M-20180724-COMP_soil_microbe_db_filtered",
                             "BLAN_033-M-20170714-COMP_soil_microbe_db_filtered",
                             "TALL_004-M-20140708-COMP_soil_microbe_db_filtered"))

is_mag_df <- bracken_with_lineage_biome %>% 
    group_by(custom_biome, biome, decimalLatitude, nlcdClass, siteID, is_alaska, sample_id, is_MAG) %>%
    summarize(sum = sum(fraction_total_reads), .groups = "drop") %>% 
    filter(is_MAG == TRUE) %>%
    left_join(read_csv(pass_filter_file) %>% filter(metric == "percent_classified")) %>%
    left_join(seq_depth_df %>% mutate(percent_classified2 = identified_reads / seq_depth)) %>%
    ungroup()

write.csv(is_mag_df, "data/classification/analysis_files/abundance_MAGs.csv", row.names = FALSE)
