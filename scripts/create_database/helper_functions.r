pacman::p_load(data.table, tidyverse, tidyfast, phyloseq, CHNOSZ)

# ============================================================================
# Configuration Helper Functions
# ============================================================================
# These functions are used by config.R and processing scripts

# Source-specific directory function
get_source_dir <- function(source_name) {
  if (!exists("GENOME_DB_DIR")) {
    stop("GENOME_DB_DIR not defined. Please source config.R first.")
  }
  file.path(GENOME_DB_DIR, source_name)
}

# Source-specific log directory function
get_source_log_dir <- function(source_name) {
  if (!exists("CREATE_DB_DIR")) {
    stop("CREATE_DB_DIR not defined. Please source config.R first.")
  }
  log_dir <- file.path(CREATE_DB_DIR, "process_genomes", source_name, "logs")
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  return(log_dir)
}

# Logging function
log_message <- function(msg, source_name = "general") {
  if (!exists("LOG_BASE_DIR")) {
    # Fallback if config.R not sourced
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"))
    return(invisible(NULL))
  }
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_dir <- if (source_name != "general") {
    get_source_log_dir(source_name)
  } else {
    LOG_BASE_DIR
  }
  log_file <- file.path(log_dir, paste0(source_name, "_", format(Sys.Date(), "%Y%m%d"), ".log"))
  message <- paste0("[", timestamp, "] ", msg)
  cat(message, "\n")
  write(message, file = log_file, append = TRUE)
}

# Check if in test mode
is_test_mode <- function() {
  if (!exists("TEST_MODE") || !exists("MAX_GENOMES")) {
    return(FALSE)
  }
  return(TEST_MODE || !is.na(MAX_GENOMES))
}

# Helper function to get source config
get_source_config <- function(source_name) {
  if (!exists("STRUO2_INPUT_DIR")) {
    stop("STRUO2_INPUT_DIR not defined. Please source config.R first.")
  }
  config_name <- paste0(toupper(source_name), "_CONFIG")
  if (exists(config_name)) {
    return(get(config_name))
  } else {
    # Return default config structure
    return(list(
      output_file = file.path(STRUO2_INPUT_DIR, paste0(tolower(source_name), "_struo.tsv"))
    ))
  }
}

# Print configuration summary (useful for debugging)
print_config <- function() {
  if (!exists("BASE_DIR")) {
    cat("Configuration not loaded. Please source config.R first.\n")
    return(invisible(NULL))
  }
  cat("=== Database Creation Configuration ===\n")
  cat("Base Directory:", BASE_DIR, "\n")
  cat("Genome Database Directory:", GENOME_DB_DIR, "\n")
  cat("Struo2 Input Directory:", STRUO2_INPUT_DIR, "\n")
  cat("Test Mode:", TEST_MODE, "\n")
  if (!is.na(MAX_GENOMES)) {
    cat("Max Genomes (test):", MAX_GENOMES, "\n")
  }
  cat("Min Completeness:", MIN_COMPLETENESS, "%\n")
  cat("Max Contamination:", MAX_CONTAMINATION, "%\n")
  cat("NCBI Tax Directory:", NCBI_TAX_DIR, "\n")
  cat("=====================================\n")
}

# Download file from GitHub if it doesn't exist locally
# GitHub URL should point to raw file content
download_from_github <- function(github_url, local_path, description = "file") {
  if (file.exists(local_path) && file.size(local_path) > 0) {
    return(local_path)
  }
  
  # Create directory if it doesn't exist
  dir.create(dirname(local_path), showWarnings = FALSE, recursive = TRUE)
  
  # Download from GitHub
  tryCatch({
    download.file(github_url, destfile = local_path, mode = "wb", quiet = TRUE)
    if (file.exists(local_path) && file.size(local_path) > 0) {
      message(paste("Downloaded", description, "from GitHub to:", local_path))
      return(local_path)
    } else {
      stop("Downloaded file is empty")
    }
  }, error = function(e) {
    if (file.exists(local_path)) unlink(local_path)
    stop(paste("Failed to download", description, "from GitHub:", e$message,
               "\nURL:", github_url,
               "\nPlease download manually and place at:", local_path))
  })
}

# ============================================================================
# Data Processing Helper Functions
# ============================================================================


# gtdb_214_metadata = fread("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/genomes/bac120_metadata_r214.tsv") %>%
# 	select(accession,checkm_completeness,checkm_contamination,gtdb_genome_representative,gtdb_taxonomy,ncbi_date,ncbi_genbank_assembly_accession,ncbi_organism_name,ncbi_taxid,ncbi_taxonomy) %>% mutate(source="GTDB_214")
#
#
# gtdb_95_metadata = fread("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/genomes/bac120_metadata_r95.tsv") %>%
# 	select(accession,checkm_completeness,checkm_contamination,gtdb_genome_representative,gtdb_taxonomy,ncbi_date,ncbi_genbank_assembly_accession,ncbi_organism_name,ncbi_taxid,ncbi_taxonomy) %>% mutate(source="GTDB_95")

# Load GTDB 207 metadata (required for taxonomy matching)
# Check if we're in a context where BASE_DIR and GENOME_DB_DIR are defined
gtdb_207_metadata <- NULL
if (exists("GENOME_DB_DIR")) {
  gtdb_207_bac_file <- file.path(GENOME_DB_DIR, "bac120_metadata_r207.tsv")
  gtdb_207_ar_file <- file.path(GENOME_DB_DIR, "ar53_metadata_r207.tsv")
  
  if (file.exists(gtdb_207_bac_file) && file.exists(gtdb_207_ar_file)) {
    gtdb_207_metadata_bac <- fread(gtdb_207_bac_file) %>%
      select(accession,checkm_completeness,checkm_contamination,gtdb_genome_representative,gtdb_taxonomy,ncbi_date,ncbi_genbank_assembly_accession,ncbi_organism_name,ncbi_taxid,ncbi_taxonomy, ncbi_genome_category) %>% mutate(source="GTDB_207")
    gtdb_207_metadata_ar <- fread(gtdb_207_ar_file) %>%
      select(accession,checkm_completeness,checkm_contamination,gtdb_genome_representative,gtdb_taxonomy,ncbi_date,ncbi_genbank_assembly_accession,ncbi_organism_name,ncbi_taxid,ncbi_taxonomy, ncbi_genome_category) %>% mutate(source="GTDB_207")
    gtdb_207_metadata <- rbind(gtdb_207_metadata_bac, gtdb_207_metadata_ar) %>% 
      mutate(is_MAG=ifelse(ncbi_genome_category=="derived from metagenome", T, F)) %>% select(-ncbi_genome_category)
  }
}

# Fallback to old hardcoded paths if GENOME_DB_DIR not available or files not found
if (is.null(gtdb_207_metadata)) {
  old_bac_path <- "/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/genomes/bac120_metadata_r207.tsv"
  old_ar_path <- "/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/genomes/ar53_metadata_r207.tsv"
  
  if (file.exists(old_bac_path) && file.exists(old_ar_path)) {
    gtdb_207_metadata_bac <- fread(old_bac_path) %>%
      select(accession,checkm_completeness,checkm_contamination,gtdb_genome_representative,gtdb_taxonomy,ncbi_date,ncbi_genbank_assembly_accession,ncbi_organism_name,ncbi_taxid,ncbi_taxonomy, ncbi_genome_category) %>% mutate(source="GTDB_207")
    gtdb_207_metadata_ar <- fread(old_ar_path) %>%
      select(accession,checkm_completeness,checkm_contamination,gtdb_genome_representative,gtdb_taxonomy,ncbi_date,ncbi_genbank_assembly_accession,ncbi_organism_name,ncbi_taxid,ncbi_taxonomy, ncbi_genome_category) %>% mutate(source="GTDB_207")
    gtdb_207_metadata <- rbind(gtdb_207_metadata_bac, gtdb_207_metadata_ar) %>% 
      mutate(is_MAG=ifelse(ncbi_genome_category=="derived from metagenome", T, F)) %>% select(-ncbi_genome_category)
  }
}

# Note: gtdb_207_metadata may still be NULL if files are not accessible
# Individual scripts should check and load via SSH if needed


# Reads kraken output file, splits taxid vector if needed, and removes duplicated reads
fread_kraken_reads = function(file.path = NULL, db_name = "db_name", sampleID = "HARV_001-M-13-7-20131122"){

	col_vec = c("classified_status","sequenceID","kraken_tax_id","sequence_length","classified_details")

reads_in <- data.table::fread(file.path, sep = "\t",
	col.names = col_vec, nThread = 8)

if (grepl("taxid",reads_in$kraken_tax_id[[1]])) {
	reads_in = reads_in %>%
		# Split named taxid column
		tidyfast::dt_separate(col = kraken_tax_id,
													into = c("taxon","tax_id"), sep="taxid") %>%
		mutate(kraken_tax_id = gsub("\\)|[[:space:]]","", tax_id),
					 taxon = gsub("[[:space:]]\\(","", taxon)) %>% select(-tax_id)
}
# Removes duplicate reads - this should not be necessary. Caused by a problem upstream??
df_return = reads_in %>%
	mutate(db_name = db_name,
				 sampleID = sampleID) %>%
	distinct(db_name,sequenceID,sampleID,.keep_all = T) %>%
	filter(!is.na(sequenceID))  %>%
	mutate(kraken_tax_id = as.character(kraken_tax_id))
return(df_return)
}


# split_taxonomy_ncbi = function(df_in, col = "taxLineage"){
# 	input_data = df_in %>% select(!!col)
split_taxonomy_ncbi = function(input_data){

	root <- str_extract(input_data, "(?<=\\|r_)[^|]+|^r_[^|]+")
	kingdom <- str_extract(input_data, "(?<=\\|k_)[^|]+|^k_[^|]+")
	domain <- str_extract(input_data, "(?<=\\|d_)[^|]+|^d_[^|]+")
	phylum <- str_extract(input_data, "(?<=\\|p_)[^|]+|^p_[^|]+")
	class <- str_extract(input_data, "(?<=\\|c_)[^|]+|^c_[^|]+")
	order <- str_extract(input_data, "(?<=\\|o_)[^|]+|^o_[^|]+")
	family <- str_extract(input_data, "(?<=\\|f_)[^|]+|^f_[^|]+")
	genus <- str_extract(input_data,    "(?<=\\|g_)[^|]+|^g_[^|]+")
	species <- str_extract(input_data, 	"(?<=\\|s_)[^|]+|^s_[^|]+")
	subspecies <- str_extract(input_data, "(?<=\\|s1_)[^|]+|^s1_[^|]+")

	lineage_df = cbind.data.frame(root, domain, kingdom, phylum, class, order, family, genus, species, subspecies)
#
# 	lineage_df = df_in %>% select(!!col) %>% separate(col,into =
# 																											 	c("root","domain",
# 																											 		#"kingdom",
# 																											 		"phylum","class","order","family","genus","species","subspecies"),
# 																											 fill = "right",
# 																											sep="\\|") %>%
# 	 mutate(across(c("root","domain","kingdom","phylum","class","order","family","genus","species","subspecies"), function(x) { gsub("s1_|s_|g_|f_|o_|c_|p_|k_|d_|r_","",x) }))
	return(lineage_df)
#	return(cbind.data.frame(df_in, lineage_df))
}


# from https://alistaire.rbind.io/blog/coalescing-joins/
coalesce_join <- function(x, y,
													by = NULL, suffix = c(".x", ".y"),
													join = dplyr::full_join, ...) {
	joined <- join(x, y, by = by, suffix = suffix, ...)
	# names of desired output
	cols <- union(names(x), names(y))

	to_coalesce <- names(joined)[!names(joined) %in% cols]
	suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
	# remove suffixes and deduplicate
	to_coalesce <- unique(substr(
		to_coalesce,
		1,
		nchar(to_coalesce) - nchar(suffix_used)
	))

	coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
		joined[[paste0(.x, suffix[1])]],
		joined[[paste0(.x, suffix[2])]]
	))
	names(coalesced) <- to_coalesce

	dplyr::bind_cols(joined, coalesced)[cols]
}


clean_ncbi <- function(NCBI_taxon_name_in) {

	# Remove prefixes
	NCBI_taxon_name = gsub("s__|g__|f__|o__|c__|p__","", NCBI_taxon_name_in)

	# Remove percent matches
	NCBI_taxon_name =  gsub("%$","", NCBI_taxon_name)
	NCBI_taxon_name =  gsub(" [[:digit:]][[:digit:]][[:digit:]].[[:digit:]]$","", NCBI_taxon_name)
	NCBI_taxon_name =  gsub(" [[:digit:]][[:digit:]].[[:digit:]][[:digit:]]$","", NCBI_taxon_name)
	NCBI_taxon_name =  gsub(" [[:digit:]][[:digit:]].[[:digit:]]$","", NCBI_taxon_name)
	NCBI_taxon_name =  gsub(" [[:digit:]].[[:digit:]][[:digit:]]$","", NCBI_taxon_name)
	NCBI_taxon_name =  gsub(" [[:digit:]].[[:digit:]]$","", NCBI_taxon_name)

	# Remove parentheses if they follow an assignment
	NCBI_taxon_name <- gsub("\\b([[:alnum:]]+)\\s*\\([^)]*\\)", "\\1", NCBI_taxon_name)

	# Remove parentheses
	NCBI_taxon_name =  gsub("\\(|\\)","", NCBI_taxon_name)

	# Remove parentheses
	NCBI_taxon_name =  gsub("\\[|\\]","", NCBI_taxon_name)

	NCBI_taxon_name
}

read_in_genomes <- function(directory, pattern = ".fa"){

# Return empty data frame if directory doesn't exist
if (!dir.exists(directory)) {
  return(data.frame(filepath = character(0), 
                    filename = character(0), 
                    user_genome = character(0),
                    stringsAsFactors = FALSE))
}

files_downloaded = list.files(directory,
														 full.names = T, #include.dirs = T,
														 recursive = T, pattern = pattern)

# Return empty data frame if no files found
if (length(files_downloaded) == 0) {
  return(data.frame(filepath = character(0), 
                    filename = character(0), 
                    user_genome = character(0),
                    stringsAsFactors = FALSE))
}

names(files_downloaded) = basename(files_downloaded)
files_downloaded <- stack(files_downloaded)
colnames(files_downloaded) = c("filepath","filename")
files_downloaded$filename = as.character(files_downloaded$filename)
files_downloaded$user_genome = gsub(pattern,"",files_downloaded$filename)
return(files_downloaded)
}



tax_id_to_ranked_lineage <- function( ncbi_tax_ids , ncbi_tax_dir =  "dir_path/to/ncbi_taxonomy_files/"){
	# function from https://www.biostars.org/p/317073/
	library(data.table)
	library(tidyverse)
	## get rankedlineage file index
	ranked_lineage_file_index  <- grep("rankedlineage.dmp" , list.files(ncbi_tax_dir , full.names = T))

	## check if the file present
	if(length(ranked_lineage_file_index) != 1){
		stop("rankedlineagess.dmp not found in given tax directory" )
	}

	## get the file and load it
	ranked_lineage_file <- list.files(ncbi_tax_dir , full.names = T)[ranked_lineage_file_index]
	ranked_lineage <- fread(ranked_lineage_file , data.table = F) %>%
		as_tibble() %>%
		select(seq(from = 1, to = 20 , by = 2)) %>%
		`colnames<-`(c("tax_id","tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom"))

	ncbi_tax_ids <- as.numeric(ncbi_tax_ids) %>% as_tibble() %>% `colnames<-`(c("query_tax_id"))
	mapped <- ranked_lineage %>%
		right_join(ncbi_tax_ids , by = c("tax_id" = "query_tax_id"))

	return(mapped)
}


source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/422f4e10d6fbb0ddffe790e3d091bb239f4f3418/scripts/process_classification.R")
add_result_columns <- function(df){


	prok_cols <- function() {
		c(any_of(c("kingdom_GTDB_207", "kingdom_PlusPF","kingdom_PlusPFP8","domain_GTDB_207", "domain_PlusPF","domain_PlusPFP8","domain_soil_microbe_db")))
	}
	viral_cols <- function() {
		c(any_of(c("kingdom_viral", "kingdom_PlusPF", "kingdom_PlusPFP8","domain_PlusPFP8","domain_PlusPF","domain_viral","domain_soil_microbe_db")))
	}

	fungal_cols <- function() {
		c(any_of(c("kingdom_custom_fungi_multi",
							 "kingdom_PlusPF","kingdom_PlusPFP8","kingdom_mycocosm_published","kingdom_mycocosm_all",
							 "domain_soil_microbe_db","domain_PlusPFP8","domain_PlusPF","domain_mycocosm_published")))
	}



	human_cols <- function() {
		c(any_of(c("species_PlusPF","species_PlusPFP8")))
	}


	all_cols <- function() {
		c(any_of(c("classified_status_GTDB_207","classified_status_PlusPF",
							 "classified_status_mycocosm_all", "classified_status_mycocosm_published",
							 "classified_status_viral","classified_status_custom_fungi_multi",
							 "classified_status_PlusPFP8","classified_status_soil_microbe_db"
							 )))
	}



	plant_cols <- function() {
		c(any_of(c("kingdom_PlusPFP8","domain_PlusPFP8")))
	}
	df = df %>% mutate(prokaryote =	ifelse(if_any(prok_cols(), ~str_detect(., "Bacteria|Archaea")), 1, 0))
	df = df %>% mutate(viral =	ifelse(if_any(viral_cols(), ~str_detect(., "Viruses")), 1, 0))
	df = df %>% mutate(fungal =	ifelse(if_any(fungal_cols(), ~str_detect(., "Fungi")), 1, 0))
	df = df %>% mutate(human =	ifelse(if_any(human_cols(), ~str_detect(., "sapiens")), 1, 0))
	df = df %>% mutate(plant =	ifelse(if_any(plant_cols(), ~str_detect(., "Plantae|plant")), 1, 0))
#	df = df %>% mutate(plant =	ifelse(grepl("Plantae", kingdom_PlusPFP8), 1, 0))

	df$prokaryote[is.na(df$prokaryote)] <- 0
	df$viral[is.na(df$viral)] <- 0
	df$fungal[is.na(df$fungal)] <- 0
	df$human[is.na(df$human)] <- 0
	df$plant[is.na(df$plant)] <- 0

	df$classified_as_bacterial_and_fungal = ifelse(df$prokaryote == 1 & df$fungal == 1, 1, 0)
	df$classified_as_bacterial_and_viral = ifelse(df$prokaryote == 1 & df$viral == 1, 1, 0)
	df$classified_as_human_and_fungal = ifelse(df$human == 1 & df$fungal == 1, 1, 0)
	df$classified_as_plant_and_bacterial = ifelse(df$prokaryote == 1 & df$plant == 1, 1, 0)

	df = df %>% mutate(classified =	ifelse(if_any(all_cols(), ~str_detect(., "C")), 1, 0))
	df$classified[is.na(df$classified)] <- 0
	df$classified_as_bacterial_and_fungal[is.na(df$classified_as_bacterial_and_fungal)] <- 0
	df$classified_as_bacterial_and_viral[is.na(df$classified_as_bacterial_and_viral)] <- 0


	df$classified_only_myco_all =
		ifelse(df$classified_status_mycocosm_all == "C" &
					 	(df$classified_status_mycocosm_published == "U" | is.na(df$classified_status_mycocosm_published)), 1, 0)
	df$classified_only_myco_all[is.na(df$classified_only_myco_all)] <- 0

 return(df)
}



# from https://github.com/larssnip/microclass/blob/3f3ac0ad10650957168e1d9d903e8eea138aae36/R/krkFun.R#L273
kmer_report <- function(kmer.tbl, taxonomy.tbl){
	kmer.tbl %>%
		filter(tax_id == 0) -> unclas.tbl
	kmer.tbl %>%
		filter(tax_id != 0) -> kmr.tbl

	utxm <- numeric(0)
	branch.mat <- attr(taxonomy.tbl, "branch.mat")
	for(i in nrow(kmr.tbl):1){
		if(!(kmr.tbl$tax_id[i] %in% utxm)){
			M <- which(branch.mat == kmr.tbl$tax_id[i], arr.ind = T)
			txm <- unique(as.integer(branch.mat[M[,1], 1:M[1,2]]))
			utxm <- c(utxm, txm)
		}
	}
	utxm <- unique(utxm[!is.na(utxm)])
	taxonomy.tbl %>%
		filter(tax_id %in% utxm) -> rep.tbl
	idx <- match(rep.tbl$tax_id, kmr.tbl$tax_id)
	idd <- which(!is.na(idx))
	rep.tbl$tax_count[idd] <- kmr.tbl$tax_count[idx[idd]]
	bss <- c(0, which(diff(rep.tbl$rank_int) <= 0), nrow(rep.tbl))

	t_count <- rep.tbl$tax_count
	for(i in length(bss):2){
		M <- which(branch.mat == rep.tbl$tax_id[bss[i]], arr.ind = T)
		bm <- branch.mat[M[1,1], 1:M[1,2]]
		bm <- bm[!is.na(bm)]
		rr <- match(bm, rep.tbl$tax_id)
		cc <- rev(cumsum(rev(t_count[rr])))
		rep.tbl$clade_count[rr] <- rep.tbl$clade_count[rr] + cc
		t_count[rr] <- 0
	}

	N.kmers <- rep.tbl$clade_count[1]
	if(nrow(unclas.tbl) > 0){
		N.kmers <- N.kmers + unclas.tbl$tax_count
		tibble(percent = 0,
					 clade_count = unclas.tbl$tax_count,
					 tax_count = unclas.tbl$tax_count,
					 rank = "U",
					 tax_id = unclas.tbl$tax_id,
					 name = "unclassified",
					 rank_int = 0) %>%
			bind_rows(rep.tbl) -> rep.tbl
	}
	rep.tbl %>%
		mutate(percent = clade_count / N.kmers) %>%
		select(-rank_int) %>%
		mutate(read_id = kmr.tbl$read_id[1]) -> rep.tbl
	return(rep.tbl)
}


kraken2_taxonomy <- function(inspect.file){
	rank.order <- c("R", "D", "K", "P", "C", "O", "F", "G", "S")
	read_kraken2_report(inspect.file) %>%
		filter(rank != "U") -> insp.tbl
	urank <- sort(unique(insp.tbl$rank))
	levs <- unlist(sapply(rank.order, function(r){urank[str_detect(urank, r)]}))
	insp.tbl %>%
		mutate(rank = factor(rank, levels = levs)) %>%
		mutate(rank_int = as.integer(rank)) %>%
		mutate(rank = as.character(rank)) %>%
		mutate(percent = 0, tax_count = 0, clade_count = 0) -> insp.tbl
	bss <- c(0, which(diff(insp.tbl$rank_int) <= 0), nrow(insp.tbl))
	branch.mat <- matrix(NA, nrow = length(bss)-1, ncol = length(levs))
	colnames(branch.mat) <- levs
	rownames(branch.mat) <- insp.tbl$tax_id[bss[-1]]
	tax <- rep(NA, length(levs))
	for(i in 2:length(bss)){
		rr <- (bss[i-1] + 1):bss[i]
		iii <- insp.tbl$rank_int[rr]
		if(min(iii) > 1){
			tax <- c(tax[1:(min(iii) - 1)], rep(NA, length(levs) - min(iii) + 1))
		}
		tax[iii] <- insp.tbl$tax_id[rr]
		branch.mat[i-1,] <- tax
	}
	attr(insp.tbl, "branch.mat") <- branch.mat
	return(insp.tbl)
}


# create sample information data.frame from NEON sample names
parseNEONsampleIDs <- function(sampleID){
	df <- data.frame(siteID = substr(sampleID, 1, 4), sampleID = sampleID, stringsAsFactors = F) %>%
		mutate(sample = sapply(strsplit(sampleID, "-GEN|-gen"),  "[[" , 1)) %>%
		mutate(geneticSampleID = sapply(strsplit(sampleID, "-DNA"),  "[[" , 1)) %>%
		mutate(sampleID = sapply(strsplit(sampleID, "-gen.fastq"),  "[[" , 1)) %>%
		mutate(dates = sapply(strsplit(sample, "-"), function(x) x[grep("[2]\\d\\d\\d\\d\\d\\d\\d", x)])) %>%
		mutate(dates = ifelse(dates == "21040514", "20140514", dates)) %>%
		mutate(asDate = as.Date(as.character(dates), "%Y%m%d")) %>%
		mutate(dateID = substr(as.character(dates), 1, 6)) %>%
		mutate(plotID = substr(sample, 1, 8)) %>%
		mutate(site_date = paste0(siteID, "-", dateID)) %>%
		mutate(horizon = ifelse(grepl("-M-", sample), "M", "O")) %>%
		mutate(without_horizon = gsub("-[M|O]-", "-", sample)) %>%
		mutate(plot_date = paste0(plotID, "-", dateID)) %>%
		as.data.frame()
	rownames(df) <- make.unique(sampleID)
	return(df)
}

