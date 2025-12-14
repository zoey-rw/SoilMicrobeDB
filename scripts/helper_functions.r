pacman::p_load(data.table, tidyverse, tidyfast, phyloseq, CHNOSZ)


# gtdb_214_metadata = fread("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/genomes/bac120_metadata_r214.tsv") %>%
# 	select(accession,checkm_completeness,checkm_contamination,gtdb_genome_representative,gtdb_taxonomy,ncbi_date,ncbi_genbank_assembly_accession,ncbi_organism_name,ncbi_taxid,ncbi_taxonomy) %>% mutate(source="GTDB_214")
#
#
# gtdb_95_metadata = fread("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/genomes/bac120_metadata_r95.tsv") %>%
# 	select(accession,checkm_completeness,checkm_contamination,gtdb_genome_representative,gtdb_taxonomy,ncbi_date,ncbi_genbank_assembly_accession,ncbi_organism_name,ncbi_taxid,ncbi_taxonomy) %>% mutate(source="GTDB_95")

gtdb_207_metadata_bac = fread("data/genome_database/bac120_metadata_r207.tsv") %>%
	select(accession,checkm_completeness,checkm_contamination,gtdb_genome_representative,gtdb_taxonomy,ncbi_date,ncbi_genbank_assembly_accession,ncbi_organism_name,ncbi_taxid,ncbi_taxonomy, ncbi_genome_category) %>% mutate(source="GTDB_207")
gtdb_207_metadata_ar = fread("data/genome_database/ar53_metadata_r207.tsv") %>%
	select(accession,checkm_completeness,checkm_contamination,gtdb_genome_representative,gtdb_taxonomy,ncbi_date,ncbi_genbank_assembly_accession,ncbi_organism_name,ncbi_taxid,ncbi_taxonomy, ncbi_genome_category) %>% mutate(source="GTDB_207")
gtdb_207_metadata = rbind(gtdb_207_metadata_bac, gtdb_207_metadata_ar) %>% 
    mutate(is_MAG=ifelse(ncbi_genome_category=="derived from metagenome", T, F)) %>% select(-ncbi_genome_category)


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

files_downloaded = list.files(directory,
														 full.names = T, #include.dirs = T,
														 recursive = T, pattern = pattern)

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


# Functions from helper_functions2.r (correlation and biome assignment functions)

# Define functions to compute correlations (this could be rewritten as one function with a second argument for variable)
cor_fun_pH <- function(df) cor.test(df$percentage, 
																		df$soilInCaClpH, 
																		method = "spearman", exact = FALSE) %>% 
	tidy()
cor_fun_nitr <- function(df) cor.test(df$percentage, 
																			df$nitrogenPercent, 
																			method = "spearman", exact = FALSE) %>% 
	tidy()
cor_fun_carbon <- function(df) cor.test(df$percentage, 
																				df$organicd13C, 
																				method = "spearman", exact = FALSE) %>% 
	tidy()
cor_fun_temperature <- function(df) cor.test(df$percentage, 
																						 df$soilTemp, 
																						 method = "spearman", exact = FALSE) %>% 
	tidy()
cor_fun_moisture <- function(df) cor.test(df$percentage, 
																					df$soilMoisture, 
																					method = "spearman", exact = FALSE) %>% 
	tidy()


assign_biome_presence = function(test_df) {
	
	# get counts of samples per biome 
	sample_counts = test_df %>% group_by(nlcdClass) %>% tally(name = "biome_count")
	
	# get counts of non-zero abundances per biome
	presence_counts = test_df %>% mutate(present = ifelse(percentage > 0, 1, 0)) %>% 
		group_by(nlcdClass, present) %>% tally(name = "presence_count")
	
	
	presence_counts_expanded = presence_counts %>% ungroup %>%  expand(nlcdClass, present)
	
	# set NA values to zero
	presence_counts_expanded = merge(presence_counts, presence_counts_expanded, all=T) %>% 
		replace_na(list(presence_count = 0))
	
	df_merged = merge(presence_counts_expanded, sample_counts, all=T) %>% 
		filter(present == 1) %>% 
		mutate(prevalence = presence_count/biome_count,
					 present_in_biome = ifelse(prevalence > .02, 1, 0)) %>% 
		select(nlcdClass,biome_count, prevalence, present_in_biome)
	
	df_out = df_merged  %>% 
		select(nlcdClass,biome_count, prevalence, present_in_biome)
	return(df_out)
}

# Load NEON soil chemistry data
# Checks local path first, then falls back to project path
load_soilChem <- function(neon_soil_file = NULL) {
    if(is.null(neon_soil_file)) {
        local_path <- "data/environmental_data/neon_soil_data_2023.rds"
        project_path <- "/projectnb/dietzelab/zrwerbin/N-cycle/neon_soil_data_2023.rds"
        
        if(file.exists(local_path)) {
            neon_soil_file <- local_path
        } else if(file.exists(project_path)) {
            neon_soil_file <- project_path
        } else {
            stop("neon_soil_data_2023.rds not found in:\n  ", local_path, "\n  ", project_path)
        }
    }
    return(readRDS(neon_soil_file))
}

# Extract genomic sample mapping (sampleID to genomicsSampleID) from soilChem
# Can take soilChem as input or load it if not provided
load_genSampleExample <- function(soilChem = NULL, neon_soil_file = NULL) {
    if(is.null(soilChem)) {
        soilChem <- load_soilChem(neon_soil_file)
    }
    if(!"sls_metagenomicsPooling" %in% names(soilChem)) {
        stop("sls_metagenomicsPooling not found in soilChem")
    }
    genomicSamples <- soilChem$sls_metagenomicsPooling %>%
        tidyr::separate(genomicsPooledIDList, into = c("first", "second", "third"), sep = "\\|", fill = "right") %>%
        dplyr::select(genomicsSampleID, first, second, third)
    genSampleExample <- genomicSamples %>%
        tidyr::pivot_longer(cols = c("first", "second", "third"), values_to = "sampleID") %>%
        dplyr::select(sampleID, genomicsSampleID) %>%
        drop_na()
    return(genSampleExample)
}

# Load and process soil cores with transformations
# Applies transformations: isForest classification, biome recoding, and compositeSampleID addition
# Checks local path first, then falls back to project path
load_soilCores <- function(neon_soil_file = NULL) {
    if(is.null(neon_soil_file)) {
        local_path <- "data/environmental_data/neon_soil_data_2023.rds"
        project_path <- "/projectnb/dietzelab/zrwerbin/N-cycle/neon_soil_data_2023.rds"
        
        if(file.exists(local_path)) {
            neon_soil_file <- local_path
        } else if(file.exists(project_path)) {
            neon_soil_file <- project_path
        } else {
            stop("neon_soil_data_2023.rds not found in:\n  ", local_path, "\n  ", project_path)
        }
    }
    soilData <- readRDS(neon_soil_file)
    soilCores <- soilData$sls_soilCoreCollection
    
    # Apply transformations from source.R
    soilCores <- soilCores %>% 
        mutate(isForest = ifelse(grepl("Forest", nlcdClass), 
                                 "forest habitat", "non forest habitat")) %>% 
        mutate(biome = recode(nlcdClass, "mixedForest" = "Forest",
                             "evergreenForest" = "Forest",
                             "deciduousForest" = "Forest",
                             "emergentHerbaceousWetlands" = "Wetlands",
                             "woodyWetlands" = "Wetlands",
                             "dwarfScrub" = "Shrubland",
                             "shrubScrub" = "Shrubland",
                             "sedgeHerbaceous" = "Herbaceous",
                             "grasslandHerbaceous" = "Herbaceous",
                             "pastureHay" = "Agricultural",
                             "cultivatedCrops" = "Agricultural"))
    
    # Add compositeSampleID from genomicSamples
    if("sls_metagenomicsPooling" %in% names(soilData)) {
        genSampleExample <- load_genSampleExample(soilChem = soilData)
        soilCores$compositeSampleID <- genSampleExample[match(soilCores$sampleID, genSampleExample$sampleID),]$genomicsSampleID
    }
    
    return(soilCores)
}

# From testing 
# (commented out example code - kept for reference)
# Combine abundance and metadata into one dataframe 
# merged_df = merge(species_abundances, 
# 									soil_metadata, by.x = "sampleID", by.y = "genomicsSampleID")
# 
# merged_df_nest = merged_df %>% 
# 	select(-c(sampleID)) %>% 
# 	group_by(db_name, taxon) %>% nest
# 
# test_df = merged_df_nest[,3][[1]][[2]]
# 
# 
# # Filter by correlation strength
# data_nest_biome <- merged_df_nest[1:100,] %>% 
# 	mutate(biome_presence = map(data, assign_biome_presence)) 
# 
# df_biome = data_nest_biome %>% 
# 	select(-data) %>% unnest(cols = c(biome_presence)) %>% ungroup 
# 	
# df_biome_wide = df_biome %>% 
# 	select(-c(biome_count, prevalence)) %>% 
# 	pivot_wider(names_from = nlcdClass, values_from = present_in_biome)
