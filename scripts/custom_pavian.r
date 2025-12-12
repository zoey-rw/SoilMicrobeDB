#custom pavian functions to work with struo kraken output

myfile="/projectnb/dietzelab/zrwerbin/soil_genome_db_output/kraken_read_assignments/pluspf/HARV_001-M-13-7-20131122.kreport"
myfile="/projectnb/dietzelab/zrwerbin/soil_genome_db_output/kraken_read_assignments/pluspf/HARV_002-O-16-1.5-20170801.kreport"
fread_report = function(myfile, sampleID = NULL, db_name = NULL, keep_taxRanks = c("D", "K",
																																									 "P", "C", "O", "F", "G", "S")){
	keep_taxRanks = c("D", "K",
										"P", "C", "O", "F", "G", "S","S1")
	report = data.table::fread(myfile, sep = "\t")
	if(ncol(report)<8) {
		col_vec <-  c("percentage", "cladeReads", "taxonReads",
									 	"taxRank", "taxID", "name")
	} else {
		col_vec <- c("percentage", "cladeReads", "taxonReads", "nKmers","n_unique_kmers",
									 	"taxRank", "taxID", "name")
}

colnames(report) <- col_vec
#report$percentage <- round(report$cladeReads/sum(report$taxonReads),6) * 100


if ('n_unique_kmers'  %in% colnames(report)){
	report$kmerpercentage <- round(report$n_unique_kmers/sum(report$n_unique_kmers,na.rm=T),6) * 100
}

report$depth <- nchar(gsub("\\S.*", "", report$name))/2
report$name <- gsub("^ *", "", report$name)
report$name <- paste(tolower(report$taxRank), report$name,
										 sep = "_")

# This is the function that takes a minute

report = report %>%
	mutate(db_name = !!db_name, sampleID = !!sampleID,
				 taxID = as.double(taxID))
return(report)
}


read_report3 <- function(myfile,collapse=TRUE,keep_taxRanks=c("D","K","P","C","O","F","G","S"),min.depth=0,filter_taxon=NULL,
												 has_header=NULL,add_taxRank_columns=FALSE) {

	first.line <- readLines(myfile,n=1)
	isASCII <-  function(txt) all(charToRaw(txt) <= as.raw(127))
	if (!isASCII(first.line)) {
		dmessage(myfile," is no valid report - not all characters are ASCII")
		return(NULL)
	}
	if (is.null(has_header)) {
		has_header <- grepl("^[a-zA-Z]",first.line)
	}

	if (has_header) {
		report <- utils::read.table(myfile,sep="\t",header = T,
																quote = "",stringsAsFactors=FALSE, comment.char="#")
		#colnames(report) <- c("percentage","cladeReads","taxonReads","taxRank","taxID","n_unique_kmers","n_kmers","perc_uniq_kmers","name")

		## harmonize column names. TODO: Harmonize them in the scripts!
		colnames(report)[colnames(report)=="clade_perc"] <- "percentage"
		colnames(report)[colnames(report)=="perc"] <- "percentage"

		colnames(report)[colnames(report)=="n_reads_clade"] <- "cladeReads"
		colnames(report)[colnames(report)=="n.clade"] <- "cladeReads"

		colnames(report)[colnames(report)=="n_reads_taxo"] <- "taxonReads"
		colnames(report)[colnames(report)=="n.stay"] <- "taxonReads"

		colnames(report)[colnames(report)=="rank"] <- "taxRank"
		colnames(report)[colnames(report)=="tax_rank"] <- "taxRank"

		colnames(report)[colnames(report)=="taxonid"] <- "taxID"
		colnames(report)[colnames(report)=="tax"] <- "taxID"

	} else {
		report <- utils::read.table(myfile,sep="\t",header = F,
																col.names = c("percentage","cladeReads","taxonReads","taxRank","taxID","name"),
																quote = "",stringsAsFactors=FALSE, comment.char="#")
	}

	report$depth <- nchar(gsub("\\S.*","",report$name))/2
	report$name <- gsub("^ *","",report$name)
	report$name <- paste(tolower(report$taxRank),report$name,sep="_")


	## Only stop at certain taxRanks
	## filter taxon and further up the tree if 'filter_taxon' is defined
	kraken.tree <- pavian:::build_kraken_tree(report)
	#report <- pavian:::collapse.taxRanks(kraken.tree,keep_taxRanks=keep_taxRanks,filter_taxon=filter_taxon)

	## Add a metaphlan-style taxon string
	# if (add_taxRank_columns) {
	# 	report[,keep_taxRanks] <- NA
	# }
	report$taxLineage = report$name
	rows_to_consider <- rep(FALSE,nrow(report))

	report$percentage <- round(report$cladeReads/sum(report$taxonReads),6) * 100

	for (column in c("taxonReads", "cladeReads"))
		if (all(floor(report[[column]]) == report[[column]]))
			report[[column]] <- as.integer(report[[column]])

	if ('n_unique_kmers'  %in% colnames(report))
		report$kmerpercentage <- round(report$n_unique_kmers/sum(report$n_unique_kmers,na.rm=T),6) * 100
	#report$taxRankperc <- 100/taxRank(report$cladeReads)

	rownames(report) <- NULL

	report
}

summarize_report_custom <- function(my_report) {
	protist_taxids <- c("-_Diplomonadida"=5738,
											"-_Amoebozoa"=554915,
											"-_Alveolata"=33630)
	rownames(my_report)[rownames(my_report) == "r_root"] <- "-_root"
	my_report <- my_report[!duplicated(my_report$name),]
	row.names(my_report) <- my_report[["name"]]
	unidentified_reads <- my_report["u_unclassified","cladeReads"]
	identified_reads <- my_report["r_root","cladeReads"]
	artificial_reads <- pavian:::zero_if_na1(my_report["s_synthetic construct",
																										 "cladeReads"])
	human_reads <- pavian:::zero_if_na1(my_report["s_Homo sapiens", "cladeReads"])
	chordate_reads <- pavian:::zero_if_na1(my_report["p_Chordata", "cladeReads"])
	root_reads <- pavian:::zero_if_na1(my_report["r_root", "taxonReads"])
	out <- data.frame(number_of_raw_reads = unidentified_reads + identified_reads,
										classified_reads = identified_reads, chordate_reads = chordate_reads,
										artificial_reads = artificial_reads, unclassified_reads = unidentified_reads,
										microbial_reads = identified_reads - chordate_reads -
											artificial_reads - root_reads, bacterial_reads = pavian:::zero_if_na1(my_report["d_Bacteria",
																																																			"cladeReads"]) + pavian:::zero_if_na1(my_report["k_Bacteria",
																																																																											"cladeReads"]), viral_reads = pavian:::zero_if_na1(my_report["d_Viruses",
																																																																																																									 "cladeReads"]) + pavian:::zero_if_na1(my_report["k_Viruses",
																																																																																																									 																								"cladeReads"]), fungal_reads = pavian:::zero_if_na1(my_report["k_Fungi",
																																																																																																									 																																																							"cladeReads"]), protozoan_reads = sum(pavian:::zero_if_na1(my_report[names(protist_taxids),
																																																																																																									 																																																																																									 "cladeReads"])))
}



fread_report2 <- function (myfile, db_name=NULL, sampleID = NULL, collapse = FALSE, keep_taxRanks = c("D", "K",
																										 "P", "C", "O", "F", "G", "S","S1"), min.depth = 0, filter_taxon = NULL,
					has_header = NULL, add_taxRank_columns = FALSE)
{
	first.line <- readLines(myfile, n = 1)
	isASCII <- function(txt) all(charToRaw(txt) <= as.raw(127))
	if (!isASCII(first.line)) {
		dmessage(myfile, " is no valid report - not all characters are ASCII")
		return(NULL)
	}
	if (is.null(has_header)) {
		has_header <- grepl("^[a-zA-Z]", first.line)
	}
	if (has_header) {
		report <- utils::read.table(myfile, sep = "\t", header = T,
																quote = "", stringsAsFactors = FALSE, comment.char = "#")
		colnames(report)[colnames(report) == "clade_perc"] <- "percentage"
		colnames(report)[colnames(report) == "perc"] <- "percentage"
		colnames(report)[colnames(report) == "n_reads_clade"] <- "cladeReads"
		colnames(report)[colnames(report) == "n.clade"] <- "cladeReads"
		colnames(report)[colnames(report) == "n_reads_taxo"] <- "taxonReads"
		colnames(report)[colnames(report) == "n.stay"] <- "taxonReads"
		colnames(report)[colnames(report) == "rank"] <- "taxRank"
		colnames(report)[colnames(report) == "tax_rank"] <- "taxRank"
		colnames(report)[colnames(report) == "taxonid"] <- "taxID"
		colnames(report)[colnames(report) == "tax"] <- "taxID"
	}
	else {
		report <- utils::read.table(myfile, sep = "\t", header = F,
																# col.names = c("percentage", "cladeReads", "taxonReads",
																# 							"taxRank", "taxID", "name"),
																quote = "", stringsAsFactors = FALSE,
																comment.char = "#")
		if(ncol(report)<8) {
			col_vec <-  c("percentage", "cladeReads", "taxonReads",
										"taxRank", "taxID", "name")
		} else {
			col_vec <- c("percentage", "cladeReads", "taxonReads", "nKmers","n_unique_kmers",
									 "taxRank", "taxID", "name")
		}

		colnames(report) <- col_vec
	}
	report$depth <- nchar(gsub("\\S.*", "", report$name))/2
	report$name <- gsub("^ *", "", report$name)
	report$name <- paste(tolower(report$taxRank), report$name,
											 sep = "_")
	kraken.tree <- pavian:::build_kraken_tree(report)
	report <- pavian:::collapse.taxRanks(kraken.tree, keep_taxRanks = keep_taxRanks,
															filter_taxon = filter_taxon)
	if (add_taxRank_columns) {
		report[, keep_taxRanks] <- NA
	}
	report$taxLineage = report$name
	rows_to_consider <- rep(FALSE, nrow(report))
	for (i in seq_len(nrow(report))) {
		if (i > 1 && report[i, "depth"] > min.depth) {
			idx <- report$depth < report[i, "depth"] & rows_to_consider
			if (!any(idx)) {
				(next)()
			}
			current.taxRank <- report[i, "taxRank"]
			my_row <- max(which(idx))
			report[i, "taxLineage"] <- paste(report[my_row, "taxLineage"],
																			 report[i, "taxLineage"], sep = "|")
			if (add_taxRank_columns) {
				if (report[my_row, "taxRank"] %in% keep_taxRanks) {
					taxRanks.cp <- keep_taxRanks[seq(from = 1,
																					 to = which(keep_taxRanks == report[my_row,
																					 																	 "taxRank"]))]
					report[i, taxRanks.cp] <- report[my_row, taxRanks.cp]
				}
				report[i, report[i, "taxRank"]] <- report[i,
																									"name"]
			}
		}
		rows_to_consider[i] <- TRUE
	}
	report <- report[report$depth >= min.depth, ]
	report$percentage <- round(report$cladeReads/sum(report$taxonReads),
														 6) * 100
	for (column in c("taxonReads", "cladeReads")) if (all(floor(report[[column]]) ==
																												report[[column]]))
		report[[column]] <- as.integer(report[[column]])
	if ("n_unique_kmers" %in% colnames(report))
		report$kmerpercentage <- round(report$n_unique_kmers/sum(report$n_unique_kmers,
																														 na.rm = T), 6) * 100
	rownames(report) <- NULL

	report = report %>%
		mutate(db_name = !!db_name, sampleID = !!sampleID,
					 taxID = as.double(taxID))
	report
}

#krk.tbl = reads_pluspf[1:2000,] %>% rename(read_id = sequenceID, tax_id = kraken_tax_id, tax_count = classified_details, length = sequence_length)
# From: https://rdrr.io/github/larssnip/microclass/src/R/krkFun.R#sym-krk2kmer
krk2kmer <- function(krk.tbl){
	krk.tbl %>%
		select(db_name, sampleID, sequenceID, classified_details) %>%
		mutate(tax_count = str_remove(classified_details, "\\|:\\| ")) %>%
		mutate(tax_id = str_remove_all(tax_count, ":[0-9]+")) %>%
		mutate(count = str_remove_all(tax_count, "[0-9]+:")) %>%
		select(-tax_count) %>%
		mutate(tax_id = str_split(tax_id, pattern = " ")) %>%
		mutate(count = str_split(count, pattern = " ")) %>%
		unnest(c(tax_id, count)) %>%
		mutate(tax_id = as.double(tax_id), count = as.integer(count)) %>%
		group_by(db_name, sampleID, sequenceID, tax_id) %>%
		summarise(tax_count = sum(count)) %>%
		ungroup() -> kmer.tbl
	return(kmer.tbl)
}



