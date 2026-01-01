#!/usr/bin/env Rscript
# 19: Misclassification analysis and impact assessment
#   - K-mer investigation for missing organisms (from script 20)
#   - Present but misclassified investigation (from script 21)
#   - Analysis across kingdoms, databases, proportions, filtering (from script 22)
#   - Impact assessment and categorization of fixable vs unresolvable (from script 23)
#
# Usage: Rscript scripts/summarize_outputs/19_misclassification_analysis.r [sampleID]
#
# Input:  
#   - {sampleID}_blast_vs_databases_comparison.csv (from script 18)
#   - taxonomic_assignment_comparison_{sampleID}.csv (from script 15)
#   - Genome table files
#   - Bracken abundance estimates
# Output: 
#   - {sampleID}_kmer_investigation_summary.txt
#   - {sampleID}_present_but_misclassified_analysis.txt
#   - {sampleID}_present_but_misclassified_detailed.csv
#   - {sampleID}_misclassification_analysis.csv (table format)
#   - {sampleID}_overall_classification_stats.csv
#   - {sampleID}_blast_verification_stats.csv
#   - {sampleID}_conflict_type_stats.csv
#   - {sampleID}_kingdom_blast_category_stats.csv
#   - {sampleID}_misclassification_impact_analysis.txt
#   - {sampleID}_misclassification_categorization.csv

library(tidyverse)
library(data.table)

# Configuration
sampleID <- if(length(commandArgs(trailingOnly = TRUE)) > 0) {
    commandArgs(trailingOnly = TRUE)[1]
} else {
    "ORNL_046-O-20170621-COMP"
}

output_dir <- "data/classification/analysis_files"
scores_dir <- "data/classification/02_bracken_output"
genome_table_dir <- "data/genome_database"

# Confidence filter thresholds
conf_thresholds <- list(
    consistency_min = 0.9,
    multiplicity_max = 2,
    entropy_max = 0.1
)

# BLAST quality thresholds
blast_quality_thresholds <- list(
    high_identity = 95,
    high_bitscore = 100,
    medium_identity = 85,
    medium_bitscore = 50
)

# Load data

# Load comparison data
comparison_file <- file.path(output_dir, paste0("taxonomic_assignment_comparison_", sampleID, ".csv"))
if(!file.exists(comparison_file)) {
    stop("Comparison file not found. Run script 15 first.")
}
comparison <- read_csv(comparison_file, show_col_types = FALSE) %>%
    mutate(read_id = as.character(read_id))

# Load BLAST comparison data
blast_comparison_file <- file.path(output_dir, paste0(sampleID, "_blast_vs_databases_comparison.csv"))
if(!file.exists(blast_comparison_file)) {
    stop("BLAST comparison file not found. Run script 18 first.")
}
blast_data <- read_csv(blast_comparison_file, show_col_types = FALSE) %>%
    filter(read_id != "read_id") %>%
    mutate(across(any_of(c("identity", "bitscore", "align_length", "evalue")), as.numeric))

# Load scores files and add to comparison
scores_files <- list(
    gtdb = file.path(scores_dir, paste0(sampleID, "_gtdb_207_unfiltered_scores.output")),
    smdb = file.path(scores_dir, paste0(sampleID, "_soil_microbe_db_scores.output")),
    pluspf = file.path(scores_dir, paste0(sampleID, "_pluspf_scores.output"))
)

for(db in names(scores_files)) {
    if(file.exists(scores_files[[db]])) {
        scores <- fread(scores_files[[db]], 
                       select = c("read_id", "n_kmers", "consistency", "multiplicity", "entropy"), 
                       showProgress = FALSE) %>%
        mutate(read_id = as.character(read_id))
        comparison <- comparison %>% 
            left_join(scores, by = "read_id", suffix = c("", paste0("_", db, "_scores")))
    }
}

# Pre-process numeric columns and add confidence flags
comparison <- comparison %>%
    mutate(across(ends_with(c("consistency", "multiplicity", "entropy", "n_kmers")), as.numeric)) %>%
    mutate(
        gtdb_is_confident = ifelse(!is.na(gtdb_consistency),
                                   gtdb_consistency >= conf_thresholds$consistency_min & 
                                   gtdb_multiplicity <= conf_thresholds$multiplicity_max & 
                                   gtdb_entropy <= conf_thresholds$entropy_max,
                                   FALSE),
        smdb_is_confident = ifelse(!is.na(smdb_consistency),
                                  smdb_consistency >= conf_thresholds$consistency_min & 
                                  smdb_multiplicity <= conf_thresholds$multiplicity_max & 
                                  smdb_entropy <= conf_thresholds$entropy_max,
                                  FALSE),
        pluspf_is_confident = ifelse(!is.na(pluspf_consistency),
                                    pluspf_consistency >= conf_thresholds$consistency_min & 
                                    pluspf_multiplicity <= conf_thresholds$multiplicity_max & 
                                    pluspf_entropy <= conf_thresholds$entropy_max,
                                    FALSE)
    )

# Load Bracken abundance data
bracken_file <- file.path("data/classification/taxonomic_rank_summaries", 
                          "soil_microbe_db_species_merged_lineage.csv")
bracken_data <- NULL
if(file.exists(bracken_file)) {
    all_bracken <- fread(bracken_file, showProgress = FALSE)
    sample_patterns <- c(
        paste0(sampleID, "-COMP_soil_microbe_db_filtered"),
        paste0(sampleID, "-COMP_soil_microbe_db"),
        paste0(sampleID, "_soil_microbe_db"),
        paste0(sampleID, "-COMP"),
        sampleID
    )
    for(pattern in sample_patterns) {
        matches <- all_bracken %>% filter(grepl(pattern, sample_id, fixed = TRUE))
        if(nrow(matches) > 0) {
            bracken_data <- matches %>%
                select(taxonomy_id, name, new_est_reads, fraction_total_reads)
            break
        }
    }
}

# Load genome tables
genome_table_files <- list(
    smdb = file.path(genome_table_dir, "soil_microbe_db_genome_table.csv"),
    gtdb = file.path(genome_table_dir, "gtdb_genome_table.csv"),
    pluspf = file.path(genome_table_dir, "pluspf_genome_table.csv")
)

genome_tables <- map(genome_table_files, ~{
    if(file.exists(.x)) read_csv(.x, show_col_types = FALSE) else NULL
}) %>% compact()

# Function to extract organism name from BLAST title (vectorized)
# Regex can handle NCBI titles with prefixes like "PREDICTED:", accession numbers, gi numbers, etc.
extract_organism_name <- function(title) {
    # Remove NCBI boilerplate patterns:
    # 1. gi|numbers|emb|accession| or gi|numbers|gb|accession| or gi|numbers|ref|accession|
    # 2. PREDICTED: prefix
    # 3. Standalone accessions like CP01234.1, NR_1234.1 at start
    clean_title <- str_remove(title, "^gi\\|[0-9]+\\|[a-z]+\\|[A-Z0-9_.]+\\|\\s*")
    clean_title <- str_remove(clean_title, "^(PREDICTED:|[A-Z]{1,2}_[0-9]+\\.[0-9]+ |[A-Z0-9]+ )")
    
    # Extract Binomial (Genus species) - look anywhere in title for Capital lowercase space lowercase pattern
    # This handles cases where organism name comes after accession numbers
    binomial <- str_extract(clean_title, "([A-Z][a-z]+ [a-z]+)")
    
    # If no binomial found, try after "MAG:" or "uncultured" patterns (vectorized)
    mag_match <- str_extract(clean_title, "(MAG:|uncultured|genome assembly, chromosome|strain|isolate)\\s+([A-Z][a-z]+ [a-z]+)")
    mag_binomial <- if_else(!is.na(mag_match), str_extract(mag_match, "([A-Z][a-z]+ [a-z]+)"), NA_character_)
    binomial <- if_else(!is.na(binomial), binomial, mag_binomial)
    
    # If still no binomial, try genus anywhere in title (vectorized)
    genus_start <- str_extract(clean_title, "^([A-Z][a-z]+)")
    genus_after_pattern <- str_extract(clean_title, "(?<=[:\\s]|^)([A-Z][a-z]+)(?=\\s)")
    genus <- if_else(!is.na(genus_start), genus_start, genus_after_pattern)
    
    # Return binomial if found, otherwise genus, otherwise NA
    if_else(!is.na(binomial), binomial, 
           if_else(!is.na(genus), genus, NA_character_))
}

# Function to search for organism in genome table (precise matching - from script 21)
search_genome_table <- function(organism_name, genome_table, db_name) {
    if(is.null(genome_table) || nrow(genome_table) == 0) {
        return(list(found = FALSE, count = 0, details = NULL, taxids = NULL))
    }
    
    # Handle empty or NA organism names (check first element if vector)
    org_check <- if(length(organism_name) > 0) organism_name[1] else NA_character_
    if(is.null(org_check) || is.na(org_check) || org_check == "" || nchar(org_check) == 0) {
        return(list(found = FALSE, count = 0, details = NULL, taxids = NULL))
    }
    
    # Use first element if vector
    org_name <- if(length(organism_name) > 0) organism_name[1] else organism_name
    
    parts <- strsplit(org_name, " ")[[1]]
    if(length(parts) == 0) {
        return(list(found = FALSE, count = 0, details = NULL, taxids = NULL))
    }
    
    genus <- parts[1]
    species <- if(length(parts) > 1) paste(parts[1], parts[2]) else NULL
    
    available_cols <- colnames(genome_table)
    matches <- data.frame()
    
    # Try exact match first
    if("ncbi_organism_name" %in% available_cols) {
        matches <- genome_table %>%
            filter(tolower(ncbi_organism_name) == tolower(org_name))
        
        if(nrow(matches) == 0 && !is.null(species)) {
            escaped_species <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", species)
            pattern <- paste0("\\b", escaped_species, "\\b")
            matches <- genome_table %>%
                filter(grepl(pattern, ncbi_organism_name, ignore.case = TRUE))
        }
    }
    
    if(nrow(matches) > 0) {
        taxids <- unique(matches$ncbi_species_taxid[!is.na(matches$ncbi_species_taxid)])
        return(list(
            found = TRUE,
            count = nrow(matches),
            details = matches %>%
                select(any_of(c("ncbi_organism_name", "ncbi_species_taxid", "source", "genus", "species"))) %>%
                head(10),
            taxids = taxids
        ))
    }
    
    return(list(found = FALSE, count = 0, details = NULL, taxids = NULL))
}

# Add n_kmers to comparison if not present
if(!("gtdb_n_kmers" %in% colnames(comparison))) {
    if(file.exists(scores_files[["gtdb"]])) {
        gtdb_scores <- fread(scores_files[["gtdb"]], select = c("read_id", "n_kmers"), showProgress = FALSE) %>%
            mutate(read_id = as.character(read_id)) %>%
            rename(gtdb_n_kmers = n_kmers)
        comparison <- comparison %>% left_join(gtdb_scores, by = "read_id")
    } else {
        comparison <- comparison %>% mutate(gtdb_n_kmers = NA)
    }
}
if(!("smdb_n_kmers" %in% colnames(comparison))) {
    if(file.exists(scores_files[["smdb"]])) {
        smdb_scores <- fread(scores_files[["smdb"]], select = c("read_id", "n_kmers"), showProgress = FALSE) %>%
            mutate(read_id = as.character(read_id)) %>%
            rename(smdb_n_kmers = n_kmers)
        comparison <- comparison %>% left_join(smdb_scores, by = "read_id")
    } else {
        comparison <- comparison %>% mutate(smdb_n_kmers = NA)
    }
}
if(!("pluspf_n_kmers" %in% colnames(comparison))) {
    if(file.exists(scores_files[["pluspf"]])) {
        pluspf_scores <- fread(scores_files[["pluspf"]], select = c("read_id", "n_kmers"), showProgress = FALSE) %>%
            mutate(read_id = as.character(read_id)) %>%
            rename(pluspf_n_kmers = n_kmers)
        comparison <- comparison %>% left_join(pluspf_scores, by = "read_id")
    } else {
        comparison <- comparison %>% mutate(pluspf_n_kmers = NA)
    }
}

# K-mer investigation for missing organisms
none_correct_kmer <- blast_data %>%
    filter(!gtdb_correct & !smdb_correct & !pluspf_correct) %>%
    mutate(organism_name = extract_organism_name(title)) %>%
    filter(!is.na(organism_name))

organism_analysis <- none_correct_kmer %>%
    group_by(organism_name, blast_category) %>%
    summarize(
        n_reads = n(),
        mean_identity = mean(identity, na.rm = TRUE),
        mean_bitscore = mean(bitscore, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(desc(n_reads)) %>%
    head(20)

# Check database presence for top organisms
detailed_analysis <- organism_analysis %>%
    filter(!is.na(organism_name) & organism_name != "" & nchar(organism_name) > 0) %>%
    rowwise() %>%
    mutate(
        smdb_check = list(search_genome_table(organism_name, genome_tables[["smdb"]], "smdb")),
        gtdb_check = list(search_genome_table(organism_name, genome_tables[["gtdb"]], "gtdb")),
        pluspf_check = list(search_genome_table(organism_name, genome_tables[["pluspf"]], "pluspf"))
    ) %>%
    ungroup() %>%
    mutate(
        smdb_present = map_lgl(smdb_check, ~ .x$found),
        smdb_count = map_int(smdb_check, ~ .x$count),
        gtdb_present = map_lgl(gtdb_check, ~ .x$found),
        gtdb_count = map_int(gtdb_check, ~ .x$count),
        pluspf_present = map_lgl(pluspf_check, ~ .x$found),
        pluspf_count = map_int(pluspf_check, ~ .x$count)
    ) %>%
    mutate(
        org_reads = map(organism_name, ~{
            none_correct_kmer %>%
                filter(organism_name == .x) %>%
        left_join(comparison, by = "read_id", suffix = c("_blast", ""))
        }),
        smd_assignments = map(org_reads, ~{
            .x %>% filter(!is.na(smdb_taxid)) %>%
                count(smdb_rank, smdb_name, sort = TRUE) %>% head(5)
        }),
        gtdb_assignments = map(org_reads, ~{
            .x %>% filter(!is.na(gtdb_taxid)) %>%
                count(gtdb_rank, gtdb_name, sort = TRUE) %>% head(5)
        }),
        pluspf_assignments = map(org_reads, ~{
            .x %>% filter(!is.na(pluspf_taxid)) %>%
                count(pluspf_rank, pluspf_name, sort = TRUE) %>% head(5)
        }),
        quality_metrics = map(org_reads, ~{
            .x %>% summarize(
                mean_smd_consistency = mean(smdb_consistency, na.rm = TRUE),
                mean_smd_n_kmers = mean(smdb_n_kmers, na.rm = TRUE),
                mean_gtdb_consistency = mean(gtdb_consistency, na.rm = TRUE),
                mean_gtdb_n_kmers = mean(gtdb_n_kmers, na.rm = TRUE),
                mean_pluspf_consistency = mean(pluspf_consistency, na.rm = TRUE),
                mean_pluspf_n_kmers = mean(pluspf_n_kmers, na.rm = TRUE)
            )
        })
    )

# Generate k-mer investigation summary
kmer_summary_lines <- c(
    "K-mer Match Investigation for BLAST-Identified Organisms",
    paste0("Sample: ", sampleID),
    paste0("Date: ", Sys.Date()),
    "",
    "Summary:",
    paste0("Total reads analyzed (no database correct): ", nrow(none_correct_kmer)),
    paste0("Unique organisms identified: ", length(unique(none_correct_kmer$organism_name))),
    "",
    "Top Organisms and Database Presence:",
    ""
)

    kmer_summary_lines <- c(kmer_summary_lines,
    map_chr(1:nrow(detailed_analysis), ~{
        row <- detailed_analysis[.x,]
        paste0(.x, ". ", row$organism_name, " (", row$blast_category, ")",
               "\n   Reads: ", row$n_reads,
               "\n   SMD present: ", ifelse(row$smdb_present, paste0("YES (", row$smdb_count, ")"), "NO"),
               "\n   GTDB present: ", ifelse(row$gtdb_present, paste0("YES (", row$gtdb_count, ")"), "NO"),
               "\n   PlusPF present: ", ifelse(row$pluspf_present, paste0("YES (", row$pluspf_count, ")"), "NO"),
               "\n")
    })
)

kmer_summary_file <- file.path(output_dir, paste0(sampleID, "_kmer_investigation_summary.txt"))
writeLines(kmer_summary_lines, kmer_summary_file)

# Present but misclassified investigation
blast_data_with_comparison <- blast_data %>%
    left_join(comparison, by = "read_id", suffix = c("_blast", ""))

none_correct_present <- blast_data_with_comparison %>%
    filter(!gtdb_correct & !smdb_correct & !pluspf_correct) %>%
    mutate(organism_name = extract_organism_name(title)) %>%
    filter(!is.na(organism_name))

# Check which organisms are present in any of the three databases
present_but_misclassified <- none_correct_present %>%
    group_by(organism_name, blast_category) %>%
    summarize(n_reads = n(), .groups = "drop") %>%
    rowwise() %>%
    mutate(
        smdb_check = list(search_genome_table(organism_name, genome_tables[["smdb"]], "smdb")),
        gtdb_check = list(search_genome_table(organism_name, genome_tables[["gtdb"]], "gtdb")),
        pluspf_check = list(search_genome_table(organism_name, genome_tables[["pluspf"]], "pluspf"))
    ) %>%
    ungroup() %>%
    mutate(
        smdb_present = map_lgl(smdb_check, ~ .x$found),
        smdb_count = map_int(smdb_check, ~ .x$count),
        smdb_taxids = map(smdb_check, ~ .x$taxids),
        gtdb_present = map_lgl(gtdb_check, ~ .x$found),
        gtdb_count = map_int(gtdb_check, ~ .x$count),
        gtdb_taxids = map(gtdb_check, ~ .x$taxids),
        pluspf_present = map_lgl(pluspf_check, ~ .x$found),
        pluspf_count = map_int(pluspf_check, ~ .x$count),
        pluspf_taxids = map(pluspf_check, ~ .x$taxids),
        any_db_present = smdb_present | gtdb_present | pluspf_present,
        all_expected_taxids = map(1:n(), ~{
            unique(c(smdb_taxids[[.x]], gtdb_taxids[[.x]], pluspf_taxids[[.x]]))
        })
    ) %>%
    filter(any_db_present == TRUE) %>%
    arrange(desc(n_reads))

if(nrow(present_but_misclassified) > 0) {
    # Detailed analysis for each organism
    present_detailed_results <- present_but_misclassified %>%
        mutate(
            org_reads_df = map(organism_name, ~{
                none_correct_present %>%
                    filter(organism_name == .x) %>%
            select(read_id, title, 
                   smdb_taxid, smdb_name, smdb_rank, smdb_consistency, smdb_multiplicity, smdb_entropy,
                   gtdb_taxid, gtdb_name, gtdb_rank, gtdb_consistency, gtdb_multiplicity, gtdb_entropy,
                   pluspf_taxid, pluspf_name, pluspf_rank, pluspf_consistency, pluspf_multiplicity, pluspf_entropy,
                   any_of(c("smdb_n_kmers", "gtdb_n_kmers", "pluspf_n_kmers", 
                           "smdb_is_confident", "gtdb_is_confident", "pluspf_is_confident")))
            }),
            total_correct = map2_int(org_reads_df, all_expected_taxids, ~{
                org_df <- .x
                expected <- .y
                correct_smd <- org_df %>% filter(smdb_taxid %in% expected) %>% pull(read_id)
                correct_gtdb <- org_df %>% filter(gtdb_taxid %in% expected) %>% pull(read_id)
                correct_pluspf <- org_df %>% filter(pluspf_taxid %in% expected) %>% pull(read_id)
                length(unique(c(correct_smd, correct_gtdb, correct_pluspf)))
            }),
            smd_assignments = map2(org_reads_df, smdb_present, ~{
                if(.y) {
                    .x %>% filter(!is.na(smdb_taxid)) %>%
                        count(smdb_rank, smdb_name, smdb_taxid, sort = TRUE) %>% head(10)
                } else data.frame()
            }),
            gtdb_assignments = map2(org_reads_df, gtdb_present, ~{
                if(.y) {
                    .x %>% filter(!is.na(gtdb_taxid)) %>%
                        count(gtdb_rank, gtdb_name, gtdb_taxid, sort = TRUE) %>% head(10)
                } else data.frame()
            }),
            pluspf_assignments = map2(org_reads_df, pluspf_present, ~{
                if(.y) {
                    .x %>% filter(!is.na(pluspf_taxid)) %>%
                        count(pluspf_rank, pluspf_name, pluspf_taxid, sort = TRUE) %>% head(10)
                } else data.frame()
            }),
            quality_metrics = map(org_reads_df, ~{
                .x %>% summarize(
                mean_smd_consistency = mean(smdb_consistency, na.rm = TRUE),
                mean_smd_n_kmers = mean(smdb_n_kmers, na.rm = TRUE),
                mean_smd_multiplicity = mean(smdb_multiplicity, na.rm = TRUE),
                mean_smd_entropy = mean(smdb_entropy, na.rm = TRUE),
                n_smd_classified = sum(!is.na(smdb_taxid)),
                n_smd_confident = sum(smdb_is_confident, na.rm = TRUE),
                mean_gtdb_consistency = mean(gtdb_consistency, na.rm = TRUE),
                mean_gtdb_n_kmers = mean(gtdb_n_kmers, na.rm = TRUE),
                mean_gtdb_multiplicity = mean(gtdb_multiplicity, na.rm = TRUE),
                mean_gtdb_entropy = mean(gtdb_entropy, na.rm = TRUE),
                n_gtdb_classified = sum(!is.na(gtdb_taxid)),
                n_gtdb_confident = sum(gtdb_is_confident, na.rm = TRUE),
                mean_pluspf_consistency = mean(pluspf_consistency, na.rm = TRUE),
                mean_pluspf_n_kmers = mean(pluspf_n_kmers, na.rm = TRUE),
                mean_pluspf_multiplicity = mean(pluspf_multiplicity, na.rm = TRUE),
                mean_pluspf_entropy = mean(pluspf_entropy, na.rm = TRUE),
                n_pluspf_classified = sum(!is.na(pluspf_taxid)),
                n_pluspf_confident = sum(pluspf_is_confident, na.rm = TRUE)
            )
            }),
            expected_taxids_str = map_chr(all_expected_taxids, ~paste(.x, collapse = ";"))
        )
    
    # Generate present but misclassified summary
    present_summary_lines <- c(
        "Organisms Present But Misclassified: Detailed Analysis",
        paste0("Sample: ", sampleID),
        paste0("Date: ", Sys.Date()),
        "",
        paste0("Total organisms present in at least one database but misclassified: ", nrow(present_detailed_results)),
        ""
    )
    
    # Generate summary lines using map
    org_summary_lines <- map_chr(1:nrow(present_detailed_results), ~{
        row <- present_detailed_results[.x,]
        present_dbs <- c()
        if(row$smdb_present) present_dbs <- c(present_dbs, paste0("SMD (", row$smdb_count, " genome(s))"))
        if(row$gtdb_present) present_dbs <- c(present_dbs, paste0("GTDB (", row$gtdb_count, " genome(s))"))
        if(row$pluspf_present) present_dbs <- c(present_dbs, paste0("PlusPF (", row$pluspf_count, " genome(s))"))
        
        lines <- c(
            paste0(.x, ". ", row$organism_name, " (", row$blast_category, ")"),
            paste0("   Reads: ", row$n_reads),
            paste0("   Present in: ", paste(present_dbs, collapse = ", ")),
            paste0("   Expected taxid(s): ", row$expected_taxids_str),
            paste0("   Reads assigned to correct taxid: ", row$total_correct, " / ", nrow(row$org_reads_df[[1]]))
        )
        
        # Add assignments
        if(row$smdb_present && nrow(row$smd_assignments[[1]]) > 0) {
            lines <- c(lines, "   SMD Assignments:",
                paste0("     ", row$smd_assignments[[1]]$smdb_rank, " | ", 
                      row$smd_assignments[[1]]$smdb_name, " (taxid: ", row$smd_assignments[[1]]$smdb_taxid, 
                      ") - ", row$smd_assignments[[1]]$n, " read(s)"))
        }
        if(row$gtdb_present && nrow(row$gtdb_assignments[[1]]) > 0) {
            lines <- c(lines, "   GTDB Assignments:",
                paste0("     ", row$gtdb_assignments[[1]]$gtdb_rank, " | ", 
                      row$gtdb_assignments[[1]]$gtdb_name, " (taxid: ", row$gtdb_assignments[[1]]$gtdb_taxid, 
                      ") - ", row$gtdb_assignments[[1]]$n, " read(s)"))
        }
        if(row$pluspf_present && nrow(row$pluspf_assignments[[1]]) > 0) {
            lines <- c(lines, "   PlusPF Assignments:",
                paste0("     ", row$pluspf_assignments[[1]]$pluspf_rank, " | ", 
                      row$pluspf_assignments[[1]]$pluspf_name, " (taxid: ", row$pluspf_assignments[[1]]$pluspf_taxid, 
                      ") - ", row$pluspf_assignments[[1]]$n, " read(s)"))
        }
        
        # Quality metrics
        qm <- row$quality_metrics[[1]]
        lines <- c(lines, "   Quality Metrics:")
        if(row$smdb_present) {
            lines <- c(lines,
                paste0("     SMD: consistency=", round(qm$mean_smd_consistency, 3),
                      ", n_kmers=", round(qm$mean_smd_n_kmers, 1),
                      ", multiplicity=", round(qm$mean_smd_multiplicity, 2),
                      ", entropy=", round(qm$mean_smd_entropy, 3)),
                paste0("     SMD classified: ", qm$n_smd_classified, " / ", nrow(row$org_reads_df[[1]]),
                      " (Confident: ", qm$n_smd_confident, ")"))
        }
        if(row$gtdb_present) {
            lines <- c(lines,
                paste0("     GTDB: consistency=", round(qm$mean_gtdb_consistency, 3),
                      ", n_kmers=", round(qm$mean_gtdb_n_kmers, 1),
                      ", multiplicity=", round(qm$mean_gtdb_multiplicity, 2),
                      ", entropy=", round(qm$mean_gtdb_entropy, 3)),
                paste0("     GTDB classified: ", qm$n_gtdb_classified, " / ", nrow(row$org_reads_df[[1]]),
                      " (Confident: ", qm$n_gtdb_confident, ")"))
        }
        if(row$pluspf_present) {
            lines <- c(lines,
                paste0("     PlusPF: consistency=", round(qm$mean_pluspf_consistency, 3),
                      ", n_kmers=", round(qm$mean_pluspf_n_kmers, 1),
                      ", multiplicity=", round(qm$mean_pluspf_multiplicity, 2),
                      ", entropy=", round(qm$mean_pluspf_entropy, 3)),
                paste0("     PlusPF classified: ", qm$n_pluspf_classified, " / ", nrow(row$org_reads_df[[1]]),
                      " (Confident: ", qm$n_pluspf_confident, ")"))
        }
        paste(lines, collapse = "\n")
    })
    
    present_summary_lines <- c(present_summary_lines, org_summary_lines, "")
    
    present_summary_file <- file.path(output_dir, paste0(sampleID, "_present_but_misclassified_analysis.txt"))
    writeLines(present_summary_lines, present_summary_file)
    
    # Save detailed CSV
    present_detailed_df <- present_detailed_results %>%
        select(organism_name, blast_category, n_reads, expected_taxids_str, 
               smdb_present, gtdb_present, pluspf_present, org_reads_df) %>%
        unnest(org_reads_df) %>%
        mutate(
            expected_taxids_list = map(expected_taxids_str, ~strsplit(.x, ";")[[1]]),
            assigned_to_correct_taxid_smd = map2_lgl(smdb_taxid, expected_taxids_list, ~.x %in% .y),
            assigned_to_correct_taxid_gtdb = map2_lgl(gtdb_taxid, expected_taxids_list, ~.x %in% .y),
            assigned_to_correct_taxid_pluspf = map2_lgl(pluspf_taxid, expected_taxids_list, ~.x %in% .y),
            assigned_to_correct_taxid_any = assigned_to_correct_taxid_smd | assigned_to_correct_taxid_gtdb | assigned_to_correct_taxid_pluspf
        ) %>%
        select(-expected_taxids_list) %>%
        rename(expected_taxid = expected_taxids_str)
    
    present_csv_file <- file.path(output_dir, paste0(sampleID, "_present_but_misclassified_detailed.csv"))
    write_csv(present_detailed_df, present_csv_file)
}

# Get sequencing depth
seq_depth_file <- file.path(output_dir, "seq_depth_df.rds")
if(!file.exists(seq_depth_file)) {
    seq_depth_file <- "data/NEON_metagenome_classification/seq_depth_df.rds"
}

total_reads <- NULL
if(file.exists(seq_depth_file)) {
    seq_depth_df <- readRDS(seq_depth_file)
    seq_depth_row <- seq_depth_df %>% filter(sampleID == !!sampleID)
    if(nrow(seq_depth_row) > 0) {
        total_reads <- as.numeric(seq_depth_row$seq_depth[1])
    }
}

# Calculate from kreport if needed
if(is.null(total_reads) || is.na(total_reads)) {
    source("scripts/custom_pavian.r")
    source("https://raw.githubusercontent.com/bhattlab/kraken2_classification/master/scripts/process_classification.R")
    
    kreport_dirs <- c("data/NEON_metagenome_classification/01_kraken_output",
                      "data/classification/01_kraken_output",
                      "data/NEON_metagenome_classification/02_bracken_output",
                      "data/classification/02_bracken_output",
                      "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/01_kraken_output",
                      "/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/02_bracken_output")
    
    kreport_file <- NULL
    for(dir in kreport_dirs[dir.exists(kreport_dirs)]) {
        pattern <- paste0("*", sampleID, "*_kraken.kreport")
        files <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
        if(length(files) > 0) {
            kreport_file <- files[1]
            break
        }
    }
    
    if(!is.null(kreport_file) && file.exists(kreport_file)) {
        tryCatch({
            my_report <- fread_report(kreport_file) %>% as.data.frame()
            rownames(my_report)[rownames(my_report) == "r_root"] <- "-_root"
            my_report <- my_report[!duplicated(my_report$name),]
            row.names(my_report) <- my_report[["name"]]
            
            unclassified_keys <- c("u_unclassified", "unclassified", "U_unclassified")
            root_keys <- c("r_root", "root", "R_root", "-_root")
            unclassified_key <- intersect(unclassified_keys, rownames(my_report))[1]
            root_key <- intersect(root_keys, rownames(my_report))[1]
            
            if(!is.na(unclassified_key) && !is.na(root_key)) {
                unclassified_reads <- as.numeric(my_report[unclassified_key, "cladeReads"])
                root_reads <- as.numeric(my_report[root_key, "cladeReads"])
                total_reads <- unclassified_reads + root_reads
            }
        }, error = function(e) {
            warning("Error calculating from kreport: ", e$message)
        })
    }
}

# Fallback to comparison file
if(is.null(total_reads) || is.na(total_reads)) {
    total_reads <- nrow(comparison)
    if(any(grepl("read_id", comparison$read_id, ignore.case = TRUE))) {
        total_reads <- total_reads - 1
    }
}

# Calculate classification rates
classification_stats <- comparison %>%
    summarize(
        gtdb_classified = sum(!is.na(gtdb_taxid)),
        smdb_classified = sum(!is.na(smdb_taxid)),
        pluspf_classified = sum(!is.na(pluspf_taxid)),
        gtdb_confident = sum(gtdb_is_confident, na.rm = TRUE),
        smdb_confident = sum(smdb_is_confident, na.rm = TRUE),
        pluspf_confident = sum(pluspf_is_confident, na.rm = TRUE)
    )

# BLAST verification results
blast_with_comparison <- blast_data %>%
    left_join(comparison, by = "read_id", suffix = c("_blast", ""))

blast_stats <- blast_with_comparison %>%
    filter(!is.na(accession) & !accession %in% c("no_hits", "blast_failed")) %>%
    summarize(
        total_blast = n(),
        gtdb_correct_count = sum(gtdb_correct, na.rm = TRUE),
        smdb_correct_count = sum(smdb_correct, na.rm = TRUE),
        pluspf_correct_count = sum(pluspf_correct, na.rm = TRUE),
        none_correct = sum(!gtdb_correct & !smdb_correct & !pluspf_correct, na.rm = TRUE),
        at_least_one_correct = sum((gtdb_correct | smdb_correct | pluspf_correct), na.rm = TRUE)
    )

# Analysis by conflict type
homo_sapiens_conflicts <- blast_with_comparison %>%
    filter(conflict_type == "homo_sapiens")

fungal_bacteria_conflicts <- blast_with_comparison %>%
    filter(conflict_type == "fungal_bacteria")

smdb_confident_fungi_verification <- blast_with_comparison %>%
    filter(conflict_type == "smdb_confident_fungi")

# Calculate conflict statistics
conflict_stats <- list()
if(nrow(homo_sapiens_conflicts) > 0) {
    conflict_stats[["homo_sapiens"]] <- homo_sapiens_conflicts %>%
        summarize(
            total = n(),
            pluspf_confident = sum(pluspf_is_confident, na.rm = TRUE),
            smdb_confident = sum(smdb_is_confident, na.rm = TRUE),
            none_correct = sum(!gtdb_correct & !smdb_correct & !pluspf_correct, na.rm = TRUE)
        )
}

if(nrow(fungal_bacteria_conflicts) > 0) {
    conflict_stats[["fungal_bacteria"]] <- fungal_bacteria_conflicts %>%
        summarize(
            total = n(),
            smdb_confident = sum(smdb_is_confident, na.rm = TRUE),
            gtdb_confident = sum(gtdb_is_confident, na.rm = TRUE),
            smdb_correct = sum(smdb_correct, na.rm = TRUE),
            gtdb_correct = sum(gtdb_correct, na.rm = TRUE),
            none_correct = sum(!gtdb_correct & !smdb_correct & !pluspf_correct, na.rm = TRUE)
        )
}

if(nrow(smdb_confident_fungi_verification) > 0) {
    conflict_stats[["smdb_confident_fungi"]] <- smdb_confident_fungi_verification %>%
        summarize(
            total = n(),
            smdb_confident = sum(smdb_is_confident, na.rm = TRUE),
            smdb_correct = sum(smdb_correct, na.rm = TRUE),
            gtdb_correct = sum(gtdb_correct, na.rm = TRUE),
            pluspf_correct = sum(pluspf_correct, na.rm = TRUE),
            none_correct = sum(!gtdb_correct & !smdb_correct & !pluspf_correct, na.rm = TRUE)
        )
}

# Analysis by kingdom/BLAST category
    kingdom_analysis <- blast_with_comparison %>%
        filter(!is.na(blast_category)) %>%
        group_by(blast_category) %>%
        summarize(
            n_reads = n(),
            gtdb_correct = sum(gtdb_correct, na.rm = TRUE),
            smdb_correct = sum(smdb_correct, na.rm = TRUE),
            pluspf_correct = sum(pluspf_correct, na.rm = TRUE),
            none_correct = sum(!gtdb_correct & !smdb_correct & !pluspf_correct, na.rm = TRUE),
            gtdb_confident = sum(gtdb_is_confident, na.rm = TRUE),
            smdb_confident = sum(smdb_is_confident, na.rm = TRUE),
            pluspf_confident = sum(pluspf_is_confident, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(desc(n_reads))
    
# Relative abundance consequences
wrong_assignments <- NULL
if(!is.null(bracken_data)) {
    none_correct_reads <- blast_with_comparison %>%
        filter(!gtdb_correct & !smdb_correct & !pluspf_correct & !is.na(blast_category))
    
        wrong_assignments <- none_correct_reads %>%
            filter(!is.na(smdb_taxid)) %>%
            count(smdb_taxid, smdb_name, sort = TRUE) %>%
        head(10) %>%
        left_join(bracken_data, by = c("smdb_taxid" = "taxonomy_id")) %>%
        mutate(
            bracken_reads = ifelse(is.na(new_est_reads), 0, new_est_reads),
            bracken_fraction = ifelse(is.na(fraction_total_reads), 0, fraction_total_reads * 100)
        )
}

# Filter effectiveness for PlusPF "Homo sapiens" reads
hs_none_correct <- NULL
hs_correlation_stats <- NULL
if(nrow(homo_sapiens_conflicts) > 0) {
    hs_none_correct <- homo_sapiens_conflicts %>%
        filter(!gtdb_correct & !smdb_correct & !pluspf_correct) %>%
        mutate(
            pluspf_consistency_num = as.numeric(pluspf_consistency),
            blast_identity_num = as.numeric(identity)
        ) %>%
        filter(!is.na(pluspf_consistency_num) & !is.na(blast_identity_num))
    
    if(nrow(hs_none_correct) > 0) {
        hs_correlation_stats <- hs_none_correct %>%
            summarize(
                correlation = cor(pluspf_consistency_num, blast_identity_num, use = "complete.obs"),
                n_reads = n(),
                miscalibrated_count = sum(pluspf_consistency_num >= 0.95 & blast_identity_num < 85, na.rm = TRUE)
            )
        
        if(hs_correlation_stats$n_reads > 3) {
            cor_test <- cor.test(hs_none_correct$pluspf_consistency_num,
                                hs_none_correct$blast_identity_num)
            hs_correlation_stats$p_value <- cor_test$p.value
        }
    }
}

# Impact assessment and categorization
blast_with_comparison_cat <- blast_data %>%
    left_join(comparison, by = "read_id", suffix = c("_blast", "")) %>%
    mutate(
        organism_name = extract_organism_name(title),
        identity_num = identity,
        bitscore_num = bitscore
    )

# Check presence for each unique organism
unique_organisms_cat <- blast_with_comparison_cat %>%
    filter(!gtdb_correct & !smdb_correct & !pluspf_correct & !is.na(organism_name)) %>%
    distinct(organism_name) %>%
    rowwise() %>%
    mutate(
        smdb_present = search_genome_table(organism_name, genome_tables[["smdb"]], "smdb")$found,
        gtdb_present = search_genome_table(organism_name, genome_tables[["gtdb"]], "gtdb")$found,
        pluspf_present = search_genome_table(organism_name, genome_tables[["pluspf"]], "pluspf")$found
    ) %>%
    ungroup()

# Categorize reads
categorized_reads <- blast_with_comparison_cat %>%
    filter(!gtdb_correct & !smdb_correct & !pluspf_correct) %>%
    left_join(unique_organisms_cat, by = "organism_name") %>%
    mutate(
        smdb_present = ifelse(is.na(smdb_present), FALSE, smdb_present),
        gtdb_present = ifelse(is.na(gtdb_present), FALSE, gtdb_present),
        pluspf_present = ifelse(is.na(pluspf_present), FALSE, pluspf_present),
        any_db_present = smdb_present | gtdb_present | pluspf_present,
        category = case_when(
            is.na(organism_name) ~ "Unclassifiable (no organism name)",
            !any_db_present ~ "Missing from databases",
            any_db_present ~ "Present but misclassified",
            TRUE ~ "Unknown"
        ),
        blast_quality = case_when(
            identity_num >= blast_quality_thresholds$high_identity & 
            bitscore_num >= blast_quality_thresholds$high_bitscore ~ "High",
            identity_num >= blast_quality_thresholds$medium_identity & 
            bitscore_num >= blast_quality_thresholds$medium_bitscore ~ "Medium",
            TRUE ~ "Low"
        ),
        was_classified = !is.na(smdb_taxid) | !is.na(gtdb_taxid) | !is.na(pluspf_taxid),
        passed_filter = (smdb_is_confident == TRUE) | (gtdb_is_confident == TRUE) | (pluspf_is_confident == TRUE)
    )

category_summary <- categorized_reads %>%
    group_by(category, blast_quality) %>%
    summarize(
        n_reads = n(),
        mean_identity = mean(identity_num, na.rm = TRUE),
        mean_bitscore = mean(bitscore_num, na.rm = TRUE),
        classified_count = sum(was_classified),
        filtered_count = sum(passed_filter),
        .groups = "drop"
    ) %>%
    arrange(desc(n_reads))

# Impact on downstream interpretation
impact_summary <- NULL
if(!is.null(bracken_data) && !is.null(total_reads)) {
    wrong_assignments_filtered <- categorized_reads %>%
        filter(passed_filter & !is.na(smdb_taxid)) %>%
        count(smdb_taxid, smdb_name, sort = TRUE) %>%
            head(20) %>%
        left_join(bracken_data, by = c("smdb_taxid" = "taxonomy_id")) %>%
            mutate(
            bracken_reads = ifelse(is.na(new_est_reads), 0, new_est_reads),
            bracken_fraction = ifelse(is.na(fraction_total_reads), 0, fraction_total_reads * 100),
            misclassification_impact = ifelse(bracken_reads > 0, (n / bracken_reads * 100), NA_real_)
        )
}

# Fixable vs unresolvable categorization
fixable_expansion <- categorized_reads %>%
    filter(category == "Missing from databases") %>%
    group_by(organism_name, blast_category) %>%
    summarize(
        n_reads = n(),
        mean_identity = mean(identity_num, na.rm = TRUE),
        mean_bitscore = mean(bitscore_num, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(desc(n_reads))

fixable_algorithm <- categorized_reads %>%
    filter(category == "Present but misclassified") %>%
    group_by(organism_name, blast_category, smdb_present, gtdb_present, pluspf_present) %>%
    summarize(
        n_reads = n(),
        mean_identity = mean(identity_num, na.rm = TRUE),
        mean_bitscore = mean(bitscore_num, na.rm = TRUE),
        classified_count = sum(was_classified),
        filtered_count = sum(passed_filter),
        .groups = "drop"
    ) %>%
    arrange(desc(n_reads))

unresolvable <- categorized_reads %>%
    filter(blast_quality == "Low") %>%
    group_by(blast_category) %>%
    summarize(
        n_reads = n(),
        mean_identity = mean(identity_num, na.rm = TRUE),
        mean_bitscore = mean(bitscore_num, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(desc(n_reads))

# Save detailed categorization
categorization_file <- file.path(output_dir, paste0(sampleID, "_misclassification_categorization.csv"))
write_csv(categorized_reads, categorization_file)

# Generate impact assessment summary
impact_summary_lines <- c(
    "Misclassification Impact Analysis",
    paste0("Sample: ", sampleID),
    paste0("Date: ", Sys.Date()),
    "",
    "Summary:",
    paste0("Total misclassified reads (none correct): ", nrow(categorized_reads)),
    "",
    "Fixable with Database Expansion:",
    paste0("Total reads: ", sum(fixable_expansion$n_reads)),
    paste0("Unique organisms: ", nrow(fixable_expansion)),
    "",
    "Top organisms:",
    if(nrow(fixable_expansion) > 0) {
    paste0("  ", fixable_expansion$organism_name[1:min(20, nrow(fixable_expansion))], 
           " (", fixable_expansion$blast_category[1:min(20, nrow(fixable_expansion))], 
               "): ", fixable_expansion$n_reads[1:min(20, nrow(fixable_expansion))], " reads")
    } else "  None",
    "",
    "Potentially Fixable with Better Algorithms:",
    paste0("Total reads: ", sum(fixable_algorithm$n_reads)),
    paste0("Unique organisms: ", nrow(fixable_algorithm)),
    "",
    "Likely Unresolvable (Low Quality):",
    paste0("Total reads: ", sum(unresolvable$n_reads)),
    "",
    "Proportions:"
)

if(!is.null(total_reads)) {
    impact_summary_lines <- c(impact_summary_lines,
        paste0("Of total sequencing depth (", total_reads, " reads):"),
        paste0("  Misclassified: ", round(nrow(categorized_reads)/total_reads*100, 4), "%"),
        paste0("  Fixable with expansion: ", round(sum(fixable_expansion$n_reads)/total_reads*100, 4), "%"),
        paste0("  Potentially fixable: ", round(sum(fixable_algorithm$n_reads)/total_reads*100, 4), "%"),
        paste0("  Likely unresolvable: ", round(sum(unresolvable$n_reads)/total_reads*100, 4), "%")
    )
}

impact_summary_file <- file.path(output_dir, paste0(sampleID, "_misclassification_impact_analysis.txt"))
writeLines(impact_summary_lines, impact_summary_file)

# Generate summary tables (CSV format instead of narrative text)

# Generate summary tables (CSV format instead of narrative text)
# 1. Overall classification statistics
overall_stats <- tibble(
    metric = c("total_reads", "gtdb_classified", "smdb_classified", "pluspf_classified",
               "gtdb_confident", "smdb_confident", "pluspf_confident"),
    count = c(total_reads,
              classification_stats$gtdb_classified,
              classification_stats$smdb_classified,
              classification_stats$pluspf_classified,
              classification_stats$gtdb_confident,
              classification_stats$smdb_confident,
              classification_stats$pluspf_confident),
    pct_of_total = c(100,
                     round(classification_stats$gtdb_classified/total_reads*100, 2),
                     round(classification_stats$smdb_classified/total_reads*100, 2),
                     round(classification_stats$pluspf_classified/total_reads*100, 2),
                     round(classification_stats$gtdb_confident/total_reads*100, 2),
                     round(classification_stats$smdb_confident/total_reads*100, 2),
                     round(classification_stats$pluspf_confident/total_reads*100, 2)),
    pct_of_classified = c(NA_real_,
                          NA_real_, NA_real_, NA_real_,
                          round(classification_stats$gtdb_confident/classification_stats$gtdb_classified*100, 2),
                          round(classification_stats$smdb_confident/classification_stats$smdb_classified*100, 2),
                          round(classification_stats$pluspf_confident/classification_stats$pluspf_classified*100, 2))
)

# 2. BLAST verification statistics
blast_verification_stats <- tibble(
    metric = c("total_blast", "none_correct", "at_least_one_correct",
               "gtdb_correct", "smdb_correct", "pluspf_correct"),
    count = c(blast_stats$total_blast,
              blast_stats$none_correct,
              blast_stats$at_least_one_correct,
              blast_stats$gtdb_correct_count,
              blast_stats$smdb_correct_count,
              blast_stats$pluspf_correct_count),
    pct_of_blast = c(100,
                     round(blast_stats$none_correct/blast_stats$total_blast*100, 2),
                     round(blast_stats$at_least_one_correct/blast_stats$total_blast*100, 2),
                     round(blast_stats$gtdb_correct_count/blast_stats$total_blast*100, 2),
                     round(blast_stats$smdb_correct_count/blast_stats$total_blast*100, 2),
                     round(blast_stats$pluspf_correct_count/blast_stats$total_blast*100, 2))
)

# 3. Conflict type statistics
conflict_stats_table <- NULL
if(length(conflict_stats) > 0) {
    conflict_stats_table <- map_dfr(names(conflict_stats), function(conflict_name) {
        stats <- conflict_stats[[conflict_name]]
        if(conflict_name == "homo_sapiens") {
            tibble(
                conflict_type = conflict_name,
                total_reads = stats$total,
                pluspf_confident = stats$pluspf_confident,
                pluspf_confident_pct = round(stats$pluspf_confident/stats$total*100, 2),
                smdb_confident = stats$smdb_confident,
                smdb_confident_pct = round(stats$smdb_confident/stats$total*100, 2),
                none_correct = stats$none_correct,
                none_correct_pct = round(stats$none_correct/stats$total*100, 2)
            )
        } else if(conflict_name == "fungal_bacteria") {
            tibble(
                conflict_type = conflict_name,
                total_reads = stats$total,
                smdb_confident = stats$smdb_confident,
                smdb_confident_pct = round(stats$smdb_confident/stats$total*100, 2),
                gtdb_confident = stats$gtdb_confident,
                gtdb_confident_pct = round(stats$gtdb_confident/stats$total*100, 2),
                smdb_correct = stats$smdb_correct,
                smdb_correct_pct = round(stats$smdb_correct/stats$total*100, 2),
                gtdb_correct = stats$gtdb_correct,
                gtdb_correct_pct = round(stats$gtdb_correct/stats$total*100, 2),
                none_correct = stats$none_correct,
                none_correct_pct = round(stats$none_correct/stats$total*100, 2)
            )
        }
    })
}

# 4. Kingdom/BLAST category analysis
kingdom_analysis_table <- kingdom_analysis %>%
    mutate(
        gtdb_correct_pct = round(gtdb_correct/n_reads*100, 2),
        smdb_correct_pct = round(smdb_correct/n_reads*100, 2),
        pluspf_correct_pct = round(pluspf_correct/n_reads*100, 2),
        none_correct_pct = round(none_correct/n_reads*100, 2),
        gtdb_confident_pct = round(gtdb_confident/n_reads*100, 2),
        smdb_confident_pct = round(smdb_confident/n_reads*100, 2),
        pluspf_confident_pct = round(pluspf_confident/n_reads*100, 2)
    )

# Save all tables to CSV files
overall_stats_file <- file.path(output_dir, paste0(sampleID, "_overall_classification_stats.csv"))
write_csv(overall_stats, overall_stats_file)
cat("\nSaved overall classification statistics to:", overall_stats_file, "\n")

blast_stats_file <- file.path(output_dir, paste0(sampleID, "_blast_verification_stats.csv"))
write_csv(blast_verification_stats, blast_stats_file)
cat("Saved BLAST verification statistics to:", blast_stats_file, "\n")

if(!is.null(conflict_stats_table) && nrow(conflict_stats_table) > 0) {
    conflict_stats_file <- file.path(output_dir, paste0(sampleID, "_conflict_type_stats.csv"))
    write_csv(conflict_stats_table, conflict_stats_file)
    cat("Saved conflict type statistics to:", conflict_stats_file, "\n")
}

if(nrow(kingdom_analysis_table) > 0) {
    kingdom_stats_file <- file.path(output_dir, paste0(sampleID, "_kingdom_blast_category_stats.csv"))
    write_csv(kingdom_analysis_table, kingdom_stats_file)
    cat("Saved kingdom/BLAST category statistics to:", kingdom_stats_file, "\n")
}

# Main summary file (overall stats as primary table)
summary_file <- file.path(output_dir, paste0(sampleID, "_misclassification_analysis.csv"))
write_csv(overall_stats %>% mutate(sample_id = sampleID, date = as.character(Sys.Date())), summary_file)

# Enhance blast_vs_databases_comparison.csv with useful categorization columns
blast_comparison_file <- file.path(output_dir, paste0(sampleID, "_blast_vs_databases_comparison.csv"))
if(file.exists(blast_comparison_file)) {
    cat("\n=== Enhancing blast_vs_databases_comparison.csv with categorization columns ===\n")
    
    # Read the cleaned comparison file
    blast_comparison <- read_csv(blast_comparison_file, show_col_types = FALSE) %>%
        mutate(read_id = as.character(read_id))
    
    # Add organism name extraction
    blast_comparison <- blast_comparison %>%
        mutate(organism_name = extract_organism_name(title))
    
    # Check database presence for unique organisms
    unique_organisms <- blast_comparison %>%
        filter(!is.na(organism_name)) %>%
        distinct(organism_name) %>%
        rowwise() %>%
        mutate(
            smdb_present = search_genome_table(organism_name, genome_tables[["smdb"]], "smdb")$found,
            gtdb_present = search_genome_table(organism_name, genome_tables[["gtdb"]], "gtdb")$found,
            pluspf_present = search_genome_table(organism_name, genome_tables[["pluspf"]], "pluspf")$found
        ) %>%
        ungroup()
    
    # Add categorization columns
    blast_comparison_enhanced <- blast_comparison %>%
        left_join(unique_organisms, by = "organism_name") %>%
        mutate(
            smdb_present = ifelse(is.na(.data$smdb_present), FALSE, .data$smdb_present),
            gtdb_present = ifelse(is.na(.data$gtdb_present), FALSE, .data$gtdb_present),
            pluspf_present = ifelse(is.na(.data$pluspf_present), FALSE, .data$pluspf_present),
            any_db_present = .data$smdb_present | .data$gtdb_present | .data$pluspf_present,
            category = case_when(
                is.na(organism_name) ~ "Unclassifiable (no organism name)",
                !any_db_present ~ "Missing from databases",
                any_db_present ~ "Present but misclassified",
                TRUE ~ "Unknown"
            ),
            blast_quality = case_when(
                identity >= blast_quality_thresholds$high_identity & 
                bitscore >= blast_quality_thresholds$high_bitscore ~ "High",
                identity >= blast_quality_thresholds$medium_identity & 
                bitscore >= blast_quality_thresholds$medium_bitscore ~ "Medium",
                TRUE ~ "Low"
            ),
            was_classified = !is.na(smdb_taxid) | !is.na(gtdb_taxid) | !is.na(pluspf_taxid),
            passed_filter = (smdb_is_confident == TRUE) | (gtdb_is_confident == TRUE) | (pluspf_is_confident == TRUE)
        ) %>%
        # Reorder columns: BLAST info, categorization, database assignments, correctness
        select(
            # BLAST information
            read_id, accession, title, evalue, bitscore, identity, align_length,
            # Conflict and categorization
            conflict_type, blast_category, reliable_hit,
            # Categorization (new)
            organism_name, category, blast_quality, 
            smdb_present, gtdb_present, pluspf_present, any_db_present,
            was_classified, passed_filter,
            # GTDB assignment
            gtdb_taxid, gtdb_name, gtdb_rank, gtdb_is_confident,
            gtdb_consistency, gtdb_multiplicity, gtdb_entropy,
            # SMD assignment
            smdb_taxid, smdb_name, smdb_rank, smdb_is_confident,
            smdb_consistency, smdb_multiplicity, smdb_entropy,
            # PlusPF assignment
            pluspf_taxid, pluspf_name, pluspf_rank, pluspf_is_confident,
            pluspf_consistency, pluspf_multiplicity, pluspf_entropy,
            # Correctness flags
            gtdb_correct, smdb_correct, pluspf_correct
        )
    
    # Write enhanced version back to the same file
    write_csv(blast_comparison_enhanced, blast_comparison_file)
    cat("Enhanced blast_vs_databases_comparison.csv with categorization columns\n")
    cat("Total columns:", ncol(blast_comparison_enhanced), "\n")
    cat("Total rows:", nrow(blast_comparison_enhanced), "\n")
}
cat("Saved main misclassification analysis table to:", summary_file, "\n")
cat("Note: Additional detailed tables saved separately (blast_verification, conflict_type, kingdom_analysis)\n")
