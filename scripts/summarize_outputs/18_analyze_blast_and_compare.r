#!/usr/bin/env Rscript
# 18: Analyze BLAST results and compare to database assignments
# Combines functionality of old scripts 18+19:
#   - Analyzes BLAST result files from script 17
#   - Categorizes BLAST hits
#   - Determines which database (if any) was correct for each BLASTed read
#
# Usage: Rscript scripts/summarize_outputs/18_analyze_blast_and_compare.r [sampleID]
#
# Input:  
#   - BLAST result files (from script 17)
#   - taxonomic_assignment_comparison_{sampleID}.csv (from script 15)
# Output: 
#   - {sampleID}_blast_vs_databases_comparison.csv
#   - {sampleID}_blast_vs_databases_summary.txt

library(tidyverse)

# Configuration
sampleID <- if(length(commandArgs(trailingOnly = TRUE)) > 0) {
    commandArgs(trailingOnly = TRUE)[1]
} else {
    "ORNL_046-O-20170621-COMP"
}

output_dir <- "data/classification/analysis_files"
comparison_file <- file.path(output_dir, paste0("taxonomic_assignment_comparison_", sampleID, ".csv"))

# Taxonomic patterns
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")
fungal_pattern <- paste0("(", paste(fungal_phyla, collapse = "|"), ")")
fungal_genera <- c("Rhizophagus", "Xylographa", "Orbilia", "Arthroderma", "Elaphomyces", 
                   "Hygrocybe", "Anguillospora", "Lipomyces", "Leotia", "Glarea", "Thuemenidium")
fungal_keywords <- paste(c("\\bfungi\\b", "\\bFungi\\b", "\\bfungus\\b", fungal_phyla, fungal_genera), collapse = "|")
bacteria_keywords <- c("\\bbacterium\\b", "\\bbacteria\\b", "\\bBacterium\\b", 
                      "Alphaproteobacteria", "Acidobacteriaceae", "MAG.*bacterium",
                      "Rhizobium", "Bradyrhizobium", "Silvibacterium", "Mycoplasmopsis",
                      "Sorangium", "Enterobacteriaceae")
bacteria_pattern <- paste(bacteria_keywords, collapse = "|")
plant_keywords <- c("\\bplant\\b", "\\bPlant\\b", "Arabis", "Oryza", "Cannabis", 
                    "Triticum", "Phyllostachys", "Rhopalocnemis")
plant_pattern <- paste(plant_keywords, collapse = "|")
animal_keywords <- c("Carcinus", "\\bcrab\\b", "Lotoria", "Tritonia", "Drosophila", 
                     "Meledella", "Melanogrammus", "Strongylocentrotus", "Nerita",
                     "Ctenopharyngodon", "Fringilla", "Ciona", "Apis", "Euprymna",
                     "\\bsnail\\b", "\\banimal\\b")
animal_pattern <- paste(animal_keywords, collapse = "|")
human_pattern <- "\\bHomo sapiens\\b|\\bhuman\\b|\\bHomo_sapiens\\b"

# Database column mappings
db_cols <- list(
    gtdb = list(name = "gtdb_name", rank = "gtdb_rank"),
    smdb = list(name = "smdb_name", rank = "smdb_rank"),
    pluspf = list(name = "pluspf_name", rank = "pluspf_rank")
)

# BLAST file configuration
blast_files <- list(
    homo_sapiens = file.path(output_dir, paste0(sampleID, "_homo_sapiens_R1_blast.txt")),
    fungal_bacteria = file.path(output_dir, paste0(sampleID, "_smdb_fungal_gtdb_bacteria_R1_blast.txt")),
    smdb_confident_fungi = file.path(output_dir, paste0(sampleID, "_smdb_confident_fungi_R1_blast.txt")),
    unclassified = file.path(output_dir, paste0(sampleID, "_unclassified_R1_blast.txt"))
)

blast_headers <- c("read_id", "accession", "title", "evalue", "bitscore", "identity", 
                   "align_length", "query_start", "query_end", "subject_start", "subject_end")

# Function to extract organism name from BLAST title (vectorized)
extract_organism_name_from_title <- function(title) {
    # Remove NCBI boilerplate patterns
    clean_title <- str_remove(title, "^gi\\|[0-9]+\\|[a-z]+\\|[A-Z0-9_.]+\\|\\s*")
    clean_title <- str_remove(clean_title, "^(PREDICTED:|[A-Z]{1,2}_[0-9]+\\.[0-9]+ |[A-Z0-9]+ )")
    
    # Extract binomial (Genus species)
    binomial <- str_extract(clean_title, "([A-Z][a-z]+ [a-z]+)")
    
    # If no binomial, try after "MAG:" or "uncultured" patterns
    mag_match <- str_extract(clean_title, "(MAG:|uncultured|genome assembly, chromosome|strain|isolate)\\s+([A-Z][a-z]+ [a-z]+)")
    mag_binomial <- if_else(!is.na(mag_match), str_extract(mag_match, "([A-Z][a-z]+ [a-z]+)"), NA_character_)
    binomial <- if_else(!is.na(binomial), binomial, mag_binomial)
    
    # If still no binomial, try genus
    genus_start <- str_extract(clean_title, "^([A-Z][a-z]+)")
    genus_after <- str_extract(clean_title, "(?<=[:\\s]|^)([A-Z][a-z]+)(?=\\s)")
    genus <- if_else(!is.na(genus_start), genus_start, genus_after)
    binomial <- if_else(!is.na(binomial), binomial, genus)
    
    return(binomial)
}

# Expanded list of common fungal genera (partial list - many more exist)
common_fungal_genera <- c(
    "Saccharomyces", "Candida", "Aspergillus", "Penicillium", "Fusarium", 
    "Trichoderma", "Neurospora", "Schizosaccharomyces", "Cryptococcus", 
    "Ustilago", "Magnaporthe", "Botrytis", "Alternaria", "Colletotrichum",
    "Verticillium", "Phytophthora", "Pythium", "Rhizopus", "Mucor",
    "Mortierella", "Glomus", "Agaricus", "Coprinus", "Pleurotus",
    "Lentinula", "Ganoderma", "Trametes", "Phanerochaete", "Schizophyllum",
    "Armillaria", "Heterobasidion", "Serpula", "Coniophora", "Postia",
    "Fomitopsis", "Piptoporus", "Laetiporus", "Grifola", "Meripilus",
    "Rhizophagus", "Xylographa", "Orbilia", "Arthroderma", "Elaphomyces",
    "Hygrocybe", "Anguillospora", "Lipomyces", "Leotia", "Glarea", "Thuemenidium",
    "Puccinia", "Uromyces", "Melampsora", "Cronartium", "Gymnosporangium",
    "Blumeria", "Erysiphe", "Podosphaera", "Uncinula", "Microsphaera",
    "Cladosporium", "Epicoccum", "Stemphylium", "Drechslera", "Bipolaris",
    "Cochliobolus", "Pyrenophora", "Leptosphaeria", "Phaeosphaeria", "Mycosphaerella",
    "Venturia", "Guignardia", "Diaporthe", "Phomopsis", "Gnomonia",
    "Sclerotinia", "Monilinia", "Botryotinia", "Dumontinia", "Ciboria",
    "Geotrichum", "Trichosporon", "Malassezia", "Pityrosporum", "Exophiala",
    "Sporothrix", "Ophiostoma", "Ceratocystis", "Grosmannia", "Leptographium",
    "Beauveria", "Metarhizium", "Cordyceps", "Isaria", "Lecanicillium",
    "Acremonium", "Fusarium", "Gibberella", "Nectria", "Hypocrea",
    "Trichoderma", "Chaetomium", "Sordaria", "Podospora", "Neurospora"
)

# Create pattern for common fungal genera
fungal_genera_pattern <- paste0("\\b(", paste(common_fungal_genera, collapse = "|"), ")\\b")

# Functions
categorize_blast_hit <- function(title, smdb_is_fungi = NA) {
    # First check for explicit keywords
    result <- case_when(
        grepl(human_pattern, title, ignore.case = TRUE) ~ "Human",
        grepl(fungal_keywords, title, ignore.case = TRUE) ~ "Fungi",
        grepl(bacteria_pattern, title, ignore.case = TRUE) ~ "Bacteria",
        grepl(plant_pattern, title, ignore.case = TRUE) ~ "Plant",
        grepl(animal_pattern, title, ignore.case = TRUE) ~ "Animal",
        TRUE ~ NA_character_
    )
    
    # If not categorized yet, try extracting organism name and checking against common fungal genera
    org_name <- extract_organism_name_from_title(title)
    result <- if_else(
        is.na(result) & !is.na(org_name) & grepl(fungal_genera_pattern, org_name, ignore.case = TRUE),
        "Fungi",
        result
    )
    
    # For reads that SMDB classified as fungi, use more lenient categorization
    # If still not categorized and SMDB says it's a fungus, and it doesn't clearly match
    # bacteria/plant/animal/human, treat as potentially fungal
    is_not_other_kingdom <- !grepl(bacteria_pattern, title, ignore.case = TRUE) &
                            !grepl(plant_pattern, title, ignore.case = TRUE) &
                            !grepl(animal_pattern, title, ignore.case = TRUE) &
                            !grepl(human_pattern, title, ignore.case = TRUE)
    
    has_binomial <- !is.na(org_name) & grepl("^[A-Z][a-z]+ [a-z]+", org_name)
    
    result <- if_else(
        is.na(result) & 
        !is.na(smdb_is_fungi) & 
        smdb_is_fungi & 
        is_not_other_kingdom & 
        has_binomial,
        "Fungi",
        result
    )
    
    # Default to Other/Unknown if still not categorized
    result <- if_else(is.na(result), "Other/Unknown", result)
    
    return(result)
}

check_db_correct <- function(blast_category, db_name, db_rank) {
    # If database didn't classify this read (name is NA), return NA (not FALSE)
    # This distinguishes "not classified" from "classified but wrong"
    case_when(
        is.na(db_name) ~ NA,
        blast_category == "Fungi" & (
            grepl(fungal_pattern, db_name, perl = TRUE) |
            (db_rank == "k" & grepl("Fungi", db_name, ignore.case = TRUE)) |
            (db_rank == "p" & grepl(fungal_pattern, db_name, perl = TRUE))
        ) ~ TRUE,
        blast_category == "Bacteria" & (
            grepl("Bacteria|bacteria", db_name) |
            (db_rank == "k" & grepl("Bacteria", db_name, ignore.case = TRUE))
        ) ~ TRUE,
        blast_category == "Human" & grepl("Homo sapiens|Homo_sapiens", db_name, ignore.case = TRUE) ~ TRUE,
        blast_category == "Plant" & grepl("plant|Plant", db_name, ignore.case = TRUE) ~ TRUE,
        blast_category == "Animal" & grepl("animal|Animal", db_name, ignore.case = TRUE) ~ TRUE,
        TRUE ~ FALSE
    )
}

# Load comparison data
if(!file.exists(comparison_file)) {
    stop("Comparison file not found: ", comparison_file, "\nPlease run script 15 first.")
}
comparison <- read_csv(comparison_file, show_col_types = FALSE)

# Process BLAST files - read all results including no_hits
all_blast_data <- map_df(names(blast_files), function(conflict_name) {
    blast_file <- blast_files[[conflict_name]]
    if(!file.exists(blast_file)) return(NULL)
    
    read_tsv(blast_file, show_col_types = FALSE, col_names = blast_headers, skip = 1) %>%
        mutate(
            conflict_type = conflict_name,
            result_type = case_when(
                accession %in% c("no_hits", "blast_failed") ~ accession,
                is.na(accession) ~ "missing",
                TRUE ~ "success"
            )
        ) %>%
        # Convert numeric columns - handle character values (e.g., "no_hits" rows have empty strings)
        mutate(across(any_of(c("bitscore", "identity", "align_length", "evalue", "query_start", "query_end", "subject_start", "subject_end")), 
                      ~ suppressWarnings(as.numeric(.x))))
}, .id = NULL)

# Deduplicate: keep best result per read_id (success > no_hits > blast_failed)
# If multiple successes, keep the one with best bitscore
all_blast_data <- all_blast_data %>%
    mutate(
        result_priority = case_when(
            result_type == "success" ~ 1,
            result_type == "no_hits" ~ 2,
            result_type == "blast_failed" ~ 3,
            TRUE ~ 4
        ),
        bitscore_rank = ifelse(result_type == "success", 
                             rank(-bitscore, ties.method = "first"), 
                             NA_real_)
    ) %>%
    group_by(read_id) %>%
    arrange(result_priority, bitscore_rank) %>%
    slice(1) %>%
    ungroup() %>%
    select(-result_priority, -bitscore_rank)

# Separate successful hits from no_hits/failed
# First join with comparison to get smdb_is_fungi for better categorization
all_results <- all_blast_data %>%
    filter(result_type == "success") %>%
    left_join(comparison, by = "read_id") %>%
    mutate(
        blast_category = categorize_blast_hit(title, smdb_is_fungi),
        reliable_hit = (align_length >= 60 & bitscore >= 50) | (bitscore >= 100),
        unreliable_hit = !reliable_hit,
        gtdb_correct = check_db_correct(blast_category, gtdb_name, gtdb_rank),
        smdb_correct = check_db_correct(blast_category, smdb_name, smdb_rank),
        pluspf_correct = check_db_correct(blast_category, pluspf_name, pluspf_rank)
    )
    
# Track no_hits separately for analysis
no_hits_data <- all_blast_data %>%
    filter(result_type %in% c("no_hits", "blast_failed")) %>%
    left_join(comparison, by = "read_id")

if(nrow(all_blast_data) == 0) {
    stop("No BLAST results found.")
}

# Save detailed comparison (successful hits only) - cleaned version for manuscript
    output_file <- file.path(output_dir, paste0(sampleID, "_blast_vs_databases_comparison.csv"))
if(nrow(all_results) > 0) {
    # Select essential columns only (remove unnecessary ones)
    all_results_cleaned <- all_results %>%
        select(
            # BLAST information
            read_id, accession, title, evalue, bitscore, identity, align_length,
            # Conflict and categorization
            conflict_type, blast_category, reliable_hit,
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
    
    write_csv(all_results_cleaned, output_file)
    cat("Saved cleaned BLAST vs databases comparison to:", output_file, "\n")
    cat("Columns:", ncol(all_results_cleaned), "\n")
}
    
# Generate summary report
    summary_lines <- c(
    "BLAST vs Database Assignments Comparison",
        paste0("Sample: ", sampleID),
        paste0("Date: ", Sys.Date()),
    ""
)

# Overall BLAST statistics
total_blast <- nrow(all_blast_data)
successful_blast <- nrow(all_results)
no_hits_count <- sum(no_hits_data$result_type == "no_hits", na.rm = TRUE)
blast_failed_count <- sum(no_hits_data$result_type == "blast_failed", na.rm = TRUE)

summary_lines <- c(summary_lines,
    "Overall BLAST Statistics:",
    paste0("  Total sequences BLASTed: ", total_blast),
    paste0("  Successful hits: ", successful_blast, 
          " (", round(successful_blast / total_blast * 100, 1), "%)"),
    paste0("  No hits: ", no_hits_count,
          " (", round(no_hits_count / total_blast * 100, 1), "%)"),
    paste0("  BLAST failed: ", blast_failed_count,
          " (", round(blast_failed_count / total_blast * 100, 1), "%)"),
    "",
    "Note: 'No hits' sequences are those that returned no significant matches",
    "in NCBI nt even with lenient expect threshold (1e-4). These may represent:",
    "  - Novel/unknown organisms not in NCBI nt",
    "  - Sequencing artifacts or contamination",
    "  - Sequences with insufficient similarity to database entries",
        ""
    )
    
    # Summary by conflict type (including no_hits)
no_hits_by_conflict <- no_hits_data %>%
    group_by(conflict_type, result_type) %>%
    summarize(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = result_type, values_from = n, values_fill = 0)

all_blast_by_conflict <- all_blast_data %>%
    group_by(conflict_type) %>%
    summarize(
        total_blast = n(),
        successful = sum(result_type == "success", na.rm = TRUE),
        no_hits = sum(result_type == "no_hits", na.rm = TRUE),
        blast_failed = sum(result_type == "blast_failed", na.rm = TRUE),
        .groups = "drop"
    )

# Summary by conflict type (successful hits only)
conflict_summaries <- all_results %>%
    group_by(conflict_type) %>%
    summarize(
        total_reads = n(),
        gtdb_correct = sum(gtdb_correct, na.rm = TRUE),
        smdb_correct = sum(smdb_correct, na.rm = TRUE),
        pluspf_correct = sum(pluspf_correct, na.rm = TRUE),
        none_correct = sum(!gtdb_correct & !smdb_correct & !pluspf_correct, na.rm = TRUE),
        multiple_correct = sum((gtdb_correct + smdb_correct + pluspf_correct) > 1, na.rm = TRUE),
        .groups = "drop"
    )

for(i in 1:nrow(conflict_summaries)) {
    row <- conflict_summaries[i,]
    conflict_data <- all_results %>% filter(conflict_type == row$conflict_type)
    conflict_blast <- all_blast_by_conflict %>% filter(conflict_type == row$conflict_type)
    
    summary_lines <- c(summary_lines,
        paste0("Conflict Type: ", row$conflict_type),
        paste0("  Total sequences BLASTed: ", conflict_blast$total_blast),
        paste0("  Successful hits: ", conflict_blast$successful,
              " (", round(conflict_blast$successful / conflict_blast$total_blast * 100, 1), "%)"),
        paste0("  No hits: ", conflict_blast$no_hits,
              " (", round(conflict_blast$no_hits / conflict_blast$total_blast * 100, 1), "%)"),
        paste0("  BLAST failed: ", conflict_blast$blast_failed,
              " (", round(conflict_blast$blast_failed / conflict_blast$total_blast * 100, 1), "%)"),
        "",
        paste0("  Valid reads (with successful hits): ", row$total_reads),
            "",
            "BLAST Category Distribution:",
            paste0("  ", names(table(conflict_data$blast_category)), ": ", 
                  table(conflict_data$blast_category), collapse = "\n"),
            "",
            "Database Correctness:",
        paste0("  GTDB correct: ", row$gtdb_correct, 
              " (", round(row$gtdb_correct / row$total_reads * 100, 1), "%)"),
        paste0("  SMD correct: ", row$smdb_correct, 
              " (", round(row$smdb_correct / row$total_reads * 100, 1), "%)"),
        paste0("  PlusPF correct: ", row$pluspf_correct, 
              " (", round(row$pluspf_correct / row$total_reads * 100, 1), "%)"),
        paste0("  None correct: ", row$none_correct,
              " (", round(row$none_correct / row$total_reads * 100, 1), "%)"),
        paste0("  Multiple correct: ", row$multiple_correct,
              " (", round(row$multiple_correct / row$total_reads * 100, 1), "%)"),
            ""
        )
        
        # Breakdown by BLAST category
    cat_summaries <- conflict_data %>%
        group_by(blast_category) %>%
        summarize(
            n = n(),
            gtdb_correct = sum(gtdb_correct, na.rm = TRUE),
            smdb_correct = sum(smdb_correct, na.rm = TRUE),
            pluspf_correct = sum(pluspf_correct, na.rm = TRUE),
            none_correct = sum(!gtdb_correct & !smdb_correct & !pluspf_correct, na.rm = TRUE),
            .groups = "drop"
        )
    
    summary_lines <- c(summary_lines, "Detailed Breakdown by BLAST Category:", "")
    for(j in 1:nrow(cat_summaries)) {
        cat_row <- cat_summaries[j,]
        summary_lines <- c(summary_lines,
            paste0("BLAST Category: ", cat_row$blast_category, " (n = ", cat_row$n, ")"),
            paste0("  GTDB correct: ", cat_row$gtdb_correct),
            paste0("  SMD correct: ", cat_row$smdb_correct),
            paste0("  PlusPF correct: ", cat_row$pluspf_correct),
            paste0("  None correct: ", cat_row$none_correct),
            ""
        )
    }
}

# PlusPF "Homo sapiens" false positives analysis
    if("homo_sapiens" %in% all_results$conflict_type) {
        homo_sapiens_data <- all_results %>%
        filter(conflict_type == "homo_sapiens") %>%
            mutate(
                smdb_classified = !is.na(smdb_taxid),
                smdb_confident = ifelse(!is.na(smdb_consistency),
                                       smdb_consistency >= 0.9 & smdb_multiplicity <= 2 & smdb_entropy <= 0.1,
                                       FALSE),
                smdb_filtered_out = !smdb_classified,
                smdb_passed_filter = smdb_classified & smdb_confident,
                smdb_failed_filter = smdb_classified & !smdb_confident
            )
        
        none_correct_data <- homo_sapiens_data %>%
            filter(!gtdb_correct & !smdb_correct & !pluspf_correct)
        
    # Extract organism names for detailed analysis
    organism_extraction_patterns <- list(
        "Carcinus aestuarii" = "Carcinus aestuarii (crab)",
        "Lotoria lotoria" = "Lotoria lotoria (snail)",
        "Nerita plicata" = "Nerita plicata (snail)",
        "Tritonia hombergii" = "Tritonia hombergii",
        "Drosophila melanogaster" = "Drosophila melanogaster (fly)",
        "Melanogrammus aeglefinus" = "Melanogrammus aeglefinus (fish)",
        "Strongylocentrotus droebachiensis" = "Strongylocentrotus droebachiensis (urchin)",
        "Fringilla coelebs" = "Fringilla coelebs (bird)",
        "Ciona intestinalis" = "Ciona intestinalis (tunicate)",
        "Euprymna scolopes" = "Euprymna scolopes (squid)",
        "Rhizophagus irregularis" = "Rhizophagus irregularis (fungus)",
        "Silvibacterium" = "Silvibacterium sp.",
        "Enterobacteriaceae" = "Enterobacteriaceae bacterium",
        "Mycoplasmopsis bovis" = "Mycoplasmopsis bovis",
        "Sorangium" = "Sorangium sp.",
        "Phyllostachys edulis" = "Phyllostachys edulis (bamboo)",
        "Rhopalocnemis phalloides" = "Rhopalocnemis phalloides (plant)",
        "Cannabis sativa" = "Cannabis sativa",
        "Triticum aestivum" = "Triticum aestivum (wheat)",
        "Apis mellifera" = "Apis mellifera (bee)",
        "Meleagris gallopavo" = "Meleagris gallopavo (turkey)",
        "Ctenopharyngodon idella" = "Ctenopharyngodon idella (fish)",
        "Meledella werneri" = "Meledella werneri",
        "Faxonius propinquus" = "Faxonius propinquus (crayfish)",
        "Pararobbsia alpina" = "Pararobbsia alpina"
    )
    
    extract_organism_name <- function(title) {
        # Try specific patterns first
        for(pattern in names(organism_extraction_patterns)) {
            if(grepl(pattern, title, ignore.case = TRUE)) {
                return(organism_extraction_patterns[[pattern]])
            }
        }
        # Try binomial extraction
        binomial <- str_extract(title, "^([A-Z][a-z]+ [a-z]+)")
        if(!is.na(binomial)) return(binomial)
        # Try genus
        genus <- str_extract(title, "^([A-Z][a-z]+)")
        if(!is.na(genus)) return(genus)
        return("Other/Unknown")
    }
    
    none_correct_with_orgs <- none_correct_data %>%
        mutate(organism_name = map_chr(title, extract_organism_name))
    
    # Generate detailed analysis
    detailed_lines <- c(
        "PlusPF 'Homo sapiens' False Positives: Detailed Analysis",
        paste0("Sample: ", sampleID),
        paste0("Date: ", Sys.Date()),
        "",
        "SMD Filter Status:",
        paste0("Total false-positive reads: ", nrow(none_correct_data)),
        paste0("SMD filtered out (not classified): ", sum(none_correct_data$smdb_filtered_out, na.rm = TRUE),
              " (", round(sum(none_correct_data$smdb_filtered_out, na.rm = TRUE) / nrow(none_correct_data) * 100, 1), "%)"),
        paste0("SMD classified but failed filter: ", sum(none_correct_data$smdb_failed_filter, na.rm = TRUE),
              " (", round(sum(none_correct_data$smdb_failed_filter, na.rm = TRUE) / nrow(none_correct_data) * 100, 1), "%)"),
        paste0("SMD classified and passed filter: ", sum(none_correct_data$smdb_passed_filter, na.rm = TRUE),
              " (", round(sum(none_correct_data$smdb_passed_filter, na.rm = TRUE) / nrow(none_correct_data) * 100, 1), "%)"),
        "",
        "BLAST Quality Metrics:",
        paste0("Mean identity: ", round(mean(none_correct_data$identity, na.rm = TRUE), 2), "%"),
        paste0("Median identity: ", round(median(none_correct_data$identity, na.rm = TRUE), 2), "%"),
        paste0("Mean bitscore: ", round(mean(none_correct_data$bitscore, na.rm = TRUE), 2)),
        paste0("Mean E-value: ", format(mean(none_correct_data$evalue, na.rm = TRUE), scientific = TRUE)),
        "",
        "BLAST Category Distribution:",
        paste0(names(table(none_correct_data$blast_category)), ": ", 
              table(none_correct_data$blast_category), collapse = "\n"),
        "",
        "Top BLAST Hit Organisms:",
        paste0(none_correct_with_orgs %>%
               count(organism_name, blast_category, sort = TRUE) %>%
               head(20) %>%
               mutate(line = paste0(organism_name, " [", blast_category, "]: ", n)) %>%
               pull(line), collapse = "\n")
    )
    
    # PlusPF quality metrics
        pluspf_quality <- none_correct_data %>%
            filter(!is.na(pluspf_consistency)) %>%
            summarize(
                mean_consistency = mean(as.numeric(pluspf_consistency), na.rm = TRUE),
                mean_multiplicity = mean(as.numeric(pluspf_multiplicity), na.rm = TRUE),
                mean_entropy = mean(as.numeric(pluspf_entropy), na.rm = TRUE),
                n_confident = sum(pluspf_is_confident, na.rm = TRUE),
                pct_confident = n_confident / n() * 100
            )
    
    detailed_lines <- c(detailed_lines,
        "",
        "PlusPF Quality Metrics:",
            paste0("Mean consistency: ", round(pluspf_quality$mean_consistency, 3)),
            paste0("Mean multiplicity: ", round(pluspf_quality$mean_multiplicity, 3)),
            paste0("Mean entropy: ", round(pluspf_quality$mean_entropy, 3)),
            paste0("Confident classifications: ", pluspf_quality$n_confident,
              " (", round(pluspf_quality$pct_confident, 1), "%)")
    )
    
    # SMD/GTDB classification of PlusPF confident reads
        pluspf_confident_human <- homo_sapiens_data %>%
            filter(pluspf_is_confident == TRUE) %>%
            mutate(
                smdb_classified = !is.na(smdb_taxid),
                gtdb_classified = !is.na(gtdb_taxid),
                smdb_confident = ifelse(!is.na(smdb_consistency),
                                       smdb_consistency >= 0.9 & smdb_multiplicity <= 2 & smdb_entropy <= 0.1,
                                       FALSE),
                gtdb_confident = ifelse(!is.na(gtdb_consistency),
                                      gtdb_consistency >= 0.9 & gtdb_multiplicity <= 2 & gtdb_entropy <= 0.1,
                                      FALSE)
            )
        
        smd_status <- pluspf_confident_human %>%
            summarize(
                total = n(),
                not_classified = sum(!smdb_classified, na.rm = TRUE),
                classified_not_confident = sum(smdb_classified & !smdb_confident, na.rm = TRUE),
                classified_confident = sum(smdb_classified & smdb_confident, na.rm = TRUE)
            )
    
        gtdb_status <- pluspf_confident_human %>%
            summarize(
                total = n(),
                not_classified = sum(!gtdb_classified, na.rm = TRUE),
                classified_not_confident = sum(gtdb_classified & !gtdb_confident, na.rm = TRUE),
                classified_confident = sum(gtdb_classified & gtdb_confident, na.rm = TRUE)
            )
    
        smd_assignments <- pluspf_confident_human %>%
            filter(smdb_classified) %>%
            count(smdb_rank, smdb_name, sort = TRUE) %>%
            head(15)
        
        gtdb_assignments <- pluspf_confident_human %>%
            filter(gtdb_classified) %>%
            count(gtdb_rank, gtdb_name, sort = TRUE) %>%
            head(15)
    
    detailed_lines <- c(detailed_lines,
        "",
        "How SMD and GTDB Classified PlusPF Confident 'Homo sapiens' Reads:",
            paste0("Total PlusPF confident 'Homo sapiens' reads: ", nrow(pluspf_confident_human)),
            "",
            "SMD Classification Status:",
            paste0("  Not classified: ", smd_status$not_classified,
                  " (", round(smd_status$not_classified / smd_status$total * 100, 1), "%)"),
            paste0("  Classified but not confident: ", smd_status$classified_not_confident,
                  " (", round(smd_status$classified_not_confident / smd_status$total * 100, 1), "%)"),
            paste0("  Classified and confident: ", smd_status$classified_confident,
                  " (", round(smd_status$classified_confident / smd_status$total * 100, 1), "%)"),
            "",
            "GTDB Classification Status:",
            paste0("  Not classified: ", gtdb_status$not_classified,
                  " (", round(gtdb_status$not_classified / gtdb_status$total * 100, 1), "%)"),
            paste0("  Classified but not confident: ", gtdb_status$classified_not_confident,
                  " (", round(gtdb_status$classified_not_confident / gtdb_status$total * 100, 1), "%)"),
            paste0("  Classified and confident: ", gtdb_status$classified_confident,
                  " (", round(gtdb_status$classified_confident / gtdb_status$total * 100, 1), "%)"),
            "",
            "Top SMD Assignments:",
        if(nrow(smd_assignments) > 0) {
            paste0("  ", smd_assignments$smdb_rank, " | ", smd_assignments$smdb_name, ": ", 
                  smd_assignments$n, collapse = "\n")
        } else "  None",
            "",
            "Top GTDB Assignments:",
        if(nrow(gtdb_assignments) > 0) {
            paste0("  ", gtdb_assignments$gtdb_rank, " | ", gtdb_assignments$gtdb_name, ": ", 
                  gtdb_assignments$n, collapse = "\n")
        } else "  None"
    )
    
    detailed_file <- file.path(output_dir, paste0(sampleID, "_homo_sapiens_false_positives_analysis.txt"))
    writeLines(detailed_lines, detailed_file)
}

# Add no_hits analysis section
if(nrow(no_hits_data) > 0) {
    summary_lines <- c(summary_lines,
        "",
        "=== No_Hits Analysis ===",
        "",
        "Sequences that returned no significant BLAST hits:",
        ""
    )
    
    # No_hits by conflict type
    no_hits_summary <- no_hits_data %>%
        group_by(conflict_type) %>%
        summarize(
            total = n(),
            no_hits = sum(result_type == "no_hits", na.rm = TRUE),
            blast_failed = sum(result_type == "blast_failed", na.rm = TRUE),
            .groups = "drop"
        )
    
    for(i in 1:nrow(no_hits_summary)) {
        row <- no_hits_summary[i,]
        summary_lines <- c(summary_lines,
            paste0("Conflict Type: ", row$conflict_type),
            paste0("  Total no_hits/failed: ", row$total),
            paste0("  No hits: ", row$no_hits,
                  " (", round(row$no_hits / row$total * 100, 1), "%)"),
            paste0("  BLAST failed: ", row$blast_failed,
                  " (", round(row$blast_failed / row$total * 100, 1), "%)"),
            ""
        )
    }
    
    # Database classification status of no_hits sequences
    no_hits_with_db <- no_hits_data %>%
        mutate(
            gtdb_classified = !is.na(gtdb_taxid),
            smdb_classified = !is.na(smdb_taxid),
            pluspf_classified = !is.na(pluspf_taxid),
            any_classified = gtdb_classified | smdb_classified | pluspf_classified
        ) %>%
        group_by(conflict_type) %>%
        summarize(
            total = n(),
            gtdb_classified = sum(gtdb_classified, na.rm = TRUE),
            smdb_classified = sum(smdb_classified, na.rm = TRUE),
            pluspf_classified = sum(pluspf_classified, na.rm = TRUE),
            any_classified = sum(any_classified, na.rm = TRUE),
            none_classified = sum(!any_classified, na.rm = TRUE),
            .groups = "drop"
        )
    
    summary_lines <- c(summary_lines,
        "Database Classification Status of No_Hits Sequences:",
        ""
    )
    
    for(i in 1:nrow(no_hits_with_db)) {
        row <- no_hits_with_db[i,]
        summary_lines <- c(summary_lines,
            paste0("Conflict Type: ", row$conflict_type),
            paste0("  Total: ", row$total),
            paste0("  GTDB classified: ", row$gtdb_classified,
                  " (", round(row$gtdb_classified / row$total * 100, 1), "%)"),
            paste0("  SMD classified: ", row$smdb_classified,
                  " (", round(row$smdb_classified / row$total * 100, 1), "%)"),
            paste0("  PlusPF classified: ", row$pluspf_classified,
                  " (", round(row$pluspf_classified / row$total * 100, 1), "%)"),
            paste0("  Any database classified: ", row$any_classified,
                  " (", round(row$any_classified / row$total * 100, 1), "%)"),
            paste0("  None classified: ", row$none_classified,
                  " (", round(row$none_classified / row$total * 100, 1), "%)"),
            ""
        )
    }
    
    summary_lines <- c(summary_lines,
        "Interpretation:",
        "  - No_hits sequences cannot be verified via BLAST, so database correctness",
        "    cannot be assessed for these sequences.",
        "  - If databases classified these sequences, they may be:",
        "    * Correct assignments to novel organisms not in NCBI nt",
        "    * Incorrect assignments that cannot be verified",
        "    * Sequences from organisms with poor NCBI nt coverage",
        "  - If no database classified these sequences, they represent:",
        "    * Truly unclassifiable sequences",
        "    * Sequencing artifacts or contamination",
        "    * Novel organisms missing from all databases",
        ""
    )
}

# Write summary
summary_file <- file.path(output_dir, paste0(sampleID, "_blast_vs_databases_summary.txt"))
writeLines(summary_lines, summary_file)
