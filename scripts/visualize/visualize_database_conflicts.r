#!/usr/bin/env Rscript
# Visualize database conflicts and BLAST verification support
# Shows when SMD assignments conflict with other databases and how BLAST verifies SMD
#
# Usage: Rscript scripts/visualize/visualize_database_conflicts.r [sampleID]
#
# Input:  
#   - taxonomic_assignment_comparison_{sampleID}.csv (from script 15)
#   - BLAST summary files (from script 18, optional)
# Output: 
#   - {sampleID}_database_conflicts.png (grouped bar chart showing conflicts and BLAST support)

library(tidyverse)
library(ggplot2)
library(patchwork)

# Sample to analyze (default to ORNL_046-O-20170621-COMP)
sampleID <- if(length(commandArgs(trailingOnly = TRUE)) > 0) {
    commandArgs(trailingOnly = TRUE)[1]
} else {
    "ORNL_046-O-20170621-COMP"
}

cat("=== Visualizing Database Conflicts and BLAST Verification ===\n")
cat("Sample:", sampleID, "\n\n")

# File paths
input_file <- file.path("data/classification/analysis_files", 
                        paste0("taxonomic_assignment_comparison_", sampleID, ".csv"))
output_dir <- "manuscript_figures"
blast_summary_dir <- "data/classification/analysis_files"

if(!file.exists(input_file)) {
    stop("Comparison file not found: ", input_file, 
         "\nPlease run script 15 first: Rscript scripts/summarize_outputs/15_compare_taxonomic_assignments.r ", sampleID)
}

cat("Reading comparison data...\n")
comparison <- read_csv(input_file, show_col_types = FALSE)

# Define fungal phyla
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")
fungal_pattern <- paste0("(", paste(fungal_phyla, collapse = "|"), ")")

# Identify confident classifications and kingdoms
comparison <- comparison %>%
    mutate(
        smdb_is_confident = smdb_consistency >= 0.9 & smdb_multiplicity <= 2 & smdb_entropy <= 0.1,
        gtdb_is_confident = gtdb_consistency >= 0.9 & gtdb_multiplicity <= 2 & gtdb_entropy <= 0.1,
        pluspf_is_confident = ifelse(!is.na(pluspf_consistency),
                                    pluspf_consistency >= 0.9 & pluspf_multiplicity <= 2 & pluspf_entropy <= 0.1,
                                    FALSE),
        smdb_is_fungi = case_when(
            smdb_rank == "k" & grepl("^k__Fungi|^Fungi", smdb_name) ~ TRUE,
            smdb_rank == "p" & grepl(fungal_pattern, smdb_name, perl = TRUE) ~ TRUE,
            grepl(paste0("^p__(", paste(fungal_phyla, collapse = "|"), ")"), smdb_name, perl = TRUE) ~ TRUE,
            grepl(paste0("^(", paste(fungal_phyla, collapse = "|"), ")"), smdb_name, perl = TRUE) ~ TRUE,
            TRUE ~ FALSE
        ),
        gtdb_is_bacteria = case_when(
            gtdb_rank == "k" & grepl("^k__Bacteria|^Bacteria|k__Bacteria", gtdb_name) ~ TRUE,
            gtdb_rank == "p" & grepl("Proteobacteria|Actinobacteria|Bacteroidota|Firmicutes", gtdb_name, ignore.case = TRUE) ~ TRUE,
            grepl("^k__Bacteria|^Bacteria", gtdb_name) ~ TRUE,
            TRUE ~ FALSE
        ),
        pluspf_is_human = grepl("Homo sapiens|Homo_sapiens", pluspf_name, ignore.case = TRUE)
    )

# Read BLAST verification results from actual BLAST files
blast_support <- list(
    homo_sapiens = NULL,
    fungal_bacteria = NULL
)

# Read BLAST results directly from files
homo_sapiens_blast_file <- file.path(blast_summary_dir, paste0(sampleID, "_homo_sapiens_R1_blast.txt"))
fungal_bacteria_blast_file <- file.path(blast_summary_dir, paste0(sampleID, "_smdb_fungal_gtdb_bacteria_R1_blast.txt"))

# Process PlusPF "Homo sapiens" BLAST results
# Prefer files with 100 reads (check file size to determine)
homo_sapiens_blast_files <- c(
    file.path(blast_summary_dir, paste0(sampleID, "_homo_sapiens_R1_blast_100.txt")),
    file.path(blast_summary_dir, paste0(sampleID, "_homo_sapiens_R1_blast.txt")),
    file.path(blast_summary_dir, paste0(sampleID, "_homo_sapiens_R1_blast_sample.txt"))
)
# Prefer files that exist and have reasonable size (indicating ~100 reads)
existing_files <- homo_sapiens_blast_files[file.exists(homo_sapiens_blast_files)]
if(length(existing_files) > 0) {
    file_sizes <- sapply(existing_files, function(f) file.info(f)$size)
    # Prefer larger files (more likely to have 100 reads)
    homo_sapiens_blast_file <- existing_files[which.max(file_sizes)]
} else {
    homo_sapiens_blast_file <- NA
}

if(!is.na(homo_sapiens_blast_file) && file.exists(homo_sapiens_blast_file)) {
    cat("Reading BLAST results for PlusPF 'Homo sapiens' conflicts from:", basename(homo_sapiens_blast_file), "\n")
    blast_results <- read_tsv(homo_sapiens_blast_file, show_col_types = FALSE,
                             col_names = c("read_id", "accession", "title", "evalue", 
                                         "bitscore", "identity", "align_length", 
                                         "query_start", "query_end", "subject_start", "subject_end"))
    
    # Filter valid results
    valid_results <- blast_results %>%
        filter(!accession %in% c("no_hits", "blast_failed") & !is.na(accession))
    
    if(nrow(valid_results) > 0) {
        # Categorize BLAST hits
        valid_results <- valid_results %>%
            mutate(
                is_fungi = grepl("fungi|Fungi|Ascomycota|Basidiomycota|Mucoromycota|fungus|Elaphomyces|Hygrocybe|Rhizophagus|Anguillospora|Lipomyces|Leotia|Glarea|Thuemenidium|Xylographa|Orbilia|Arthroderma", title, ignore.case = TRUE),
                is_human = grepl("Homo sapiens|human", title, ignore.case = TRUE),
                is_bacteria = grepl("bacterium|bacteria|Bacterium|Alphaproteobacteria|Acidobacteriaceae", title, ignore.case = TRUE),
                is_plant = grepl("plant|Plant|Arabis|Oryza|Cannabis", title, ignore.case = TRUE),
                is_animal = grepl("Carcinus|crab|Eulemur|Luscinia|Melanogrammus|snail", title, ignore.case = TRUE)
            )
        
        fungi_count <- sum(valid_results$is_fungi, na.rm = TRUE)
        human_count <- sum(valid_results$is_human, na.rm = TRUE)
        total_blast <- nrow(blast_results)
        
        blast_support$homo_sapiens <- list(
            total = total_blast,
            valid = nrow(valid_results),
            fungi = fungi_count,
            human = human_count,
            bacteria = sum(valid_results$is_bacteria, na.rm = TRUE),
            plant = sum(valid_results$is_plant, na.rm = TRUE),
            animal = sum(valid_results$is_animal, na.rm = TRUE),
            pct_fungi = ifelse(nrow(valid_results) > 0, fungi_count / nrow(valid_results) * 100, 0)
        )
        cat("  Total BLASTed:", total_blast, "| Valid hits:", nrow(valid_results), 
            "| Fungi:", fungi_count, "| Human:", human_count, "\n")
    }
}

# Process SMD fungal → GTDB bacteria BLAST results
if(file.exists(fungal_bacteria_blast_file)) {
    cat("Reading BLAST results for SMD fungal → GTDB bacteria conflicts...\n")
    blast_results <- read_tsv(fungal_bacteria_blast_file, show_col_types = FALSE,
                             col_names = c("read_id", "accession", "title", "evalue", 
                                         "bitscore", "identity", "align_length", 
                                         "query_start", "query_end", "subject_start", "subject_end"))
    
    # Filter valid results
    valid_results <- blast_results %>%
        filter(!accession %in% c("no_hits", "blast_failed") & !is.na(accession))
    
    if(nrow(valid_results) > 0) {
        # Categorize BLAST hits
        valid_results <- valid_results %>%
            mutate(
                is_fungi = grepl("fungi|Fungi|Ascomycota|Basidiomycota|Mucoromycota|fungus|Elaphomyces|Hygrocybe|Rhizophagus|Anguillospora|Lipomyces|Leotia|Glarea|Thuemenidium|Xylographa|Orbilia|Arthroderma", title, ignore.case = TRUE),
                is_bacteria = grepl("bacterium|bacteria|Bacterium|Alphaproteobacteria|Acidobacteriaceae|MAG.*bacterium", title, ignore.case = TRUE)
            )
        
        fungi_count <- sum(valid_results$is_fungi, na.rm = TRUE)
        bacteria_count <- sum(valid_results$is_bacteria, na.rm = TRUE)
        total_blast <- nrow(blast_results)
        
        blast_support$fungal_bacteria <- list(
            total = total_blast,
            valid = nrow(valid_results),
            fungi = fungi_count,
            bacteria = bacteria_count,
            pct_fungi = ifelse(nrow(valid_results) > 0, fungi_count / nrow(valid_results) * 100, 0)
        )
        cat("  Total BLASTed:", total_blast, "| Valid hits:", nrow(valid_results), 
            "| Fungi:", fungi_count, "| Bacteria:", bacteria_count, "\n")
    }
}

# Only analyze conflicts that were actually BLASTed
cat("\n=== Identifying BLASTed Conflicts ===\n")

# Get read IDs that were BLASTed
homo_sapiens_read_ids <- NULL
fungal_bacteria_read_ids <- NULL

if(!is.null(blast_support$homo_sapiens) && file.exists(homo_sapiens_blast_file)) {
    blast_results <- read_tsv(homo_sapiens_blast_file, show_col_types = FALSE,
                             col_names = c("read_id", "accession", "title", "evalue", 
                                         "bitscore", "identity", "align_length", 
                                         "query_start", "query_end", "subject_start", "subject_end"))
    homo_sapiens_read_ids <- unique(blast_results$read_id)
    cat("PlusPF 'Homo sapiens' reads BLASTed:", length(homo_sapiens_read_ids), "\n")
}

if(!is.null(blast_support$fungal_bacteria) && file.exists(fungal_bacteria_blast_file)) {
    blast_results <- read_tsv(fungal_bacteria_blast_file, show_col_types = FALSE,
                             col_names = c("read_id", "accession", "title", "evalue", 
                                         "bitscore", "identity", "align_length", 
                                         "query_start", "query_end", "subject_start", "subject_end"))
    fungal_bacteria_read_ids <- unique(blast_results$read_id)
    cat("SMD fungal → GTDB bacteria reads BLASTed:", length(fungal_bacteria_read_ids), "\n")
}

# Build visualization data only for BLASTed conflicts
viz_data <- data.frame()

# PlusPF "Homo sapiens" conflicts (only those BLASTed)
if(!is.null(homo_sapiens_read_ids) && !is.null(blast_support$homo_sapiens)) {
    pluspf_human <- comparison %>%
        filter(read_id %in% homo_sapiens_read_ids & !is.na(pluspf_taxid) & pluspf_is_human)
    
    if(nrow(pluspf_human) > 0) {
        # Count SMD fungal assignments in BLASTed reads
        smd_fungi_in_human <- pluspf_human %>%
            filter(smdb_is_fungi) %>%
            nrow()
        
        # Get BLAST verification counts
        blast_verified_fungi <- blast_support$homo_sapiens$fungi
        blast_verified_total <- blast_support$homo_sapiens$valid
        
        # Show all BLASTed conflicts (PlusPF says human, but they're not)
        # Count what SMD assigns and what BLAST verifies
        smd_eukaryote_in_human <- pluspf_human %>%
            filter(grepl("Eukaryota|eukaryota", smdb_name) | 
                   (smdb_rank == "k" & grepl("Eukaryota", smdb_name)) |
                   smdb_is_fungi) %>%
            nrow()
        
        # Show the conflict - PlusPF says human, but BLAST/SMD show otherwise
        viz_data <- bind_rows(
            viz_data,
            data.frame(
                conflict_type = "PlusPF: Homo sapiens\nvs\nSMD: Fungi/Eukaryota",
                total_conflicts = nrow(pluspf_human),
                smd_fungi_count = smd_fungi_in_human,
                blast_verified_fungi = blast_verified_fungi,
                blast_verified_total = blast_verified_total
            )
        )
        cat("PlusPF 'Homo sapiens' reads BLASTed:", nrow(pluspf_human), "\n")
        cat("  SMD assigns as fungi/eukaryotes:", smd_eukaryote_in_human, "\n")
        cat("  BLAST verified as fungi:", blast_verified_fungi, 
            "(of", blast_verified_total, "valid BLAST hits)\n")
    }
}

# SMD fungal → GTDB bacteria conflicts (only those BLASTed)
if(!is.null(fungal_bacteria_read_ids) && !is.null(blast_support$fungal_bacteria)) {
    smdb_confident_fungi <- comparison %>%
        filter(read_id %in% fungal_bacteria_read_ids & 
               smdb_is_confident & smdb_is_fungi & !is.na(smdb_taxid))
    
    fungal_bacteria_mismatch <- smdb_confident_fungi %>%
        filter(!is.na(gtdb_taxid) & gtdb_is_bacteria)
    
    if(nrow(fungal_bacteria_mismatch) > 0) {
        viz_data <- bind_rows(
            viz_data,
            data.frame(
                conflict_type = "SMD: Fungi\nvs\nGTDB: Bacteria",
                total_conflicts = nrow(fungal_bacteria_mismatch),
                smd_fungi_count = nrow(fungal_bacteria_mismatch),
                blast_verified_fungi = blast_support$fungal_bacteria$fungi,
                blast_verified_total = blast_support$fungal_bacteria$valid
            )
        )
    }
}

if(nrow(viz_data) > 0) {
    # Prepare data for visualization with both conflict counts and BLAST verification
    plot_data <- viz_data %>%
        mutate(
            conflict_label = conflict_type,
            conflicts_identified = total_conflicts,
            blast_verified = blast_verified_fungi
        ) %>%
        select(conflict_label, conflicts_identified, blast_verified) %>%
        pivot_longer(cols = c(conflicts_identified, blast_verified), 
                    names_to = "category", values_to = "count") %>%
        mutate(
            category_label = case_when(
                category == "conflicts_identified" ~ "Conflicts identified",
                category == "blast_verified" ~ "BLAST verified (Fungi)"
            )
        )
    
    # Create grouped bar chart
    p1 <- ggplot(plot_data, aes(x = conflict_label, y = count, fill = category_label)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7) +
        geom_text(aes(label = count), 
                 position = position_dodge(width = 0.7), 
                 vjust = -0.3, size = 3.5) +
        scale_fill_manual(
            name = "",
            values = c("Conflicts identified" = "#E31A1C",
                      "BLAST verified (Fungi)" = "#33A02C")
        ) +
        labs(
            x = "Conflict Type",
            y = "Number of Reads",
            title = "Database Conflicts: SMD Assignments vs Other Databases\n(BLAST Verification Supports SMD Fungal Assignments)",
            subtitle = if(any(viz_data$blast_verified_total > 0 & viz_data$blast_verified_total < viz_data$total_conflicts)) {
                paste0("Note: BLAST verification based on sample of ", 
                      sum(viz_data$blast_verified_total), " reads")
            } else {
                ""
            }
        ) +
        theme_bw(base_size = 12) +
        theme(
            legend.position = "bottom",
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 10)
        )
    
    # Save visualization
    output_file <- file.path(output_dir, paste0(sampleID, "_database_conflicts.png"))
    ggsave(output_file, p1, width = 10, height = 6, dpi = 300)
    cat("\nSaved visualization to:", output_file, "\n")
    cat("\nSummary:\n")
    for(i in 1:nrow(viz_data)) {
        cat("  ", viz_data$conflict_type[i], ":\n")
        cat("    Conflicts identified:", viz_data$total_conflicts[i], "\n")
        cat("    BLAST verified as fungi:", viz_data$blast_verified_fungi[i], 
            sprintf("(%.1f%% of BLASTed reads)\n", 
                   viz_data$blast_verified_fungi[i] / viz_data$blast_verified_total[i] * 100))
    }
} else {
    cat("No BLASTed conflicts found to visualize.\n")
    cat("Make sure BLAST results are available by running scripts 16-18.\n")
}

cat("\n=== Visualization complete ===\n")

