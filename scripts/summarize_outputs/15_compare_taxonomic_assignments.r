#!/usr/bin/env Rscript
# 15: Compare taxonomic assignments between GTDB, SoilMicrobeDB, and PlusPF
# Performs multiple cross-database analyses:
#   - Analysis 1: PlusPF "Homo sapiens" misclassifications (saves read IDs for BLAST)
#   - Analysis 2: GTDB confident fungal reads → SMD assignments
#   - Analysis 3: SMD confident fungal reads → GTDB bacteria mismatches (saves read IDs for BLAST)
#   - Analysis 4: SMD confident bacterial reads → GTDB assignments
#   - Analysis 5: Overall taxonomic rank distribution
#
# Usage: Rscript scripts/summarize_outputs/15_compare_taxonomic_assignments.r [sampleID]
#
# Input:  *_scores.output files (from 02_run_architeuthis.sh)
# Output: 
#   - taxonomic_assignment_comparison_{sampleID}.csv (full comparison)
#   - taxonomic_assignment_comparison_{sampleID}_summary.txt (summary report)
#   - {sampleID}_homo_sapiens_read_ids.txt (PlusPF "Homo sapiens" read IDs for BLAST, if found)
#   - {sampleID}_smdb_fungal_gtdb_bacteria_read_ids.txt (SMD fungal → GTDB bacteria read IDs for BLAST, if found)
#
# Workflow for Misclassification Verification:
#   1. Run this script (15) → identifies BOTH misclassification types and generates read ID lists:
#      - {sampleID}_homo_sapiens_read_ids.txt (PlusPF "Homo sapiens" misclassifications)
#      - {sampleID}_smdb_fungal_gtdb_bacteria_read_ids.txt (SMD fungal → GTDB bacteria mismatches)
#   2. Run script 16 → automatically extracts reads from both lists to FASTA format
#      Usage: python3 scripts/summarize_outputs/16_extract_reads_to_fasta.py <sampleID> <R1.fastq.gz> <R2.fastq.gz>
#   3. Run script 17 → automatically BLASTs all FASTA files
#      Usage: python3 scripts/summarize_outputs/17_blast_reads_ncbi.py <sampleID> [sample_size]
#   4. Run script 18 → automatically analyzes all BLAST results
#      Usage: Rscript scripts/summarize_outputs/18_analyze_blast_results.r <sampleID>

library(tidyverse)
library(data.table)

# Sample to analyze (default to ORNL_046-O-20170621-COMP)
sampleID <- if(length(commandArgs(trailingOnly = TRUE)) > 0) {
    commandArgs(trailingOnly = TRUE)[1]
} else {
    "ORNL_046-O-20170621-COMP"
}

cat("=== Comparing Taxonomic Assignments ===\n")
cat("Sample:", sampleID, "\n\n")

# File paths
scores_dir <- "data/classification/02_bracken_output"
gtdb_file <- file.path(scores_dir, paste0(sampleID, "_gtdb_207_unfiltered_scores.output"))
smdb_file <- file.path(scores_dir, paste0(sampleID, "_soil_microbe_db_scores.output"))
pluspf_file <- file.path(scores_dir, paste0(sampleID, "_pluspf_scores.output"))

# Check if files exist
if(!file.exists(gtdb_file)) {
    stop("GTDB file not found: ", gtdb_file)
}
if(!file.exists(smdb_file)) {
    stop("SoilMicrobeDB file not found: ", smdb_file)
}
if(!file.exists(pluspf_file)) {
    cat("Warning: PlusPF file not found: ", pluspf_file, "\n")
    cat("Continuing without PlusPF...\n\n")
    pluspf_file <- NULL
}

cat("Reading scores files...\n")
gtdb_scores <- fread(gtdb_file, showProgress = FALSE)
smdb_scores <- fread(smdb_file, showProgress = FALSE)
if(!is.null(pluspf_file)) {
    pluspf_scores <- fread(pluspf_file, showProgress = FALSE)
} else {
    pluspf_scores <- NULL
}

cat("GTDB reads:", nrow(gtdb_scores), "\n")
cat("SMD reads:", nrow(smdb_scores), "\n")
if(!is.null(pluspf_scores)) {
    cat("PlusPF reads:", nrow(pluspf_scores), "\n")
}
cat("\n")

# Filter to only classified reads (n_kmers > 0)
gtdb_classified <- gtdb_scores %>% filter(n_kmers > 0)
smdb_classified <- smdb_scores %>% filter(n_kmers > 0)
if(!is.null(pluspf_scores)) {
    pluspf_classified <- pluspf_scores %>% filter(n_kmers > 0)
} else {
    pluspf_classified <- NULL
}

cat("GTDB classified reads:", nrow(gtdb_classified), "\n")
cat("SMD classified reads:", nrow(smdb_classified), "\n")
if(!is.null(pluspf_classified)) {
    cat("PlusPF classified reads:", nrow(pluspf_classified), "\n")
}
cat("\n")

# Define fungal phyla
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")
fungal_pattern <- paste0("(", paste(fungal_phyla, collapse = "|"), ")")

# Identify fungi in each database
# Fungi can be identified by:
# 1. rank == "k" and name contains "Fungi" (kingdom level)
# 2. rank == "p" and name contains fungal phylum names (phylum level)
# 3. For other ranks, check if name starts with phylum prefix (e.g., "p__Basidiomycota")
# Note: We avoid matching species names that contain "fungi" but are bacteria

identify_fungi <- function(df) {
    df %>%
        mutate(
            is_fungi = case_when(
                # Kingdom level
                rank == "k" & grepl("^k__Fungi|^Fungi", name) ~ TRUE,
                # Phylum level - match phylum names
                rank == "p" & grepl(fungal_pattern, name, perl = TRUE) ~ TRUE,
                # Other ranks - check if name contains phylum prefix (e.g., "p__Basidiomycota")
                grepl(paste0("^p__(", paste(fungal_phyla, collapse = "|"), ")"), name, perl = TRUE) ~ TRUE,
                # Or if name contains phylum name at start (for non-prefixed formats)
                grepl(paste0("^(", paste(fungal_phyla, collapse = "|"), ")"), name, perl = TRUE) ~ TRUE,
                TRUE ~ FALSE
            )
        )
}

gtdb_classified <- identify_fungi(gtdb_classified)
smdb_classified <- identify_fungi(smdb_classified)
if(!is.null(pluspf_classified)) {
    pluspf_classified <- identify_fungi(pluspf_classified)
}

# Identify bacteria in GTDB
identify_bacteria <- function(df) {
    df %>%
        mutate(
            is_bacteria = case_when(
                # Kingdom level
                rank == "k" & grepl("^k__Bacteria|^Bacteria|k__Bacteria", name) ~ TRUE,
                # Phylum level - common bacterial phyla
                rank == "p" & grepl("Proteobacteria|Actinobacteria|Bacteroidota|Firmicutes|Acidobacteria|Chloroflexota|Planctomycetota|Verrucomicrobiota|Cyanobacteria|Desulfobacterota|Myxococcota|Nitrospirota|Spirochaetota|Thermotogota|Deinococcota|Synergistota|Fusobacteriota|Campilobacterota|Patescibacteria|Elusimicrobiota|Armatimonadota|Chlorobi|Gemmatimonadota|Lentisphaerota|Thermodesulfobacteriota|Dependentiae|Omnitrophota|Fibrobacterota|Deferribacterota|Aquificota|Caldisericota|Coprothermobacterota|Dictyoglomota|Kiritimatiellota|Rhodothermaeota|Sumerlaeota|Tectomicrobia|Thermodesulfobiota|Thermosulfidibacterota|WPS-2|Zixibacteria", name, ignore.case = TRUE) ~ TRUE,
                # Other ranks - check if name contains kingdom prefix
                grepl("^k__Bacteria|^Bacteria", name) ~ TRUE,
                # Check if name contains bacterial phylum prefix
                grepl(paste0("^p__(Proteobacteria|Actinobacteria|Bacteroidota|Firmicutes|Acidobacteria|Chloroflexota|Planctomycetota|Verrucomicrobiota|Cyanobacteria|Desulfobacterota|Myxococcota|Nitrospirota|Spirochaetota|Thermotogota|Deinococcota|Synergistota|Fusobacteriota|Campilobacterota|Patescibacteria|Elusimicrobiota|Armatimonadota|Chlorobi|Gemmatimonadota|Lentisphaerota|Thermodesulfobacteriota|Dependentiae|Omnitrophota|Fibrobacterota|Deferribacterota|Aquificota|Caldisericota|Coprothermobacterota|Dictyoglomota|Kiritimatiellota|Rhodothermaeota|Sumerlaeota|Tectomicrobia|Thermodesulfobiota|Thermosulfidibacterota|WPS-2|Zixibacteria)"), name, perl = TRUE) ~ TRUE,
                TRUE ~ FALSE
            )
        )
}

gtdb_classified <- identify_bacteria(gtdb_classified)

cat("=== Fungal Read Counts ===\n")
cat("GTDB fungal reads:", sum(gtdb_classified$is_fungi), "\n")
cat("SMD fungal reads:", sum(smdb_classified$is_fungi), "\n")
if(!is.null(pluspf_classified)) {
    cat("PlusPF fungal reads:", sum(pluspf_classified$is_fungi), "\n")
}
cat("\n")

# Define "confident" classification
# Using filter thresholds: consistency >= 0.9, multiplicity <= 2, entropy <= 0.1
gtdb_classified <- gtdb_classified %>%
    mutate(
        is_confident = consistency >= 0.9 & multiplicity <= 2 & entropy <= 0.1
    )

smdb_classified <- smdb_classified %>%
    mutate(
        is_confident = consistency >= 0.9 & multiplicity <= 2 & entropy <= 0.1
    )
if(!is.null(pluspf_classified)) {
    pluspf_classified <- pluspf_classified %>%
        mutate(
            is_confident = consistency >= 0.9 & multiplicity <= 2 & entropy <= 0.1
        )
}

cat("=== Confident Classifications ===\n")
cat("GTDB confident reads:", sum(gtdb_classified$is_confident), "\n")
cat("SMD confident reads:", sum(smdb_classified$is_confident), "\n")
if(!is.null(pluspf_classified)) {
    cat("PlusPF confident reads:", sum(pluspf_classified$is_confident), "\n")
}
cat("GTDB confident fungal reads:", sum(gtdb_classified$is_confident & gtdb_classified$is_fungi), "\n")
cat("SMD confident fungal reads:", sum(smdb_classified$is_confident & smdb_classified$is_fungi), "\n")
if(!is.null(pluspf_classified)) {
    cat("PlusPF confident fungal reads:", sum(pluspf_classified$is_confident & pluspf_classified$is_fungi), "\n")
}
cat("\n")

# Merge on read_id to compare assignments
comparison <- gtdb_classified %>%
    select(read_id, 
           gtdb_taxid = taxid, 
           gtdb_name = name, 
           gtdb_rank = rank,
           gtdb_is_fungi = is_fungi,
           gtdb_is_bacteria = is_bacteria,
           gtdb_is_confident = is_confident,
           gtdb_consistency = consistency,
           gtdb_multiplicity = multiplicity,
           gtdb_entropy = entropy) %>%
    full_join(
        smdb_classified %>%
            select(read_id,
                   smdb_taxid = taxid,
                   smdb_name = name,
                   smdb_rank = rank,
                   smdb_is_fungi = is_fungi,
                   smdb_is_confident = is_confident,
                   smdb_consistency = consistency,
                   smdb_multiplicity = multiplicity,
                   smdb_entropy = entropy),
        by = "read_id"
    )

# Add PlusPF if available
if(!is.null(pluspf_classified)) {
    comparison <- comparison %>%
        full_join(
            pluspf_classified %>%
                select(read_id,
                       pluspf_taxid = taxid,
                       pluspf_name = name,
                       pluspf_rank = rank,
                       pluspf_is_fungi = is_fungi,
                       pluspf_is_confident = is_confident,
                       pluspf_consistency = consistency,
                       pluspf_multiplicity = multiplicity,
                       pluspf_entropy = entropy),
            by = "read_id"
        )
}

cat("=== Read Overlap ===\n")
cat("Reads in both databases:", sum(!is.na(comparison$gtdb_taxid) & !is.na(comparison$smdb_taxid)), "\n")
cat("Reads only in GTDB:", sum(!is.na(comparison$gtdb_taxid) & is.na(comparison$smdb_taxid)), "\n")
cat("Reads only in SMD:", sum(is.na(comparison$gtdb_taxid) & !is.na(comparison$smdb_taxid)), "\n\n")

# Analysis 1: Investigate "Homo sapiens" assignments in PlusPF
cat("=== ANALYSIS 1: PlusPF 'Homo sapiens' Assignments Investigation ===\n")
if(!is.null(pluspf_classified)) {
    pluspf_human <- comparison %>%
        filter(!is.na(pluspf_taxid) & 
               grepl("Homo sapiens|Homo_sapiens", pluspf_name, ignore.case = TRUE))
    
    cat("PlusPF reads assigned to 'Homo sapiens':", nrow(pluspf_human), "\n")
    
    if(nrow(pluspf_human) > 0) {
        # What are these reads assigned to in SMD?
        cat("\n--- SMD assignments for PlusPF 'Homo sapiens' reads ---\n")
        smdb_assignments_human <- pluspf_human %>%
            filter(!is.na(smdb_taxid)) %>%
            count(smdb_rank, smdb_name, sort = TRUE) %>%
            mutate(pct = n / nrow(pluspf_human) * 100)
        
        cat("  Reads with SMD assignment:", sum(!is.na(pluspf_human$smdb_taxid)), "\n")
        cat("  Reads without SMD assignment:", sum(is.na(pluspf_human$smdb_taxid)), "\n")
        
        if(nrow(smdb_assignments_human) > 0) {
            cat("\nTop 30 SMD assignments for PlusPF 'Homo sapiens' reads:\n")
            print(head(smdb_assignments_human, 30))
            
            # Check if any are actually fungi
            smdb_fungi_count <- pluspf_human %>%
                filter(!is.na(smdb_taxid) & smdb_is_fungi) %>%
                nrow()
            cat("\n  Reads that are fungi in SMD:", smdb_fungi_count, 
                sprintf("(%.2f%%)", smdb_fungi_count / nrow(pluspf_human) * 100), "\n")
        }
        
        # What are these reads assigned to in GTDB?
        cat("\n--- GTDB assignments for PlusPF 'Homo sapiens' reads ---\n")
        gtdb_assignments_human <- pluspf_human %>%
            filter(!is.na(gtdb_taxid)) %>%
            count(gtdb_rank, gtdb_name, sort = TRUE) %>%
            mutate(pct = n / nrow(pluspf_human) * 100)
        
        cat("  Reads with GTDB assignment:", sum(!is.na(pluspf_human$gtdb_taxid)), "\n")
        cat("  Reads without GTDB assignment:", sum(is.na(pluspf_human$gtdb_taxid)), "\n")
        
        if(nrow(gtdb_assignments_human) > 0) {
            cat("\nTop 30 GTDB assignments for PlusPF 'Homo sapiens' reads:\n")
            print(head(gtdb_assignments_human, 30))
        }
        
        # Quality metrics for these reads in PlusPF
        cat("\n--- PlusPF quality metrics for 'Homo sapiens' reads ---\n")
        quality_summary <- pluspf_human %>%
            summarize(
                n_reads = n(),
                mean_consistency = mean(pluspf_consistency, na.rm = TRUE),
                mean_multiplicity = mean(pluspf_multiplicity, na.rm = TRUE),
                mean_entropy = mean(pluspf_entropy, na.rm = TRUE),
                n_confident = sum(pluspf_is_confident, na.rm = TRUE),
                pct_confident = n_confident / n_reads * 100
            )
        print(quality_summary)
        
        # Check if these reads are confident in other databases
        cat("\n--- Cross-database confidence comparison ---\n")
        cross_db_conf <- pluspf_human %>%
            summarize(
                n_total = n(),
                pluspf_confident = sum(pluspf_is_confident, na.rm = TRUE),
                smdb_confident = sum(smdb_is_confident, na.rm = TRUE),
                gtdb_confident = sum(gtdb_is_confident, na.rm = TRUE),
                smdb_classified = sum(!is.na(smdb_taxid)),
                gtdb_classified = sum(!is.na(gtdb_taxid))
            )
        print(cross_db_conf)
        
        # Save read IDs for BLAST verification
        output_dir <- "data/classification/analysis_files"
        read_ids_file <- file.path(output_dir, paste0(sampleID, "_homo_sapiens_read_ids.txt"))
        
        pluspf_human %>%
            select(read_id) %>%
            write_tsv(read_ids_file, col_names = FALSE)
        
        cat("\nSaved PlusPF 'Homo sapiens' read IDs to:", read_ids_file, "\n")
        cat("Total read IDs:", nrow(pluspf_human), "\n")
    } else {
        # Initialize empty variable if no PlusPF "Homo sapiens" reads found
        pluspf_human <- data.frame()
    }
} else {
    cat("PlusPF data not available\n")
    pluspf_human <- data.frame()
}

# Analysis 2: Reads confidently assigned to fungi in GTDB - what are they in SMD?
cat("\n=== ANALYSIS 2: GTDB Confident Fungal Reads → SMD Assignments ===\n")
gtdb_confident_fungi <- comparison %>%
    filter(gtdb_is_confident & gtdb_is_fungi & !is.na(gtdb_taxid))

cat("GTDB confident fungal reads:", nrow(gtdb_confident_fungi), "\n")

if(nrow(gtdb_confident_fungi) > 0) {
    # What are these reads assigned to in SMD?
    smdb_assignments <- gtdb_confident_fungi %>%
        filter(!is.na(smdb_taxid)) %>%
        count(smdb_is_fungi, smdb_rank, smdb_name, sort = TRUE) %>%
        mutate(pct = n / nrow(gtdb_confident_fungi) * 100)
    
    cat("\nSMD assignments for GTDB confident fungal reads:\n")
    cat("  Reads with SMD assignment:", sum(!is.na(gtdb_confident_fungi$smdb_taxid)), "\n")
    cat("  Reads without SMD assignment:", sum(is.na(gtdb_confident_fungi$smdb_taxid)), "\n")
    
    if(nrow(smdb_assignments) > 0) {
        cat("\nTop SMD assignments:\n")
        print(head(smdb_assignments, 20))
        
        # Summary by kingdom/fungi status
        smdb_summary <- gtdb_confident_fungi %>%
            filter(!is.na(smdb_taxid)) %>%
            count(smdb_is_fungi) %>%
            mutate(pct = n / sum(n) * 100)
        
        cat("\nSMD assignment summary (fungi vs non-fungi):\n")
        print(smdb_summary)
    }
}

# Analysis 3: SMD confident fungal reads assigned to bacteria in GTDB
cat("\n=== ANALYSIS 3: SMD Confident Fungal Reads → GTDB Bacteria Mismatches ===\n")
smdb_confident_fungi <- comparison %>%
    filter(smdb_is_confident & smdb_is_fungi & !is.na(smdb_taxid))

cat("SMD confident fungal reads:", nrow(smdb_confident_fungi), "\n")

if(nrow(smdb_confident_fungi) > 0) {
    # Find reads assigned to bacteria in GTDB
    fungal_bacteria_mismatch <- smdb_confident_fungi %>%
        filter(!is.na(gtdb_taxid) & gtdb_is_bacteria)
    
    cat("SMD confident fungal reads assigned to bacteria in GTDB:", nrow(fungal_bacteria_mismatch), "\n")
    cat("Percentage:", sprintf("%.2f%%", nrow(fungal_bacteria_mismatch) / nrow(smdb_confident_fungi) * 100), "\n")
    
    if(nrow(fungal_bacteria_mismatch) > 0) {
        # GTDB rank distribution
        rank_dist <- fungal_bacteria_mismatch %>%
            count(gtdb_rank, sort = TRUE) %>%
            mutate(pct = n / nrow(fungal_bacteria_mismatch) * 100)
        
        cat("\nGTDB rank distribution for mismatched reads:\n")
        print(rank_dist)
        
        # Top GTDB assignments
        gtdb_assignments <- fungal_bacteria_mismatch %>%
            count(gtdb_rank, gtdb_name, sort = TRUE) %>%
            mutate(pct = n / nrow(fungal_bacteria_mismatch) * 100) %>%
            head(20)
        
        cat("\nTop GTDB assignments:\n")
        print(gtdb_assignments)
        
        # Top SMD fungal taxa being misclassified
        smdb_fungal_taxa <- fungal_bacteria_mismatch %>%
            count(smdb_rank, smdb_name, sort = TRUE) %>%
            mutate(pct = n / nrow(fungal_bacteria_mismatch) * 100) %>%
            head(20)
        
        cat("\nTop SMD fungal taxa being misclassified:\n")
        print(smdb_fungal_taxa)
        
        # Quality metrics comparison
        quality_comparison <- fungal_bacteria_mismatch %>%
            summarize(
                n_reads = n(),
                smdb_mean_consistency = mean(smdb_consistency, na.rm = TRUE),
                smdb_mean_multiplicity = mean(smdb_multiplicity, na.rm = TRUE),
                smdb_mean_entropy = mean(smdb_entropy, na.rm = TRUE),
                smdb_confident_pct = sum(smdb_is_confident, na.rm = TRUE) / n_reads * 100,
                gtdb_mean_consistency = mean(gtdb_consistency, na.rm = TRUE),
                gtdb_mean_multiplicity = mean(gtdb_multiplicity, na.rm = TRUE),
                gtdb_mean_entropy = mean(gtdb_entropy, na.rm = TRUE),
                gtdb_confident_pct = sum(gtdb_is_confident, na.rm = TRUE) / n_reads * 100
            )
        
        cat("\nQuality metrics comparison:\n")
        print(quality_comparison)
        
        # Save read IDs for BLAST verification
        read_ids_file <- file.path(output_dir, paste0(sampleID, "_smdb_fungal_gtdb_bacteria_read_ids.txt"))
        
        fungal_bacteria_mismatch %>%
            select(read_id) %>%
            write_tsv(read_ids_file, col_names = FALSE)
        
        cat("\nSaved SMD fungal → GTDB bacteria mismatch read IDs to:", read_ids_file, "\n")
        cat("Total read IDs:", nrow(fungal_bacteria_mismatch), "\n")
    }
} else {
    # Initialize empty variables if no SMD confident fungi found
    smdb_confident_fungi <- data.frame()
    fungal_bacteria_mismatch <- data.frame()
}

# Analysis 3b: Sample 500 confident SMD fungal reads for BLAST verification
cat("\n=== ANALYSIS 3b: Sample SMD Confident Fungal Reads for BLAST Verification ===\n")
if(exists("smdb_confident_fungi") && is.data.frame(smdb_confident_fungi) && nrow(smdb_confident_fungi) > 0) {
    # Sample 500 reads (or all if fewer than 500)
    sample_size <- min(500, nrow(smdb_confident_fungi))
    smdb_confident_fungi_sample <- smdb_confident_fungi %>%
        sample_n(sample_size)
    
    cat("Sampled", nrow(smdb_confident_fungi_sample), "SMD confident fungal reads for BLAST verification\n")
    cat("(from", nrow(smdb_confident_fungi), "total confident fungal reads)\n")
    
    # Save read IDs for BLAST verification
    smdb_fungi_verification_file <- file.path(output_dir, paste0(sampleID, "_smdb_confident_fungi_read_ids.txt"))
    
    smdb_confident_fungi_sample %>%
        select(read_id) %>%
        write_tsv(smdb_fungi_verification_file, col_names = FALSE)
    
    cat("Saved SMD confident fungal read IDs to:", smdb_fungi_verification_file, "\n")
    cat("Total read IDs:", nrow(smdb_confident_fungi_sample), "\n")
} else {
    cat("No SMD confident fungal reads found. Skipping verification sample.\n")
}

# Analysis 4: Reads confidently assigned to bacteria in SMD - what are they in GTDB?
cat("\n=== ANALYSIS 4: SMD Confident Bacterial Reads → GTDB Assignments ===\n")
smdb_confident_bacteria <- comparison %>%
    filter(smdb_is_confident & 
           !smdb_is_fungi & 
           !is.na(smdb_taxid) &
           (smdb_rank == "k" & grepl("Bacteria|bacteria", smdb_name) |
            grepl("^[sgfoc]__", smdb_name)))  # Species, genus, family, order, class level

cat("SMD confident bacterial reads (approx):", nrow(smdb_confident_bacteria), "\n")

if(nrow(smdb_confident_bacteria) > 0) {
    # Sample a subset for analysis (files can be very large)
    if(nrow(smdb_confident_bacteria) > 10000) {
        smdb_confident_bacteria_sample <- smdb_confident_bacteria %>% sample_n(10000)
        cat("  (Sampling 10,000 reads for analysis)\n")
    } else {
        smdb_confident_bacteria_sample <- smdb_confident_bacteria
    }
    
    gtdb_assignments_bact <- smdb_confident_bacteria_sample %>%
        filter(!is.na(gtdb_taxid)) %>%
        count(gtdb_rank, sort = TRUE) %>%
        mutate(pct = n / nrow(smdb_confident_bacteria_sample) * 100)
    
    cat("\nGTDB rank distribution for SMD confident bacterial reads:\n")
    print(gtdb_assignments_bact)
}

# Analysis 5: Overall taxonomic rank distribution
cat("\n=== ANALYSIS 5: Overall Taxonomic Rank Distribution ===\n")
gtdb_ranks <- gtdb_classified %>%
    count(rank, sort = TRUE) %>%
    mutate(pct = n / nrow(gtdb_classified) * 100,
           db = "GTDB")

smdb_ranks <- smdb_classified %>%
    count(rank, sort = TRUE) %>%
    mutate(pct = n / nrow(smdb_classified) * 100,
           db = "SMD")

rank_comparison <- bind_rows(gtdb_ranks, smdb_ranks) %>%
    pivot_wider(names_from = db, values_from = c(n, pct), values_fill = 0)

cat("Rank distribution comparison:\n")
print(rank_comparison)

# Save detailed comparison for specific reads
output_dir <- "data/classification/analysis_files"
output_file <- file.path(output_dir, paste0("taxonomic_assignment_comparison_", sampleID, ".csv"))
write_csv(comparison, output_file)

# Generate summary report
summary_file <- file.path(output_dir, paste0("taxonomic_assignment_comparison_", sampleID, "_summary.txt"))

homo_sapiens_read_ids_file <- file.path(output_dir, paste0(sampleID, "_homo_sapiens_read_ids.txt"))
fungal_bacteria_read_ids_file <- file.path(output_dir, paste0(sampleID, "_smdb_fungal_gtdb_bacteria_read_ids.txt"))

summary_lines <- c(
    paste0("Taxonomic Assignment Comparison: GTDB vs SoilMicrobeDB vs PlusPF"),
    paste0("==================================================================="),
    "",
    paste0("Sample: ", sampleID),
    paste0("Date: ", Sys.Date()),
    "",
    "=== Analysis Summary ===",
    ""
)

# Check for PlusPF "Homo sapiens" misclassifications (reuse variable from Analysis 1)
if(exists("pluspf_human") && is.data.frame(pluspf_human) && nrow(pluspf_human) > 0) {
    summary_lines <- c(summary_lines,
        "Analysis 1: PlusPF 'Homo sapiens' Assignments",
        paste0("  Total reads: ", nrow(pluspf_human)),
        paste0("  Assigned in SMD: ", sum(!is.na(pluspf_human$smdb_taxid))),
        paste0("  Assigned in GTDB: ", sum(!is.na(pluspf_human$gtdb_taxid))),
        paste0("  Confident in PlusPF: ", sum(pluspf_human$pluspf_is_confident, na.rm = TRUE)),
        paste0("  Confident in SMD: ", sum(pluspf_human$smdb_is_confident, na.rm = TRUE)),
        paste0("  Confident in GTDB: ", sum(pluspf_human$gtdb_is_confident, na.rm = TRUE)),
        "",
        paste0("  Read IDs saved to: ", homo_sapiens_read_ids_file),
        ""
    )
}

# Check for SMD fungal → GTDB bacteria mismatches (reuse variables from Analysis 3)
if(exists("smdb_confident_fungi") && exists("fungal_bacteria_mismatch") && 
   is.data.frame(fungal_bacteria_mismatch) && nrow(fungal_bacteria_mismatch) > 0) {
    summary_lines <- c(summary_lines,
        "Analysis 3: SMD Fungal → GTDB Bacteria Mismatches",
        paste0("  SMD confident fungal reads: ", nrow(smdb_confident_fungi)),
        paste0("  Assigned to bacteria in GTDB: ", nrow(fungal_bacteria_mismatch)),
        paste0("  Percentage: ", sprintf("%.2f%%", nrow(fungal_bacteria_mismatch) / nrow(smdb_confident_fungi) * 100)),
        "",
        paste0("  Read IDs saved to: ", fungal_bacteria_read_ids_file),
        ""
    )
}

read_id_files <- c()
if(file.exists(homo_sapiens_read_ids_file)) {
    read_id_files <- c(read_id_files, homo_sapiens_read_ids_file)
}
if(file.exists(fungal_bacteria_read_ids_file)) {
    read_id_files <- c(read_id_files, fungal_bacteria_read_ids_file)
}

summary_lines <- c(summary_lines,
    "=== Output Files ===",
    paste0("  Full comparison: ", output_file),
    paste0("  Summary report: ", summary_file),
    "",
    "=== Next Steps ===",
    "  Both misclassification types identified above can be verified via BLAST:",
    "    1. Run script 16: python3 scripts/summarize_outputs/16_extract_reads_to_fasta.py \\",
    "         {sampleID} <R1.fastq.gz> <R2.fastq.gz>",
    "    2. Run script 17: python3 scripts/summarize_outputs/17_blast_reads_ncbi.py \\",
    "         {sampleID} [sample_size]",
    "    3. Run script 18: Rscript scripts/summarize_outputs/18_analyze_blast_results.r \\",
    "         {sampleID}",
    ""
)

if(length(read_id_files) > 0) {
    summary_lines <- c(summary_lines,
        "  Read ID files generated:",
        paste0("    - ", read_id_files, collapse = "\n")
    )
} else {
    summary_lines <- c(summary_lines,
        "  Read ID files: (none - no mismatches found)"
    )
}

writeLines(summary_lines, summary_file)
cat("\n=== Summary saved to:", summary_file, "===\n")
