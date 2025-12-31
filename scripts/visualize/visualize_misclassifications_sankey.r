#!/usr/bin/env Rscript
# Visualize misclassifications using Sankey diagrams
# Shows incorrect vs correct proportions for:
#   - PlusPF "Homo sapiens" misclassifications → SMD/GTDB assignments
#   - SMD confident fungal reads → GTDB bacteria mismatches
#
# Usage: Rscript scripts/visualize/visualize_misclassifications_sankey.r [sampleID]
#
# Input:  taxonomic_assignment_comparison_{sampleID}.csv (from script 15)
# Output: 
#   - {sampleID}_homo_sapiens_sankey.html (Sankey diagram for PlusPF misclassifications)
#   - {sampleID}_fungal_bacteria_sankey.html (Sankey diagram for SMD fungal → GTDB bacteria mismatches)

library(tidyverse)
library(plotly)
library(htmlwidgets)

# Try to load kaleido for static image export (optional)
has_kaleido <- requireNamespace("kaleido", quietly = TRUE)
if(!has_kaleido) {
    cat("Note: kaleido package not available. Only HTML files will be generated.\n")
    cat("To install: install.packages('kaleido')\n")
}


# Sample to analyze (default to ORNL_046-O-20170621-COMP)
sampleID <- if(length(commandArgs(trailingOnly = TRUE)) > 0) {
    commandArgs(trailingOnly = TRUE)[1]
} else {
    "ORNL_046-O-20170621-COMP"
}

cat("=== Visualizing Misclassifications with Sankey Diagrams ===\n")
cat("Sample:", sampleID, "\n\n")

# File paths
input_file <- file.path("data/classification/analysis_files", 
                        paste0("taxonomic_assignment_comparison_", sampleID, ".csv"))
output_dir <- "manuscript_figures"

if(!file.exists(input_file)) {
    stop("Comparison file not found: ", input_file, 
         "\nPlease run script 15 first: Rscript scripts/summarize_outputs/15_compare_taxonomic_assignments.r ", sampleID)
}

cat("Reading comparison data...\n")
comparison <- read_csv(input_file, show_col_types = FALSE)

# Define fungal phyla for identification
fungal_phyla <- c("Ascomycota", "Basidiomycota", "Blastocladiomycota", 
                  "Chytridiomycota", "Cryptomycota", "Mucoromycota", 
                  "Microsporidia", "Olpidiomycota", "Zoopagomycota")
fungal_pattern <- paste0("(", paste(fungal_phyla, collapse = "|"), ")")

# Helper function to create Sankey diagram data
create_sankey_data <- function(source_nodes, target_nodes, values, 
                                source_labels, target_labels) {
    # Create node mapping
    all_nodes <- unique(c(source_nodes, target_nodes))
    node_map <- setNames(0:(length(all_nodes) - 1), all_nodes)
    
    # Create links
    links <- data.frame(
        source = node_map[source_nodes],
        target = node_map[target_nodes],
        value = values
    )
    
    # Create node labels
    node_labels <- c(
        setNames(source_labels, all_nodes[1:length(unique(source_nodes))]),
        setNames(target_labels, all_nodes[(length(unique(source_nodes)) + 1):length(all_nodes)])
    )
    node_labels <- node_labels[all_nodes]
    
    list(
        links = links,
        nodes = data.frame(label = node_labels),
        node_map = node_map
    )
}

# Analysis 1: PlusPF "Homo sapiens" misclassifications
cat("\n=== Creating Sankey diagram for PlusPF 'Homo sapiens' misclassifications ===\n")
pluspf_human <- comparison %>%
    filter(!is.na(pluspf_taxid) & 
           grepl("Homo sapiens|Homo_sapiens", pluspf_name, ignore.case = TRUE))

if(nrow(pluspf_human) > 0) {
    cat("PlusPF 'Homo sapiens' reads:", nrow(pluspf_human), "\n")
    
    # Categorize SMD assignments
    smdb_categories <- pluspf_human %>%
        mutate(
            smdb_category = case_when(
                is.na(smdb_taxid) ~ "Not classified in SMD",
                smdb_is_fungi ~ "Fungi (SMD)",
                grepl("Eukaryota|eukaryota", smdb_name) | 
                (smdb_rank == "k" & grepl("Eukaryota", smdb_name)) ~ "Eukaryota (SMD)",
                grepl("Bacteria|bacteria", smdb_name) | 
                (smdb_rank == "k" & grepl("Bacteria", smdb_name)) ~ "Bacteria (SMD)",
                TRUE ~ "Other (SMD)"
            )
        ) %>%
        count(smdb_category) %>%
        mutate(
            source = "PlusPF: Homo sapiens",
            pct = n / nrow(pluspf_human) * 100
        )
    
    # Categorize GTDB assignments
    gtdb_categories <- pluspf_human %>%
        mutate(
            gtdb_category = case_when(
                is.na(gtdb_taxid) ~ "Not classified in GTDB",
                grepl("Eukaryota|eukaryota", gtdb_name) | 
                (gtdb_rank == "k" & grepl("Eukaryota", gtdb_name)) ~ "Eukaryota (GTDB)",
                grepl("Bacteria|bacteria", gtdb_name) | 
                (gtdb_rank == "k" & grepl("Bacteria", gtdb_name)) ~ "Bacteria (GTDB)",
                TRUE ~ "Other (GTDB)"
            )
        ) %>%
        count(gtdb_category) %>%
        mutate(
            source = "PlusPF: Homo sapiens",
            pct = n / nrow(pluspf_human) * 100
        )
    
    # Create Sankey data for SMD
    source_nodes_smd <- rep("PlusPF: Homo sapiens", nrow(smdb_categories))
    target_nodes_smd <- smdb_categories$smdb_category
    values_smd <- smdb_categories$n
    
    # Create Sankey data for GTDB
    source_nodes_gtdb <- rep("PlusPF: Homo sapiens", nrow(gtdb_categories))
    target_nodes_gtdb <- gtdb_categories$gtdb_category
    values_gtdb <- gtdb_categories$n
    
    # Combine for visualization
    all_nodes <- unique(c("PlusPF: Homo sapiens", 
                          unique(smdb_categories$smdb_category),
                          unique(gtdb_categories$gtdb_category)))
    node_map <- setNames(0:(length(all_nodes) - 1), all_nodes)
    
    # Create links
    links_smd <- data.frame(
        source = node_map["PlusPF: Homo sapiens"],
        target = node_map[target_nodes_smd],
        value = values_smd,
        label = paste0(round(smdb_categories$pct, 1), "%")
    )
    
    links_gtdb <- data.frame(
        source = node_map["PlusPF: Homo sapiens"],
        target = node_map[target_nodes_gtdb],
        value = values_gtdb,
        label = paste0(round(gtdb_categories$pct, 1), "%")
    )
    
    # Create Sankey diagram
    fig <- plot_ly(
        type = "sankey",
        orientation = "h",
        node = list(
            label = all_nodes,
            pad = 15,
            thickness = 20,
            line = list(color = "black", width = 0.5),
            color = rep("rgba(200, 200, 200, 0.8)", length(all_nodes))
        ),
        link = list(
            source = c(links_smd$source, links_gtdb$source),
            target = c(links_smd$target, links_gtdb$target),
            value = c(links_smd$value, links_gtdb$value),
            color = c(rep("rgba(31, 119, 180, 0.4)", nrow(links_smd)),
                     rep("rgba(255, 127, 14, 0.4)", nrow(links_gtdb)))
        )
    ) %>%
        layout(
            title = list(
                text = paste0("PlusPF 'Homo sapiens' Misclassifications: Assignments in SMD and GTDB<br>",
                            "Sample: ", sampleID, " (n = ", nrow(pluspf_human), " reads)"),
                x = 0.5,
                xanchor = "center"
            ),
            font = list(size = 12)
        )
    
    output_file_html <- file.path(output_dir, paste0(sampleID, "_homo_sapiens_sankey.html"))
    htmlwidgets::saveWidget(fig, output_file_html, selfcontained = TRUE)
    cat("Saved Sankey diagram (HTML) to:", output_file_html, "\n")
    
    # Also save as static image if kaleido is available
    if(has_kaleido) {
        output_file_png <- file.path(output_dir, paste0(sampleID, "_homo_sapiens_sankey.png"))
        tryCatch({
            kaleido::kaleido(fig, file = output_file_png, width = 1200, height = 800)
            cat("Saved Sankey diagram (PNG) to:", output_file_png, "\n")
        }, error = function(e) {
            cat("Note: PNG export failed:", conditionMessage(e), "\n")
        })
    }
} else {
    cat("No PlusPF 'Homo sapiens' reads found.\n")
}

# Analysis 2: SMD confident fungal reads → GTDB bacteria mismatches
cat("\n=== Creating Sankey diagram for SMD fungal → GTDB bacteria mismatches ===\n")

# Identify confident classifications
comparison <- comparison %>%
    mutate(
        smdb_is_confident = smdb_consistency >= 0.9 & smdb_multiplicity <= 2 & smdb_entropy <= 0.1,
        gtdb_is_confident = gtdb_consistency >= 0.9 & gtdb_multiplicity <= 2 & gtdb_entropy <= 0.1
    )

# Identify fungi in SMD
comparison <- comparison %>%
    mutate(
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
        )
    )

smdb_confident_fungi <- comparison %>%
    filter(smdb_is_confident & smdb_is_fungi & !is.na(smdb_taxid))

if(nrow(smdb_confident_fungi) > 0) {
    cat("SMD confident fungal reads:", nrow(smdb_confident_fungi), "\n")
    
    # Find reads assigned to bacteria in GTDB
    fungal_bacteria_mismatch <- smdb_confident_fungi %>%
        filter(!is.na(gtdb_taxid) & gtdb_is_bacteria)
    
    cat("Assigned to bacteria in GTDB:", nrow(fungal_bacteria_mismatch), 
        sprintf("(%.2f%%)\n", nrow(fungal_bacteria_mismatch) / nrow(smdb_confident_fungi) * 100))
    
    if(nrow(fungal_bacteria_mismatch) > 0) {
        # Categorize GTDB bacterial assignments
        gtdb_bact_categories <- fungal_bacteria_mismatch %>%
            mutate(
                gtdb_category = case_when(
                    is.na(gtdb_taxid) ~ "Not classified in GTDB",
                    gtdb_rank == "k" ~ "Bacteria (Kingdom, GTDB)",
                    gtdb_rank == "p" ~ paste0("Bacteria (Phylum: ", gtdb_name, ", GTDB)"),
                    gtdb_rank == "c" ~ paste0("Bacteria (Class: ", gtdb_name, ", GTDB)"),
                    gtdb_rank == "o" ~ paste0("Bacteria (Order: ", gtdb_name, ", GTDB)"),
                    gtdb_rank == "f" ~ paste0("Bacteria (Family: ", gtdb_name, ", GTDB)"),
                    gtdb_rank == "g" ~ paste0("Bacteria (Genus: ", gtdb_name, ", GTDB)"),
                    gtdb_rank == "s" ~ paste0("Bacteria (Species: ", gtdb_name, ", GTDB)"),
                    TRUE ~ paste0("Bacteria (", gtdb_rank, ": ", gtdb_name, ", GTDB)")
                )
            ) %>%
            count(gtdb_category, sort = TRUE) %>%
            mutate(
                source = "SMD: Confident Fungi",
                pct = n / nrow(fungal_bacteria_mismatch) * 100
            )
        
        # Limit to top categories for readability
        if(nrow(gtdb_bact_categories) > 10) {
            top_categories <- head(gtdb_bact_categories, 9)
            other_count <- sum(tail(gtdb_bact_categories, -9)$n)
            other_pct <- sum(tail(gtdb_bact_categories, -9)$pct)
            gtdb_bact_categories <- bind_rows(
                top_categories,
                data.frame(
                    gtdb_category = "Bacteria (Other, GTDB)",
                    n = other_count,
                    source = "SMD: Confident Fungi",
                    pct = other_pct
                )
            )
        }
        
        # Also show reads correctly classified as fungi in GTDB
        correctly_classified <- smdb_confident_fungi %>%
            filter(!is.na(gtdb_taxid) & !gtdb_is_bacteria) %>%
            mutate(
                gtdb_category = case_when(
                    grepl("Eukaryota|eukaryota", gtdb_name) | 
                    (gtdb_rank == "k" & grepl("Eukaryota", gtdb_name)) ~ "Eukaryota/Fungi (GTDB)",
                    TRUE ~ "Other (GTDB)"
                )
            ) %>%
            count(gtdb_category) %>%
            mutate(
                source = "SMD: Confident Fungi",
                pct = n / nrow(smdb_confident_fungi) * 100
            )
        
        # Create Sankey data
        all_nodes <- unique(c("SMD: Confident Fungi",
                              gtdb_bact_categories$gtdb_category,
                              correctly_classified$gtdb_category))
        node_map <- setNames(0:(length(all_nodes) - 1), all_nodes)
        
        # Links for mismatches (bacteria)
        links_mismatch <- data.frame(
            source = node_map["SMD: Confident Fungi"],
            target = node_map[gtdb_bact_categories$gtdb_category],
            value = gtdb_bact_categories$n,
            label = paste0(round(gtdb_bact_categories$pct, 1), "%")
        )
        
        # Links for correct classifications
        links_correct <- data.frame(
            source = node_map["SMD: Confident Fungi"],
            target = node_map[correctly_classified$gtdb_category],
            value = correctly_classified$n,
            label = paste0(round(correctly_classified$pct, 1), "%")
        )
        
        # Create node colors
        node_colors <- rep("rgba(200, 200, 200, 0.8)", length(all_nodes))
        node_colors[all_nodes == "SMD: Confident Fungi"] <- "rgba(31, 119, 180, 0.8)"
        
        # Create Sankey diagram
        fig <- plot_ly(
            type = "sankey",
            orientation = "h",
            node = list(
                label = all_nodes,
                pad = 15,
                thickness = 20,
                line = list(color = "black", width = 0.5),
                color = node_colors
            ),
            link = list(
                source = c(links_mismatch$source, links_correct$source),
                target = c(links_mismatch$target, links_correct$target),
                value = c(links_mismatch$value, links_correct$value),
                color = c(rep("rgba(227, 26, 28, 0.4)", nrow(links_mismatch)),
                         rep("rgba(51, 160, 44, 0.4)", nrow(links_correct)))
            )
        ) %>%
            layout(
                title = list(
                    text = paste0("SMD Confident Fungal Reads: GTDB Assignments<br>",
                                "Sample: ", sampleID, 
                                " (n = ", nrow(smdb_confident_fungi), " SMD confident fungi, ",
                                nrow(fungal_bacteria_mismatch), " mismatched as bacteria)"),
                    x = 0.5,
                    xanchor = "center"
                ),
                font = list(size = 12)
            )
        
        output_file_html <- file.path(output_dir, paste0(sampleID, "_fungal_bacteria_sankey.html"))
        htmlwidgets::saveWidget(fig, output_file_html, selfcontained = TRUE)
        cat("Saved Sankey diagram (HTML) to:", output_file_html, "\n")
        
        # Also save as static image if kaleido is available
        if(has_kaleido) {
            output_file_png <- file.path(output_dir, paste0(sampleID, "_fungal_bacteria_sankey.png"))
            tryCatch({
                kaleido::kaleido(fig, file = output_file_png, width = 1200, height = 800)
                cat("Saved Sankey diagram (PNG) to:", output_file_png, "\n")
            }, error = function(e) {
                cat("Note: PNG export failed:", conditionMessage(e), "\n")
            })
        }
    } else {
        cat("No SMD fungal → GTDB bacteria mismatches found.\n")
    }
} else {
    cat("No SMD confident fungal reads found.\n")
}

cat("\n=== Visualization complete ===\n")
