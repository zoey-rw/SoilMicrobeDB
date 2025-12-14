# Query FUNguild API for functional guild assignments
# Uses all fungal genera from metagenome classification

library(tidyverse)
library(data.table)
library(jsonlite)

# Read in genus-level abundances to get all fungal genera
bracken_genus <- fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_genus_merged_lineage.csv", nThread = 8)

# Get unique fungal genera (those with Eukaryota in lineage)
unique_genera <- bracken_genus %>% 
    select(lineage, name) %>% 
    distinct()

genus_fungi <- unique_genera[grepl("Eukary", unique_genera$lineage),]$name

cat("Found", length(genus_fungi), "fungal genera to query FUNguild\n")
cat("Querying FUNguild API...\n")

# Query FUNguild using taxon names
results_list_full <- lapply(genus_fungi, function(taxon) {
    tryCatch({
        jsonlite::fromJSON(paste0("https://www.mycoportal.org/funguild/services/api/db_return.php?qDB=funguild_db&qField=taxon&qText=", taxon))
    }, error = function(e) {
        cat("  Warning: Failed to query", taxon, "-", e$message, "\n")
        return(NULL)
    })
})

# Remove NULL results
results_list_full <- results_list_full[!sapply(results_list_full, is.null)]

if(length(results_list_full) == 0) {
    stop("❌ ERROR: No FUNguild results returned!\n",
         "   Check internet connection and FUNguild API availability.")
}

funguild_results_full <- do.call(rbind, results_list_full) %>% 
    as.data.frame() %>% 
    select(taxon, trophicMode, guild, confidenceRanking) 

# Assign functional groups based on guild assignments
funguild_results_full <- funguild_results_full %>%  
    mutate(functional_group = 
           ifelse(grepl("Ecto", guild), 
                  "Ectomycorrhizae", 
                  ifelse(grepl("Arbusc", guild), 
                         "Arbuscular mycorrhizae", 
                         ifelse(grepl("Pathogen", guild), 
                                "Pathogen",
                                ifelse(grepl("Endophyte", guild), 
                                       "Endophyte",
                                       ifelse(grepl("Sapro", guild), 
                                              "Other saprotroph",
                                              NA))))))

funguild_results_full$name <- tolower(funguild_results_full$taxon)

write_csv(funguild_results_full, "data/classification/analysis_files/FUNguild_assignments.csv")
cat("✅ Saved FUNguild assignments to: data/classification/analysis_files/FUNguild_assignments.csv\n")
cat("   Total assignments:", nrow(funguild_results_full), "\n")
