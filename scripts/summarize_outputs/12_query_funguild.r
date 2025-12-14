# Query FUNguild API for functional guild assignments of fungal genera
# Input: Genus-level abundances from soil_microbe_db_genus_merged_lineage.csv
#        Filters to genera with Eukaryota in lineage
# Output: FUNguild_assignments.csv with trophic mode, guild, and functional group assignments
#         Functional groups: Ectomycorrhizae, Arbuscular mycorrhizae, Pathogen, Endophyte, Other saprotroph
#         Queries FUNguild API: https://www.mycoportal.org/funguild/services/api/

library(tidyverse)
library(data.table)
library(jsonlite)

bracken_genus <- fread("data/classification/taxonomic_rank_summaries/soil_microbe_db_genus_merged_lineage.csv", nThread = 8)

genus_fungi <- bracken_genus %>% 
    select(lineage, name) %>% 
    distinct() %>%
    filter(grepl("Eukary", lineage)) %>%
    pull(name)

results_list_full <- lapply(genus_fungi, function(taxon) {
    tryCatch({
        jsonlite::fromJSON(paste0("https://www.mycoportal.org/funguild/services/api/db_return.php?qDB=funguild_db&qField=taxon&qText=", taxon))
    }, error = function(e) NULL)
})

results_list_full <- results_list_full[!sapply(results_list_full, is.null)]

if(length(results_list_full) == 0) {
    stop("No FUNguild results returned")
}

funguild_results_full <- do.call(rbind, results_list_full) %>% 
    as.data.frame() %>% 
    select(taxon, trophicMode, guild, confidenceRanking) %>%
    mutate(functional_group = 
           ifelse(grepl("Ecto", guild), "Ectomycorrhizae", 
                  ifelse(grepl("Arbusc", guild), "Arbuscular mycorrhizae", 
                         ifelse(grepl("Pathogen", guild), "Pathogen",
                                ifelse(grepl("Endophyte", guild), "Endophyte",
                                       ifelse(grepl("Sapro", guild), "Other saprotroph", NA)))))),
           name = tolower(taxon))

write_csv(funguild_results_full, "data/classification/analysis_files/FUNguild_assignments.csv")
