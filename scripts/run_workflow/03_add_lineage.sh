#!/bin/bash
# Merge Bracken outputs and add lineage information
# Processes all databases automatically at species, genus, domain, and phylum levels
#
# Usage: bash scripts/run_workflow/03_add_lineage.sh
#
# Input:  *.b2 files from Step 2
# Output: *_merged_lineage.csv files

module load miniconda
conda activate struo2
cd /projectnb/frpmars/soil_microbe_db/

# Set up paths - check both new and old locations for compatibility
if [ -d "/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/02_bracken_output" ]; then
    bracken_dir=/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/02_bracken_output
    summary_dir=/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files
elif [ -d "/projectnb/frpmars/soil_microbe_db/data/classification/02_bracken_output" ]; then
    bracken_dir=/projectnb/frpmars/soil_microbe_db/data/classification/02_bracken_output
    summary_dir=/projectnb/frpmars/soil_microbe_db/data/classification/taxonomic_rank_summaries
else
    echo "ERROR: Could not find bracken output directory!"
    echo "Checked: /projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/02_bracken_output"
    echo "Checked: /projectnb/frpmars/soil_microbe_db/data/classification/02_bracken_output"
    exit 1
fi

# Create summary directory if it doesn't exist (flat structure, no subdirectories)
mkdir -p ${summary_dir}

# Function to merge and add lineage
merge_and_add_lineage() {
    local DBNAME=$1
    local DBDIR=$2
    local DB_taxonomy_dir=$3
    local b2_pattern=$4
    local rank=$5  # species, genus, domain, or phylum
    
    # All ranks use consistent naming: ${DBNAME}_${rank}_merged.csv
    # Flat structure (no subdirectories)
    merged_file=${summary_dir}/${DBNAME}_${rank}_merged.csv
    lineage_file=${summary_dir}/${DBNAME}_${rank}_merged_lineage.csv
    
    echo "Processing $DBNAME at $rank level..."
    
    # Expand glob pattern to find actual files
    # Files may have single or double underscores after database name (e.g., _filtered or __filtered)
    # Create alternate patterns to try both single and double underscores
    # Examples: *soil_microbe_db_filtered.b2 -> *soil_microbe_db__filtered.b2
    #           *pluspf__filtered.b2 -> *pluspf_filtered.b2
    b2_pattern_alt=""
    if [[ "$b2_pattern" == *"${DBNAME}_"* ]]; then
        # Pattern has single underscore after DBNAME, try double underscore
        b2_pattern_alt=$(echo "$b2_pattern" | sed "s/${DBNAME}_/${DBNAME}__/")
    elif [[ "$b2_pattern" == *"${DBNAME}__"* ]]; then
        # Pattern has double underscore after DBNAME, try single underscore
        b2_pattern_alt=$(echo "$b2_pattern" | sed "s/${DBNAME}__/${DBNAME}_/")
    elif [[ "$b2_pattern" == *"${DBNAME}.b2" ]]; then
        # For patterns like *gtdb_207_unfiltered.b2, add double underscore before .b2
        b2_pattern_alt=$(echo "$b2_pattern" | sed "s/${DBNAME}\.b2/${DBNAME}__.b2/")
    fi
    
    # Search for both patterns (single and double underscore) if alternate exists
    if [ -n "$b2_pattern_alt" ] && [ "$b2_pattern_alt" != "$b2_pattern" ]; then
        b2_files=$(find ${bracken_dir} \( -name "${b2_pattern}" -o -name "${b2_pattern_alt}" \) -type f 2>/dev/null | sort)
    else
        b2_files=$(find ${bracken_dir} -name "${b2_pattern}" -type f 2>/dev/null | sort)
    fi
    
    if [ -z "$b2_files" ]; then
        echo "  WARNING: No files found matching pattern ${b2_pattern} in ${bracken_dir}"
        # Diagnostic: show what .b2 files actually exist
        all_b2_files=$(find ${bracken_dir} -name "*.b2" -type f 2>/dev/null | head -5)
        if [ -n "$all_b2_files" ]; then
            echo "  Diagnostic: Found .b2 files in directory (showing first 5 examples):"
            echo "$all_b2_files" | while read -r file; do
                echo "    $(basename "$file")"
            done
        else
            echo "  Diagnostic: No .b2 files found in ${bracken_dir} at all"
        fi
        
        # Check if there are filtered_kraken.kreport files that could be used to generate .b2 files
        # This helps identify samples that need Bracken to be run
        if [ "$rank" = "species" ]; then
            kreport_pattern="*${DBNAME}_filtered_kraken.kreport"
        elif [ "$rank" = "genus" ]; then
            kreport_pattern="*${DBNAME}_genus_filtered_kraken.kreport"
        elif [ "$rank" = "domain" ]; then
            kreport_pattern="*${DBNAME}_domain_filtered_kraken.kreport"
        elif [ "$rank" = "phylum" ]; then
            kreport_pattern="*${DBNAME}_phylum_filtered_kraken.kreport"
        else
            kreport_pattern="*${DBNAME}*filtered_kraken.kreport"
        fi
        
        # Try alternate patterns for different naming conventions
        if [ "$DBNAME" = "gtdb_207_unfiltered" ]; then
            kreport_pattern_alt="*gtdb_207_unfiltered*filtered_kraken.kreport"
        elif [ "$DBNAME" = "gtdb_207" ]; then
            kreport_pattern_alt="*gtdb_207*filtered_kraken.kreport"
        elif [ "$DBNAME" = "pluspf" ]; then
            kreport_pattern_alt="*pluspf*filtered_kraken.kreport"
        else
            kreport_pattern_alt="$kreport_pattern"
        fi
        
        potential_kreports=$(find ${bracken_dir} \( -name "${kreport_pattern}" -o -name "${kreport_pattern_alt}" \) -type f 2>/dev/null | head -5)
        if [ -n "$potential_kreports" ]; then
            echo "  Note: Found filtered_kraken.kreport files that could be used to generate .b2 files:"
            echo "$potential_kreports" | while read -r file; do
                echo "    $(basename "$file")"
            done
            echo "  These samples may need Bracken to be run (script 02_run_architeuthis.sh) to generate .b2 files"
        fi
        
        echo "  Skipping merge for $DBNAME at $rank level"
        return
    fi
    
    # Count files found
    file_count=$(echo "$b2_files" | wc -l)
    echo "  Found $file_count file(s) matching pattern"
    
    # Check if merge is needed (file doesn't exist OR newer .b2 files exist)
    needs_merge=false
    
    if [ ! -f "$lineage_file" ]; then
        needs_merge=true
        echo "  Merged file doesn't exist, will create"
    else
        # Check if any .b2 file is newer than the merged file
        # Suppress all output except the modification time number
        merged_file_time=$(stat -f %m "$lineage_file" 2>&1 | grep -oE '^[0-9]+$' || stat -c %Y "$lineage_file" 2>&1 | grep -oE '^[0-9]+$' || echo 0)
        
        for b2_file in $b2_files; do
            b2_file_time=$(stat -f %m "$b2_file" 2>&1 | grep -oE '^[0-9]+$' || stat -c %Y "$b2_file" 2>&1 | grep -oE '^[0-9]+$' || echo 0)
            if [ "$b2_file_time" -gt "$merged_file_time" ]; then
                needs_merge=true
                echo "  Found newer .b2 file(s), will re-merge"
                break
            fi
        done
        
        # Check if existing file has incorrect sample_ids (ending with _gtdb_207_filtered instead of _gtdb_207)
        # This can happen if files were created before we removed the normalization logic
        if [ "$needs_merge" = false ] && [ "$DBNAME" = "gtdb_207" ]; then
            if Rscript -e "
                library(data.table)
                if(file.exists('$lineage_file')) {
                    df <- fread('$lineage_file', nrows = 1000)
                    if('sample_id' %in% names(df)) {
                        has_incorrect <- any(grepl('_gtdb_207_filtered$', df\$sample_id))
                        if(has_incorrect) {
                            cat('INCORRECT_IDS')
                        }
                    }
                }
            " 2>/dev/null | grep -q "INCORRECT_IDS"; then
                needs_merge=true
                echo "  Found incorrect sample_ids (ending with _gtdb_207_filtered), will regenerate"
            fi
        fi
        
        # Check if there are .b2 files with sample_ids not in the merged CSV
        # This catches cases where merged CSV exists but doesn't include all .b2 files
        if [ "$needs_merge" = false ] && [ -f "$lineage_file" ] && [ -n "$b2_files" ]; then
            # Get sample_ids from merged CSV
            merged_sample_count=$(Rscript -e "
                library(data.table)
                if(file.exists('$lineage_file')) {
                    df <- fread('$lineage_file', select = 'sample_id', showProgress = FALSE)
                    cat(length(unique(df\$sample_id)))
                } else {
                    cat('0')
                }
            " 2>/dev/null)
            
            # Count unique .b2 files (each .b2 file represents one sample)
            b2_file_count=$(echo "$b2_files" | wc -l | tr -d ' ')
            
            # If we have more .b2 files than unique samples in merged CSV, we need to merge
            if [ "$b2_file_count" -gt "${merged_sample_count:-0}" ]; then
                needs_merge=true
                echo "  Found $b2_file_count .b2 file(s) but merged CSV only has ${merged_sample_count} unique sample(s), will re-merge"
            fi
        fi
        
        if [ "$needs_merge" = false ]; then
            echo "  $lineage_file already exists and is up-to-date, skipping merge"
        fi
    fi
    
    # Merge bracken files if needed
    if [ "$needs_merge" = true ]; then
        # Check if existing merged file exists (from previous runs with old .b2 files)
        if [ -f "$merged_file" ]; then
            echo "  Existing merged file found, will combine with new .b2 files"
            
            # Identify new .b2 files (newer than existing merged file)
            # Suppress all output except the modification time number
            merged_file_time=$(stat -f %m "$merged_file" 2>&1 | grep -oE '^[0-9]+$' || stat -c %Y "$merged_file" 2>&1 | grep -oE '^[0-9]+$' || echo 0)
            new_b2_files=""
            old_b2_files=""
            
            for b2_file in $b2_files; do
                b2_file_time=$(stat -f %m "$b2_file" 2>&1 | grep -oE '^[0-9]+$' || stat -c %Y "$b2_file" 2>&1 | grep -oE '^[0-9]+$' || echo 0)
                if [ "$b2_file_time" -gt "$merged_file_time" ]; then
                    new_b2_files="$new_b2_files $b2_file"
                else
                    old_b2_files="$old_b2_files $b2_file"
                fi
            done
            
            # If there are new .b2 files, merge them and combine with existing merged file
            if [ -n "$new_b2_files" ]; then
                echo "  Found $(echo $new_b2_files | wc -w | tr -d ' ') new .b2 file(s) to merge"
                temp_merged_file="${merged_file}.tmp_new"
                
                # Merge only new .b2 files
                architeuthis merge -o "$temp_merged_file" $new_b2_files
                
                # Combine existing merged file with new merged file using R
                Rscript -e "
                library(data.table)
                old_merged <- fread('$merged_file')
                new_merged <- fread('$temp_merged_file')
                
                # Combine by summing reads for matching sample_id + name + taxonomy_id
                # Group by sample_id, name, taxonomy_id, taxonomy_lvl
                combined <- rbindlist(list(old_merged, new_merged))
                combined[, kraken_assigned_reads := sum(kraken_assigned_reads, na.rm=TRUE), 
                         by=.(sample_id, name, taxonomy_id, taxonomy_lvl)]
                combined[, added_reads := sum(added_reads, na.rm=TRUE), 
                         by=.(sample_id, name, taxonomy_id, taxonomy_lvl)]
                combined[, new_est_reads := sum(new_est_reads, na.rm=TRUE), 
                         by=.(sample_id, name, taxonomy_id, taxonomy_lvl)]
                
                # Remove duplicates and recalculate fraction_total_reads
                combined <- unique(combined, by=c('sample_id', 'name', 'taxonomy_id', 'taxonomy_lvl'))
                
                # Recalculate fraction_total_reads per sample
                combined[, total_reads_sample := sum(new_est_reads), by=sample_id]
                combined[, fraction_total_reads := new_est_reads / total_reads_sample]
                combined[, total_reads_sample := NULL]
                
                # Normalize sample_id to remove rank suffixes for consistency across ranks
                # Remove patterns like: _genus, _phylum, _domain, _class, _order, _family
                # But keep _filtered or __filtered
                combined[, sample_id := gsub('_(genus|phylum|domain|class|order|family)_filtered', '_filtered', sample_id)]
                combined[, sample_id := gsub('__(genus|phylum|domain|class|order|family)_filtered', '__filtered', sample_id)]
                
                fwrite(combined, '$merged_file', sep=',')
                "
                
                # Clean up temporary file
                rm -f "$temp_merged_file"
                
                echo "  ✓ Combined existing merged file with new .b2 files"
            else
                echo "  No new .b2 files found (all are older than merged file)"
            fi
        else
            # No existing merged file, merge all .b2 files normally
            architeuthis merge -o $merged_file $b2_files
            
            # Normalize sample_id immediately after initial merge
            Rscript -e "
            library(data.table)
            df <- fread('$merged_file')
            
            # Remove rank suffixes from sample_id to make them consistent across ranks
            # Remove patterns like: _genus, _phylum, _domain, _class, _order, _family
            # But keep _filtered or __filtered
            df[, sample_id := gsub('_(genus|phylum|domain|class|order|family)_filtered', '_filtered', sample_id)]
            df[, sample_id := gsub('__(genus|phylum|domain|class|order|family)_filtered', '__filtered', sample_id)]
            
            fwrite(df, '$merged_file', sep=',')
            "
        fi
        
        # Normalize sample_id to remove rank suffixes for consistency across ranks (if not already done)
        # This ensures species, genus, phylum, etc. all have the same sample_id format
        # Pattern: SAMPLE_db_rank_filtered -> SAMPLE_db_filtered
        # Pattern: SAMPLE_db__rank_filtered -> SAMPLE_db__filtered
        echo "  Normalizing sample_id to remove rank suffixes..."
        Rscript -e "
        library(data.table)
        df <- fread('$merged_file')
        
        # Remove rank suffixes from sample_id to make them consistent across ranks
        # Remove patterns like: _genus, _phylum, _domain, _class, _order, _family
        # But keep _filtered or __filtered
        df[, sample_id := gsub('_(genus|phylum|domain|class|order|family)_filtered', '_filtered', sample_id)]
        df[, sample_id := gsub('__(genus|phylum|domain|class|order|family)_filtered', '__filtered', sample_id)]
        
        fwrite(df, '$merged_file', sep=',')
        "
        
        # Run lineage on the (possibly combined) merged file
        architeuthis lineage $merged_file --db $DBDIR --data-dir $DB_taxonomy_dir --out $lineage_file
        
        if [ -f "$lineage_file" ]; then
            echo "  ✓ Successfully created/updated $lineage_file"
        else
            echo "  WARNING: $lineage_file not created"
        fi
    fi
}

# Define databases and their paths
declare -A DB_PATHS
declare -A DB_TAXONOMY_DIRS

# soil_microbe_db
DB_PATHS[soil_microbe_db]="/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2"
DB_TAXONOMY_DIRS[soil_microbe_db]="/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/"

# gtdb_207
DB_PATHS[gtdb_207]="/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2"
DB_TAXONOMY_DIRS[gtdb_207]="/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxonomy/"

# gtdb_207_unfiltered
DB_PATHS[gtdb_207_unfiltered]="/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2"
DB_TAXONOMY_DIRS[gtdb_207_unfiltered]="/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/taxonomy"

# pluspf
DB_PATHS[pluspf]="/projectnb/frpmars/soil_microbe_db/databases/pluspf"
DB_TAXONOMY_DIRS[pluspf]="/projectnb/microbiome/ref_db/NCBI-taxdump"

# Define taxonomic ranks to process
ranks=("species" "genus" "domain" "phylum")

# Process each database at each taxonomic rank
for DBNAME in "${!DB_PATHS[@]}"; do
    DBDIR=${DB_PATHS[$DBNAME]}
    DB_taxonomy_dir=${DB_TAXONOMY_DIRS[$DBNAME]}
    
    for rank in "${ranks[@]}"; do
        # Determine file pattern based on database and rank
        # Patterns must match actual file naming conventions
        case "$rank" in
            "species")
                if [ "$DBNAME" = "soil_microbe_db" ]; then
                    b2_pattern="*${DBNAME}_filtered.b2"
                elif [ "$DBNAME" = "gtdb_207_unfiltered" ]; then
                    b2_pattern="*gtdb_207_unfiltered*filtered.b2"
                elif [ "$DBNAME" = "gtdb_207" ]; then
                    b2_pattern="*gtdb_207__filtered.b2"
                elif [ "$DBNAME" = "pluspf" ]; then
                    b2_pattern="*pluspf__filtered.b2"
                else
                    b2_pattern="*${DBNAME}_filtered.b2"
                fi
                ;;
            "genus")
                if [ "$DBNAME" = "soil_microbe_db" ]; then
                    b2_pattern="*${DBNAME}_genus_filtered.b2"
                elif [ "$DBNAME" = "gtdb_207_unfiltered" ]; then
                    b2_pattern="*gtdb_207_unfiltered_genus_filtered.b2"
                elif [ "$DBNAME" = "gtdb_207" ]; then
                    b2_pattern="*gtdb_207__genus_filtered.b2"
                elif [ "$DBNAME" = "pluspf" ]; then
                    b2_pattern="*pluspf__genus_filtered.b2"
                else
                    b2_pattern="*${DBNAME}_genus_filtered.b2"
                fi
                ;;
            "domain")
                if [ "$DBNAME" = "soil_microbe_db" ]; then
                    b2_pattern="*${DBNAME}_domain_filtered.b2"
                elif [ "$DBNAME" = "gtdb_207_unfiltered" ]; then
                    b2_pattern="*gtdb_207_unfiltered_domain_filtered.b2"
                elif [ "$DBNAME" = "gtdb_207" ]; then
                    b2_pattern="*gtdb_207__domain_filtered.b2"
                elif [ "$DBNAME" = "pluspf" ]; then
                    b2_pattern="*pluspf__domain_filtered.b2"
                else
                    b2_pattern="*${DBNAME}_domain_filtered.b2"
                fi
                ;;
            "phylum")
                if [ "$DBNAME" = "soil_microbe_db" ]; then
                    b2_pattern="*${DBNAME}_phylum_filtered.b2"
                elif [ "$DBNAME" = "gtdb_207_unfiltered" ]; then
                    b2_pattern="*gtdb_207_unfiltered_phylum_filtered.b2"
                elif [ "$DBNAME" = "gtdb_207" ]; then
                    b2_pattern="*gtdb_207__phylum_filtered.b2"
                elif [ "$DBNAME" = "pluspf" ]; then
                    b2_pattern="*pluspf__phylum_filtered.b2"
                else
                    b2_pattern="*${DBNAME}_phylum_filtered.b2"
                fi
                ;;
        esac
        
        merge_and_add_lineage "$DBNAME" "$DBDIR" "$DB_taxonomy_dir" "$b2_pattern" "$rank"
    done
done

echo "All databases processed at species, genus, domain, and phylum levels."
echo ""
echo "Note: Cleanup of individual .b2 files is handled by scripts/run_workflow/05_cleanup_intermediate_files.sh"
echo "      This script will delete .b2 files only after verifying merged CSV files exist."
