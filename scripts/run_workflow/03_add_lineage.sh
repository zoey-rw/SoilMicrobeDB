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
    b2_files=$(find ${bracken_dir} -name "${b2_pattern}" -type f 2>/dev/null | sort)
    
    if [ -z "$b2_files" ]; then
        echo "  WARNING: No files found matching pattern ${b2_pattern} in ${bracken_dir}"
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
        merged_file_time=$(stat -f %m "$lineage_file" 2>/dev/null || stat -c %Y "$lineage_file" 2>/dev/null || echo 0)
        
        for b2_file in $b2_files; do
            b2_file_time=$(stat -f %m "$b2_file" 2>/dev/null || stat -c %Y "$b2_file" 2>/dev/null || echo 0)
            if [ "$b2_file_time" -gt "$merged_file_time" ]; then
                needs_merge=true
                echo "  Found newer .b2 file(s), will re-merge"
                break
            fi
        done
        
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
            merged_file_time=$(stat -f %m "$merged_file" 2>/dev/null || stat -c %Y "$merged_file" 2>/dev/null || echo 0)
            new_b2_files=""
            old_b2_files=""
            
            for b2_file in $b2_files; do
                b2_file_time=$(stat -f %m "$b2_file" 2>/dev/null || stat -c %Y "$b2_file" 2>/dev/null || echo 0)
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
        fi
        
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
        # Determine file pattern based on rank
        case "$rank" in
            "species")
                b2_pattern="*${DBNAME}_filtered.b2"
                ;;
            "genus")
                b2_pattern="*${DBNAME}_genus_filtered.b2"
                ;;
            "domain")
                b2_pattern="*${DBNAME}_domain_filtered.b2"
                ;;
            "phylum")
                b2_pattern="*${DBNAME}_phylum_filtered.b2"
                ;;
        esac
        
        merge_and_add_lineage "$DBNAME" "$DBDIR" "$DB_taxonomy_dir" "$b2_pattern" "$rank"
    done
done

echo "All databases processed at species, genus, domain, and phylum levels."
echo ""
echo "Note: Cleanup of individual .b2 files is handled by scripts/run_workflow/05_cleanup_intermediate_files.sh"
echo "      This script will delete .b2 files only after verifying merged CSV files exist."
