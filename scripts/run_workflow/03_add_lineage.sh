#!/bin/bash
# After running Kraken, Architeuthis, and Bracken
# Merge Bracken outputs with "architeuthis merge"
# Then, lineage info can be added to the merged file with "architeuthis lineage"

module load miniconda
conda activate struo2
cd /projectnb/frpmars/soil_microbe_db/

bracken_dir=/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/02_bracken_output
summary_dir=/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files

# Function to merge and add lineage
merge_and_add_lineage() {
    local DBNAME=$1
    local DBDIR=$2
    local DB_taxonomy_dir=$3
    local b2_pattern=$4
    local rank=$5  # species, genus, domain, or phylum
    
    # Handle naming convention: phylum files use "filtered_phylum", others use rank name
    if [ "$rank" = "phylum" ]; then
        merged_file=${summary_dir}/${DBNAME}_filtered_phylum_merged.csv
        lineage_file=${summary_dir}/${DBNAME}_filtered_phylum_merged_lineage.csv
    else
        merged_file=${summary_dir}/${DBNAME}_${rank}_merged.csv
        lineage_file=${summary_dir}/${DBNAME}_${rank}_merged_lineage.csv
    fi
    
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
    
    # Merge bracken files
    if [ ! -f "$lineage_file" ]; then
        architeuthis merge -o $merged_file $b2_files
        architeuthis lineage $merged_file --db $DBDIR --data-dir $DB_taxonomy_dir --out $lineage_file
        
        if [ -f "$lineage_file" ]; then
            echo "  ✓ Successfully created $lineage_file"
        else
            echo "  WARNING: $lineage_file not created"
        fi
    else
        echo "  $lineage_file already exists, skipping merge"
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
