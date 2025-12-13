#!/bin/bash
# After running Kraken, Architeuthis, and Bracken
# Merge Bracken outputs with "architeuthis merge"
# Then, lineage info can be added to the merged file with "architeuthis lineage"
# Note: Cleanup of .b2 files is handled by scripts/run_workflow/05_cleanup_intermediate_files.sh

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
    
    # Merge bracken files
    if [ ! -f "$lineage_file" ]; then
        architeuthis merge -o $merged_file ${bracken_dir}/${b2_pattern}
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

# Process each database at species level
DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/
merge_and_add_lineage $DBNAME $DBDIR $DB_taxonomy_dir "*${DBNAME}_filtered.b2" "species"

DBNAME=gtdb_207
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxonomy/
merge_and_add_lineage $DBNAME $DBDIR $DB_taxonomy_dir "*gtdb_207_filtered.b2" "species"

DBNAME=gtdb_207_unfiltered
DBDIR=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2
DB_taxonomy_dir=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/taxonomy
merge_and_add_lineage $DBNAME $DBDIR $DB_taxonomy_dir "*gtdb_207_unfiltered.b2" "species"

DBNAME=pluspf
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/pluspf
DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump
merge_and_add_lineage $DBNAME $DBDIR $DB_taxonomy_dir "*pluspf_filtered.b2" "species"

# Process each database at phylum level
DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/
merge_and_add_lineage $DBNAME $DBDIR $DB_taxonomy_dir "*${DBNAME}_phylum_filtered.b2" "phylum"

DBNAME=gtdb_207
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxonomy/
merge_and_add_lineage $DBNAME $DBDIR $DB_taxonomy_dir "*gtdb_207_phylum_filtered.b2" "phylum"

DBNAME=gtdb_207_unfiltered
DBDIR=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2
DB_taxonomy_dir=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/taxonomy
merge_and_add_lineage $DBNAME $DBDIR $DB_taxonomy_dir "*gtdb_207_unfiltered_phylum_filtered.b2" "phylum"

DBNAME=pluspf
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/pluspf
DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump
merge_and_add_lineage $DBNAME $DBDIR $DB_taxonomy_dir "*pluspf_phylum_filtered.b2" "phylum"

echo "All databases processed at species and phylum levels."
echo ""
echo "Note: Cleanup of individual .b2 files is handled by scripts/run_workflow/05_cleanup_intermediate_files.sh"
echo "      This script will delete .b2 files only after verifying merged CSV files exist."
