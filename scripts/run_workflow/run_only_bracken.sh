#!/bin/bash
# Run Bracken at phylum level for specified database
# Then merge and add lineage
# Note: Cleanup of .b2 files is handled by scripts/run_workflow/05_cleanup_intermediate_files.sh

module load miniconda
conda activate struo2
cd /projectnb/frpmars/soil_microbe_db/

architeuthis_dir_path=/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/02_bracken_output
kraken_dir_path=/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/01_kraken_output
bracken_outdir_path=/projectnb/dietzelab/zrwerbin/soil_genome_db_output/bracken_phylum
summary_dir=/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/summary_files

# Set database - change this for different databases
DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/

time_with_seconds=$(date +%T)
echo "Beginning phylum-level Bracken for $DBNAME at: $time_with_seconds"

# Create output directory if it doesn't exist
mkdir -p $bracken_outdir_path

for samp_file in ${kraken_dir_path}/*.kreport; do
	samp_file_basename="$(basename -- $samp_file)"
	samp_name=${samp_file_basename::-15}
	
	# Only process files for the current database
	if [[ $samp_name != *$DBNAME* ]]; then
		continue
	fi
	
	echo "Processing $samp_name..."
	
	ARCHITEUTHIS_FILTERED_REPORT=${architeuthis_dir_path}/${samp_name}_filtered_kraken.kreport
	BRACKEN_OUTPUT_PHYLUM=${bracken_outdir_path}/${samp_name}_phylum_filtered.b2
	
	if test -e $BRACKEN_OUTPUT_PHYLUM; then
		echo "  $BRACKEN_OUTPUT_PHYLUM exists; skipping."
	else
		if test -e $ARCHITEUTHIS_FILTERED_REPORT; then
			/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_PHYLUM -l P
		else
			echo "  WARNING: $ARCHITEUTHIS_FILTERED_REPORT not found, skipping $samp_name"
		fi
	fi
done

time_with_seconds=$(date +%T)
echo "Finished phylum-level Bracken at: $time_with_seconds"

# Merge phylum files and add lineage
merged_file=${summary_dir}/${DBNAME}_filtered_phylum_merged.csv
lineage_merged_file=${summary_dir}/${DBNAME}_filtered_phylum_merged_lineage.csv

echo "Merging phylum files..."
if [ ! -f "$lineage_merged_file" ]; then
	architeuthis merge -o $merged_file ${bracken_outdir_path}/*${DBNAME}_phylum_filtered.b2
	architeuthis lineage $merged_file --db $DBDIR --data-dir $DB_taxonomy_dir --out $lineage_merged_file
	
	if [ -f "$lineage_merged_file" ]; then
		echo "✓ Successfully created $lineage_merged_file"
	else
		echo "WARNING: $lineage_merged_file not created"
	fi
else
	echo "$lineage_merged_file already exists, skipping merge"
fi

echo "Phylum analysis complete for $DBNAME"
echo ""
echo "Note: Cleanup of individual .b2 files is handled by scripts/run_workflow/05_cleanup_intermediate_files.sh"
echo "      This script will delete .b2 files only after verifying merged CSV files exist."
